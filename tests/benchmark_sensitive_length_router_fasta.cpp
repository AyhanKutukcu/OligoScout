#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/sensitive_candidate_search.hpp"
#include "primerpair/sensitive_cost_estimator.hpp"
#include "primerpair/sensitive_primer_search.hpp"

namespace {

std::string load_fasta(
    const std::string& path
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open FASTA."
        );
    }

    std::string sequence;
    std::string line;

    while (std::getline(input, line)) {

        if (
            line.empty() ||
            line.front() == '>'
        ) {
            continue;
        }

        for (const char raw : line) {

            if (
                std::isspace(
                    static_cast<unsigned char>(
                        raw
                    )
                )
            ) {
                continue;
            }

            const char base =
                static_cast<char>(
                    std::toupper(
                        static_cast<unsigned char>(
                            raw
                        )
                    )
                );

            switch (base) {

                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':
                    sequence.push_back(
                        base
                    );
                    break;

                case 'R':
                case 'Y':
                case 'S':
                case 'W':
                case 'K':
                case 'M':
                case 'B':
                case 'D':
                case 'H':
                case 'V':
                    sequence.push_back(
                        'N'
                    );
                    break;

                default:
                    throw std::runtime_error(
                        "Unsupported FASTA nucleotide."
                    );
            }
        }
    }

    if (sequence.empty()) {
        throw std::runtime_error(
            "Empty FASTA."
        );
    }

    return sequence;
}


bool all_acgt(
    const std::string_view sequence
) {
    return std::all_of(
        sequence.begin(),
        sequence.end(),
        [](
            const char base
        ) {
            return
                base == 'A' ||
                base == 'C' ||
                base == 'G' ||
                base == 'T';
        }
    );
}


double percentile(
    std::vector<double> values,
    const double fraction
) {
    if (values.empty()) {
        return 0.0;
    }

    std::sort(
        values.begin(),
        values.end()
    );

    const double scaled =
        fraction *
        static_cast<double>(
            values.size() - 1
        );

    const std::size_t lower =
        static_cast<std::size_t>(
            scaled
        );

    const std::size_t upper =
        std::min(
            lower + 1,
            values.size() - 1
        );

    const double weight =
        scaled -
        static_cast<double>(
            lower
        );

    return
        values.at(lower) *
            (1.0 - weight)
        +
        values.at(upper) *
            weight;
}


double mean(
    const std::vector<double>& values
) {
    if (values.empty()) {
        return 0.0;
    }

    long double total = 0.0L;

    for (const double value : values) {
        total += value;
    }

    return static_cast<double>(
        total /
        static_cast<long double>(
            values.size()
        )
    );
}


std::vector<std::string>
collect_primers(
    const std::string& reference,
    const std::size_t primer_length,
    const std::size_t primer_count
) {
    if (
        reference.size() <
        primer_length
    ) {
        throw std::runtime_error(
            "Reference shorter than primer."
        );
    }

    constexpr std::uint64_t
        stride = 104729ULL;

    const std::uint64_t available =
        static_cast<std::uint64_t>(
            reference.size() -
            primer_length +
            1
        );

    std::vector<std::string>
        primers;

    primers.reserve(
        primer_count
    );

    std::unordered_set<std::string>
        seen;

    seen.reserve(
        primer_count * 2
    );

    /*
     * Different deterministic starting point for
     * every primer length, while keeping the run
     * completely reproducible.
     */
    std::uint64_t cursor =
        (
            static_cast<std::uint64_t>(
                primer_length
            )
            *
            1000003ULL
        )
        %
        available;

    for (
        std::size_t attempt = 0;
        attempt < reference.size() &&
        primers.size() < primer_count;
        ++attempt
    ) {
        cursor =
            (
                cursor +
                stride
            )
            %
            available;

        const std::string_view candidate(
            reference.data() + cursor,
            primer_length
        );

        if (!all_acgt(candidate)) {
            continue;
        }

        const std::string primer(
            candidate
        );

        if (
            !seen.insert(
                primer
            ).second
        ) {
            continue;
        }

        primers.push_back(
            primer
        );
    }

    if (
        primers.size() !=
        primer_count
    ) {
        throw std::runtime_error(
            "Could not collect requested primers."
        );
    }

    return primers;
}

}  // namespace


int main(
    int argc,
    char* argv[]
) {
    try {
        if (
            argc < 2 ||
            argc > 5
        ) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " [primer_count=256]"
                << " [repeats=3]"
                << " [sa_rate=8]\n";

            return 2;
        }

        constexpr std::size_t
            min_primer_length = 18;

        constexpr std::size_t
            max_primer_length = 35;

        constexpr std::size_t
            budget = 3;

        const std::size_t primer_count =
            argc >= 3
                ?
                static_cast<std::size_t>(
                    std::stoull(
                        argv[2]
                    )
                )
                :
                256;

        const std::size_t repeats =
            argc >= 4
                ?
                static_cast<std::size_t>(
                    std::stoull(
                        argv[3]
                    )
                )
                :
                3;

        const std::size_t sa_rate =
            argc >= 5
                ?
                static_cast<std::size_t>(
                    std::stoull(
                        argv[4]
                    )
                )
                :
                8;

        if (
            primer_count == 0 ||
            repeats == 0 ||
            sa_rate == 0
        ) {
            throw std::invalid_argument(
                "primer_count, repeats and "
                "sa_rate must be > 0."
            );
        }

        const std::string reference =
            load_fasta(
                argv[1]
            );

        const primerpair::PackedReference
            packed(
                reference
            );

        /*
         * One index is reused for every primer
         * length. SA8 is the current balanced
         * training profile unless overridden.
         */
        const primerpair::BidirectionalFMIndex
            index(
                reference,
                sa_rate
            );

        const primerpair::
            SensitivePrimerSearchEngine
                exhaustive_engine(
                    index
                );

        const primerpair::
            SensitiveCandidateSearchEngine
                candidate_engine(
                    index,
                    packed
                );

        const primerpair::
            SensitiveCandidateCostEstimator
                estimator(
                    index
                );


        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "# reference_bp\t"
            << reference.size()
            << '\n';

        std::cout
            << "# primer_length_range\t"
            << min_primer_length
            << '-'
            << max_primer_length
            << '\n';

        std::cout
            << "# primer_count_per_length\t"
            << primer_count
            << '\n';

        std::cout
            << "# repeats\t"
            << repeats
            << '\n';

        std::cout
            << "# budget\t"
            << budget
            << '\n';

        std::cout
            << "# sa_rate\t"
            << sa_rate
            << '\n';

        std::cout
            << "# PRIMER_COLUMNS\t"
            << "length"
            << "\tprimer_index"
            << "\tmax_seed_occurrences"
            << "\ttotal_seed_occurrences"
            << "\testimator_us"
            << "\texhaustive_us"
            << "\tcandidate_us"
            << "\twinner"
            << "\texhaustive_hits"
            << "\tcandidate_hits"
            << "\tequality"
            << '\n';


        std::size_t
            total_equality_violations = 0;


        for (
            std::size_t primer_length =
                min_primer_length;
            primer_length <=
                max_primer_length;
            ++primer_length
        ) {
            const auto primers =
                collect_primers(
                    reference,
                    primer_length,
                    primer_count
                );

            std::vector<double>
                exhaustive_medians;

            std::vector<double>
                candidate_medians;

            exhaustive_medians.reserve(
                primers.size()
            );

            candidate_medians.reserve(
                primers.size()
            );

            std::uint64_t
                exhaustive_hits_total = 0;

            std::uint64_t
                candidate_hits_total = 0;

            std::size_t
                equality_violations = 0;

            std::size_t
                candidate_wins = 0;

            std::size_t
                exhaustive_wins = 0;


            for (
                std::size_t primer_index = 0;
                primer_index < primers.size();
                ++primer_index
            ) {
                const auto& primer =
                    primers.at(
                        primer_index
                    );

                primerpair::
                    SensitiveCandidateCostEstimate
                        cost{};

                std::vector<double>
                    estimator_times;

                estimator_times.reserve(
                    repeats
                );

                for (
                    std::size_t repeat = 0;
                    repeat < repeats;
                    ++repeat
                ) {
                    const auto start =
                        std::chrono::
                            steady_clock::now();

                    cost =
                        estimator.estimate_k3(
                            primer
                        );

                    const auto stop =
                        std::chrono::
                            steady_clock::now();

                    estimator_times.push_back(
                        static_cast<double>(
                            std::chrono::
                                duration_cast<
                                    std::chrono::
                                        nanoseconds
                                >(
                                    stop - start
                                ).count()
                        )
                        /
                        1000.0
                    );
                }

                const double estimator_us =
                    percentile(
                        estimator_times,
                        0.50
                    );

                std::vector<double>
                    exhaustive_times;

                std::vector<double>
                    candidate_times;

                exhaustive_times.reserve(
                    repeats
                );

                candidate_times.reserve(
                    repeats
                );

                primerpair::
                    SensitivePrimerSearchResult
                        exhaustive_result;

                primerpair::
                    SensitivePrimerSearchResult
                        candidate_result;


                for (
                    std::size_t repeat = 0;
                    repeat < repeats;
                    ++repeat
                ) {
                    auto run_exhaustive =
                        [&]() {
                            const auto start =
                                std::chrono::
                                    steady_clock::now();

                            exhaustive_result =
                                exhaustive_engine.search(
                                    primer,
                                    budget
                                );

                            const auto stop =
                                std::chrono::
                                    steady_clock::now();

                            exhaustive_times.push_back(
                                static_cast<double>(
                                    std::chrono::
                                        duration_cast<
                                            std::chrono::
                                                nanoseconds
                                        >(
                                            stop - start
                                        ).count()
                                )
                                /
                                1000.0
                            );
                        };


                    auto run_candidate =
                        [&]() {
                            const auto start =
                                std::chrono::
                                    steady_clock::now();

                            candidate_result =
                                candidate_engine.search(
                                    primer,
                                    budget
                                );

                            const auto stop =
                                std::chrono::
                                    steady_clock::now();

                            candidate_times.push_back(
                                static_cast<double>(
                                    std::chrono::
                                        duration_cast<
                                            std::chrono::
                                                nanoseconds
                                        >(
                                            stop - start
                                        ).count()
                                )
                                /
                                1000.0
                            );
                        };


                    /*
                     * Alternate timing order.
                     */
                    if (
                        (
                            primer_index +
                            repeat +
                            primer_length
                        )
                        %
                        2
                        ==
                        0
                    ) {
                        run_exhaustive();
                        run_candidate();

                    } else {

                        run_candidate();
                        run_exhaustive();
                    }
                }


                const double exhaustive_us =
                    percentile(
                        exhaustive_times,
                        0.50
                    );

                const double candidate_us =
                    percentile(
                        candidate_times,
                        0.50
                    );

                exhaustive_medians.push_back(
                    exhaustive_us
                );

                candidate_medians.push_back(
                    candidate_us
                );

                exhaustive_hits_total +=
                    exhaustive_result
                        .hit_count();

                candidate_hits_total +=
                    candidate_result
                        .hit_count();

                const bool equality =
                    candidate_result.hits ==
                    exhaustive_result.hits;

                if (!equality) {
                    ++equality_violations;
                    ++total_equality_violations;
                }

                const bool candidate_faster =
                    candidate_us <
                    exhaustive_us;

                if (candidate_faster) {
                    ++candidate_wins;

                } else {

                    ++exhaustive_wins;
                }


                std::cout
                    << "PRIMER\t"
                    << primer_length
                    << '\t'
                    << primer_index
                    << '\t'
                    << cost.max_seed_occurrences
                    << '\t'
                    << cost.total_seed_occurrences
                    << '\t'
                    << estimator_us
                    << '\t'
                    << exhaustive_us
                    << '\t'
                    << candidate_us
                    << '\t'
                    << (
                        candidate_faster
                            ? "CANDIDATE"
                            : "EXHAUSTIVE"
                    )
                    << '\t'
                    << exhaustive_result.hit_count()
                    << '\t'
                    << candidate_result.hit_count()
                    << '\t'
                    << (
                        equality
                            ? "YES"
                            : "NO"
                    )
                    << '\n';
            }


            const double exhaustive_mean =
                mean(
                    exhaustive_medians
                );

            const double candidate_mean =
                mean(
                    candidate_medians
                );

            const double exhaustive_median =
                percentile(
                    exhaustive_medians,
                    0.50
                );

            const double candidate_median =
                percentile(
                    candidate_medians,
                    0.50
                );

            const double exhaustive_p95 =
                percentile(
                    exhaustive_medians,
                    0.95
                );

            const double candidate_p95 =
                percentile(
                    candidate_medians,
                    0.95
                );

            const double mean_speedup =
                candidate_mean > 0.0
                    ?
                    exhaustive_mean /
                    candidate_mean
                    :
                    0.0;

            const bool all_checks =
                equality_violations == 0
                &&
                exhaustive_hits_total ==
                    candidate_hits_total;


            std::cout
                << "LENGTH_SUMMARY\t"
                << primer_length
                << '\t'
                << "left_seed="
                << (
                    primer_length /
                    2
                )
                << '\t'
                << "right_seed="
                << (
                    primer_length -
                    primer_length /
                    2
                )
                << '\t'
                << "candidate_wins="
                << candidate_wins
                << '\t'
                << "exhaustive_wins="
                << exhaustive_wins
                << '\t'
                << "equality_violations="
                << equality_violations
                << '\t'
                << "exhaustive_mean_us="
                << exhaustive_mean
                << '\t'
                << "candidate_mean_us="
                << candidate_mean
                << '\t'
                << "mean_speedup="
                << mean_speedup
                << '\t'
                << "exhaustive_median_us="
                << exhaustive_median
                << '\t'
                << "candidate_median_us="
                << candidate_median
                << '\t'
                << "exhaustive_p95_us="
                << exhaustive_p95
                << '\t'
                << "candidate_p95_us="
                << candidate_p95
                << '\t'
                << "all_checks="
                << (
                    all_checks
                        ? "YES"
                        : "NO"
                )
                << '\n';


            if (!all_checks) {
                throw std::logic_error(
                    "Length-aware Candidate search "
                    "differs from exhaustive reference."
                );
            }
        }


        std::cout
            << "TOTAL_EQUALITY_VIOLATIONS\t"
            << total_equality_violations
            << '\n';

        std::cout
            << "ALL_CHECKS\t"
            << (
                total_equality_violations == 0
                    ? "YES"
                    : "NO"
            )
            << '\n';

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
