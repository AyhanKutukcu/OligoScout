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
#include <utility>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/sensitive_adaptive_search.hpp"
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
                    static_cast<unsigned char>(raw)
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
                    sequence.push_back(base);
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
                    sequence.push_back('N');
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


double median(
    std::vector<double> values
) {
    if (values.empty()) {
        return 0.0;
    }

    std::sort(
        values.begin(),
        values.end()
    );

    const std::size_t n =
        values.size();

    if (n % 2 == 1) {
        return values.at(
            n / 2
        );
    }

    return (
        values.at(n / 2 - 1)
        +
        values.at(n / 2)
    ) /
    2.0;
}


struct Sample {
    std::string primer;

    primerpair::
        SensitiveCandidateCostEstimate
            cost{};
};


std::vector<std::string>
collect_pool(
    const std::string& reference,
    const std::size_t primer_length,
    const std::size_t pool_size
) {
    if (
        reference.size() <
        primer_length
    ) {
        throw std::runtime_error(
            "Reference shorter than primer."
        );
    }

    const std::uint64_t available =
        static_cast<std::uint64_t>(
            reference.size() -
            primer_length +
            1
        );

    /*
     * Coprime-style deterministic stride.
     */
    constexpr std::uint64_t
        stride = 104729ULL;

    std::uint64_t cursor =
        (
            1000003ULL *
            static_cast<std::uint64_t>(
                primer_length
            )
        )
        %
        available;

    std::vector<std::string>
        primers;

    primers.reserve(
        pool_size
    );

    std::unordered_set<std::string>
        seen;

    seen.reserve(
        pool_size * 2
    );

    for (
        std::size_t attempt = 0;
        attempt < reference.size() &&
        primers.size() < pool_size;
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

        std::string primer(
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
            std::move(primer)
        );
    }

    if (
        primers.size() !=
        pool_size
    ) {
        throw std::runtime_error(
            "Could not collect requested primer pool."
        );
    }

    return primers;
}


template <class Function>
double timed_us(
    Function&& function
) {
    const auto start =
        std::chrono::
            steady_clock::now();

    function();

    const auto stop =
        std::chrono::
            steady_clock::now();

    return static_cast<double>(
        std::chrono::
            duration_cast<
                std::chrono::
                    nanoseconds
            >(
                stop - start
            ).count()
    ) /
    1000.0;
}

}  // namespace


int main(
    int argc,
    char* argv[]
) {
    try {
        if (
            argc < 2 ||
            argc > 6
        ) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " [pool_size=4096]"
                << " [top_n=16]"
                << " [repeats=3]"
                << " [sa_rate=8]\n";

            return 2;
        }

        const std::size_t pool_size =
            argc >= 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 4096;

        const std::size_t top_n =
            argc >= 4
                ? static_cast<std::size_t>(
                      std::stoull(argv[3])
                  )
                : 16;

        const std::size_t repeats =
            argc >= 5
                ? static_cast<std::size_t>(
                      std::stoull(argv[4])
                  )
                : 3;

        const std::size_t sa_rate =
            argc >= 6
                ? static_cast<std::size_t>(
                      std::stoull(argv[5])
                  )
                : 8;

        if (
            pool_size == 0 ||
            top_n == 0 ||
            top_n > pool_size ||
            repeats == 0 ||
            sa_rate == 0
        ) {
            throw std::invalid_argument(
                "Invalid benchmark parameters."
            );
        }

        constexpr std::size_t
            min_length = 23;

        constexpr std::size_t
            max_length = 35;

        constexpr std::size_t
            mismatch_budget = 3;


        const std::string reference =
            load_fasta(
                argv[1]
            );

        const primerpair::PackedReference
            packed(
                reference
            );

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
            SensitiveAdaptiveSearchEngine
                adaptive_engine(
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
            << "# pool_size\t"
            << pool_size
            << '\n';

        std::cout
            << "# top_n\t"
            << top_n
            << '\n';

        std::cout
            << "# repeats\t"
            << repeats
            << '\n';

        std::cout
            << "# sa_rate\t"
            << sa_rate
            << '\n';


        std::size_t
            total_equality_violations = 0;

        std::size_t
            total_candidate_slower = 0;


        for (
            std::size_t length = min_length;
            length <= max_length;
            ++length
        ) {
            const auto pool =
                collect_pool(
                    reference,
                    length,
                    pool_size
                );

            std::vector<Sample>
                samples;

            samples.reserve(
                pool.size()
            );

            /*
             * Estimator is used only to identify
             * deliberately difficult real genomic
             * primers for this stress benchmark.
             *
             * Production 23..35 routing does not
             * execute this estimator.
             */
            for (const auto& primer : pool) {

                samples.push_back(
                    Sample{
                        primer,
                        estimator.estimate_k3(
                            primer
                        )
                    }
                );
            }


            std::sort(
                samples.begin(),
                samples.end(),
                [](
                    const Sample& lhs,
                    const Sample& rhs
                ) {
                    if (
                        lhs.cost
                            .max_seed_occurrences
                        !=
                        rhs.cost
                            .max_seed_occurrences
                    ) {
                        return
                            lhs.cost
                                .max_seed_occurrences
                            >
                            rhs.cost
                                .max_seed_occurrences;
                    }

                    return
                        lhs.cost
                            .total_seed_occurrences
                        >
                        rhs.cost
                            .total_seed_occurrences;
                }
            );


            std::vector<double>
                exhaustive_times;

            std::vector<double>
                candidate_times;

            std::size_t
                equality_violations = 0;

            std::size_t
                candidate_slower = 0;

            double
                worst_candidate_over_exhaustive =
                    0.0;


            for (
                std::size_t rank = 0;
                rank < top_n;
                ++rank
            ) {
                const auto& sample =
                    samples.at(
                        rank
                    );

                std::vector<double>
                    exhaustive_repeat_times;

                std::vector<double>
                    candidate_repeat_times;

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
                    if (
                        (
                            rank +
                            repeat +
                            length
                        )
                        %
                        2
                        ==
                        0
                    ) {
                        exhaustive_repeat_times.push_back(
                            timed_us(
                                [&]() {
                                    exhaustive_result =
                                        exhaustive_engine.search(
                                            sample.primer,
                                            mismatch_budget
                                        );
                                }
                            )
                        );

                        candidate_repeat_times.push_back(
                            timed_us(
                                [&]() {
                                    candidate_result =
                                        candidate_engine.search(
                                            sample.primer,
                                            mismatch_budget
                                        );
                                }
                            )
                        );

                    } else {

                        candidate_repeat_times.push_back(
                            timed_us(
                                [&]() {
                                    candidate_result =
                                        candidate_engine.search(
                                            sample.primer,
                                            mismatch_budget
                                        );
                                }
                            )
                        );

                        exhaustive_repeat_times.push_back(
                            timed_us(
                                [&]() {
                                    exhaustive_result =
                                        exhaustive_engine.search(
                                            sample.primer,
                                            mismatch_budget
                                        );
                                }
                            )
                        );
                    }
                }


                const double exhaustive_us =
                    median(
                        exhaustive_repeat_times
                    );

                const double candidate_us =
                    median(
                        candidate_repeat_times
                    );

                exhaustive_times.push_back(
                    exhaustive_us
                );

                candidate_times.push_back(
                    candidate_us
                );


                const bool equality =
                    exhaustive_result.hits ==
                    candidate_result.hits;

                if (!equality) {
                    ++equality_violations;
                    ++total_equality_violations;
                }


                const auto adaptive =
                    adaptive_engine.search(
                        sample.primer,
                        mismatch_budget
                    );

                const bool adaptive_direct_candidate =
                    (
                        adaptive.backend ==
                        primerpair::
                            SensitiveAdaptiveBackend::
                                Candidate
                    )
                    &&
                    !adaptive.estimator_used;

                const bool adaptive_equality =
                    adaptive.search_result.hits ==
                    exhaustive_result.hits;

                if (!adaptive_equality) {
                    ++equality_violations;
                    ++total_equality_violations;
                }


                const double ratio =
                    exhaustive_us > 0.0
                        ? candidate_us /
                          exhaustive_us
                        : 0.0;

                worst_candidate_over_exhaustive =
                    std::max(
                        worst_candidate_over_exhaustive,
                        ratio
                    );

                const bool slower =
                    candidate_us >
                    exhaustive_us;

                if (slower) {
                    ++candidate_slower;
                    ++total_candidate_slower;
                }


                std::cout
                    << "STRESS\t"
                    << length
                    << '\t'
                    << rank
                    << '\t'
                    << sample.cost
                        .max_seed_occurrences
                    << '\t'
                    << sample.cost
                        .total_seed_occurrences
                    << '\t'
                    << exhaustive_result.hit_count()
                    << '\t'
                    << exhaustive_us
                    << '\t'
                    << candidate_us
                    << '\t'
                    << ratio
                    << '\t'
                    << (
                        slower
                            ? "EXHAUSTIVE"
                            : "CANDIDATE"
                    )
                    << '\t'
                    << (
                        equality
                            ? "YES"
                            : "NO"
                    )
                    << '\t'
                    << (
                        adaptive_direct_candidate
                            ? "YES"
                            : "NO"
                    )
                    << '\t'
                    << (
                        adaptive_equality
                            ? "YES"
                            : "NO"
                    )
                    << '\n';
            }


            const double exhaustive_mean =
                mean(
                    exhaustive_times
                );

            const double candidate_mean =
                mean(
                    candidate_times
                );

            const double speedup =
                candidate_mean > 0.0
                    ? exhaustive_mean /
                      candidate_mean
                    : 0.0;


            std::cout
                << "LENGTH_STRESS_SUMMARY\t"
                << length
                << '\t'
                << "selected="
                << top_n
                << '\t'
                << "max_seed_top="
                << samples.front()
                    .cost
                    .max_seed_occurrences
                << '\t'
                << "candidate_slower="
                << candidate_slower
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
                << "speedup="
                << speedup
                << '\t'
                << "worst_candidate_over_exhaustive="
                << worst_candidate_over_exhaustive
                << '\t'
                << "all_checks="
                << (
                    equality_violations == 0
                        ? "YES"
                        : "NO"
                )
                << '\n';
        }


        std::cout
            << "TOTAL_EQUALITY_VIOLATIONS\t"
            << total_equality_violations
            << '\n';

        std::cout
            << "TOTAL_CANDIDATE_SLOWER\t"
            << total_candidate_slower
            << '\n';

        std::cout
            << "ALL_CHECKS\t"
            << (
                total_equality_violations == 0
                    ? "YES"
                    : "NO"
            )
            << '\n';


        return
            total_equality_violations == 0
                ? 0
                : 1;

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
