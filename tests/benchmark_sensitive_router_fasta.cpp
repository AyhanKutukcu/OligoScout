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

    return
        (
            values.at(n / 2 - 1) +
            values.at(n / 2)
        ) /
        2.0;
}

}  // namespace

int main(
    int argc,
    char* argv[]
) {
    try {
        if (argc < 2 || argc > 4) {

            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " [primer_count]"
                << " [repeats]\n";

            return 2;
        }

        constexpr std::size_t
            primer_length = 20;

        constexpr std::size_t
            budget = 3;

        constexpr std::uint64_t
            stride = 104729ULL;

        const std::size_t primer_count =
            argc >= 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 256;

        const std::size_t repeats =
            argc >= 4
                ? static_cast<std::size_t>(
                      std::stoull(argv[3])
                  )
                : 3;

        if (
            primer_count == 0 ||
            repeats == 0
        ) {
            throw std::invalid_argument(
                "primer_count and repeats "
                "must be > 0."
            );
        }

        const std::string reference =
            load_fasta(
                argv[1]
            );

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

        std::uint64_t cursor = 0;

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
                ) %
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

        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
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
            << "# reference_bp\t"
            << reference.size()
            << '\n';

        std::cout
            << "# primer_count\t"
            << primers.size()
            << '\n';

        std::cout
            << "# primer_length\t"
            << primer_length
            << '\n';

        std::cout
            << "# budget\t"
            << budget
            << '\n';

        std::cout
            << "# repeats\t"
            << repeats
            << '\n';

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "primer_index"
            << '\t'
            << "primer"
            << '\t'
            << "forward_left_occurrences"
            << '\t'
            << "forward_right_occurrences"
            << '\t'
            << "reverse_left_occurrences"
            << '\t'
            << "reverse_right_occurrences"
            << '\t'
            << "total_seed_occurrences"
            << '\t'
            << "max_seed_occurrences"
            << '\t'
            << "estimator_median_us"
            << '\t'
            << "exhaustive_median_us"
            << '\t'
            << "candidate_median_us"
            << '\t'
            << "candidate_speedup"
            << '\t'
            << "winner"
            << '\t'
            << "hit_count"
            << '\t'
            << "all_checks"
            << '\n';

        for (
            std::size_t primer_index = 0;
            primer_index < primers.size();
            ++primer_index
        ) {
            const auto& primer =
                primers.at(
                    primer_index
                );

            std::vector<double>
                estimator_times;

            std::vector<double>
                exhaustive_times;

            std::vector<double>
                candidate_times;

            estimator_times.reserve(
                repeats
            );

            exhaustive_times.reserve(
                repeats
            );

            candidate_times.reserve(
                repeats
            );

            primerpair::
                SensitiveCandidateCostEstimate
                    stable_estimate{};

            bool have_estimate = false;

            std::size_t hit_count = 0;

            for (
                std::size_t repeat = 0;
                repeat < repeats;
                ++repeat
            ) {
                primerpair::
                    SensitivePrimerSearchResult
                        exhaustive_result;

                primerpair::
                    SensitivePrimerSearchResult
                        candidate_result;

                primerpair::
                    SensitiveCandidateCostEstimate
                        current_estimate;

                double estimator_us = 0.0;
                double exhaustive_us = 0.0;
                double candidate_us = 0.0;

                auto run_estimator =
                    [&]() {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        current_estimate =
                            estimator.estimate_k3(
                                primer
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        estimator_us =
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0;
                    };

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

                        exhaustive_us =
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0;
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

                        candidate_us =
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0;
                    };

                /*
                 * Rotate execution order to reduce
                 * systematic cache-order bias.
                 */
                switch (
                    (
                        primer_index +
                        repeat
                    ) %
                    3
                ) {
                    case 0:
                        run_estimator();
                        run_exhaustive();
                        run_candidate();
                        break;

                    case 1:
                        run_candidate();
                        run_estimator();
                        run_exhaustive();
                        break;

                    default:
                        run_exhaustive();
                        run_candidate();
                        run_estimator();
                        break;
                }

                if (
                    exhaustive_result.hits !=
                    candidate_result.hits
                ) {
                    throw std::logic_error(
                        "Candidate differs from "
                        "exhaustive reference."
                    );
                }

                if (!have_estimate) {

                    stable_estimate =
                        current_estimate;

                    have_estimate = true;

                } else if (
                    stable_estimate !=
                    current_estimate
                ) {
                    throw std::logic_error(
                        "Cost estimate is not deterministic."
                    );
                }

                hit_count =
                    exhaustive_result.hit_count();

                estimator_times.push_back(
                    estimator_us
                );

                exhaustive_times.push_back(
                    exhaustive_us
                );

                candidate_times.push_back(
                    candidate_us
                );
            }

            const double estimator_median =
                median(
                    estimator_times
                );

            const double exhaustive_median =
                median(
                    exhaustive_times
                );

            const double candidate_median =
                median(
                    candidate_times
                );

            const double candidate_speedup =
                candidate_median > 0.0
                    ? exhaustive_median /
                      candidate_median
                    : 0.0;

            const bool candidate_wins =
                candidate_median <
                exhaustive_median;

            std::cout
                << primer_index
                << '\t'
                << primer
                << '\t'
                << stable_estimate
                    .forward_left_occurrences
                << '\t'
                << stable_estimate
                    .forward_right_occurrences
                << '\t'
                << stable_estimate
                    .reverse_left_occurrences
                << '\t'
                << stable_estimate
                    .reverse_right_occurrences
                << '\t'
                << stable_estimate
                    .total_seed_occurrences
                << '\t'
                << stable_estimate
                    .max_seed_occurrences
                << '\t'
                << estimator_median
                << '\t'
                << exhaustive_median
                << '\t'
                << candidate_median
                << '\t'
                << candidate_speedup
                << '\t'
                << (
                    candidate_wins
                        ? "CANDIDATE"
                        : "EXHAUSTIVE"
                )
                << '\t'
                << hit_count
                << '\t'
                << "YES"
                << '\n';
        }

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
