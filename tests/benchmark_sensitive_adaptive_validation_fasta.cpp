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
#include "primerpair/sensitive_adaptive_search.hpp"
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
                        static_cast<unsigned char>(raw)
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

std::unordered_set<std::string>
load_excluded_primers(
    const std::string& path
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open training TSV."
        );
    }

    std::unordered_set<std::string>
        excluded;

    std::string line;

    while (std::getline(input, line)) {

        if (
            line.empty() ||
            line.front() == '#'
        ) {
            continue;
        }

        const std::size_t first_tab =
            line.find('\t');

        if (
            first_tab ==
            std::string::npos
        ) {
            continue;
        }

        const std::size_t second_tab =
            line.find(
                '\t',
                first_tab + 1
            );

        if (
            second_tab ==
            std::string::npos
        ) {
            continue;
        }

        const std::string primer =
            line.substr(
                first_tab + 1,
                second_tab -
                first_tab -
                1
            );

        if (primer == "primer") {
            continue;
        }

        excluded.insert(
            primer
        );
    }

    if (excluded.empty()) {
        throw std::runtime_error(
            "No excluded primers loaded."
        );
    }

    return excluded;
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

    const double position =
        fraction *
        static_cast<double>(
            values.size() - 1
        );

    const std::size_t lo =
        static_cast<std::size_t>(
            position
        );

    const std::size_t hi =
        std::min(
            lo + 1,
            values.size() - 1
        );

    const double weight =
        position -
        static_cast<double>(
            lo
        );

    return
        values[lo] *
        (1.0 - weight)
        +
        values[hi] *
        weight;
}

double mean(
    const std::vector<double>& values
) {
    double total = 0.0;

    for (const double value : values) {
        total += value;
    }

    return
        values.empty()
            ? 0.0
            : total /
              static_cast<double>(
                  values.size()
              );
}

}  // namespace

int main(
    int argc,
    char* argv[]
) {
    try {
        if (
            argc < 3 ||
            argc > 5
        ) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " <training_router_tsv>"
                << " [primer_count]"
                << " [repeats]\n";

            return 2;
        }

        constexpr std::size_t
            primer_length = 20;

        constexpr std::size_t
            budget = 3;

        /*
         * Deliberately different from the
         * training benchmark stride.
         */
        constexpr std::uint64_t
            validation_stride =
                130363ULL;

        constexpr std::uint64_t
            initial_offset =
                7919ULL;

        const std::size_t primer_count =
            argc >= 4
                ? static_cast<std::size_t>(
                      std::stoull(argv[3])
                  )
                : 256;

        const std::size_t repeats =
            argc >= 5
                ? static_cast<std::size_t>(
                      std::stoull(argv[4])
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

        const auto excluded =
            load_excluded_primers(
                argv[2]
            );

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

        std::uint64_t cursor =
            initial_offset %
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
                    validation_stride
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
                excluded.contains(
                    primer
                )
            ) {
                continue;
            }

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
                "Could not collect requested "
                "independent validation primers."
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
            SensitiveAdaptiveSearchEngine
                adaptive_engine(
                    index,
                    packed
                );

        std::vector<double>
            exhaustive_medians;

        std::vector<double>
            adaptive_medians;

        exhaustive_medians.reserve(
            primers.size()
        );

        adaptive_medians.reserve(
            primers.size()
        );

        std::uint64_t
            candidate_routes = 0;

        std::uint64_t
            exhaustive_routes = 0;

        std::uint64_t
            equality_violations = 0;

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "# reference_bp\t"
            << reference.size()
            << '\n';

        std::cout
            << "# excluded_training_primers\t"
            << excluded.size()
            << '\n';

        std::cout
            << "# validation_primers\t"
            << primers.size()
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
            << "# threshold\t"
            << adaptive_engine
                .k3_max_seed_threshold()
            << '\n';

        std::cout
            << "primer_index"
            << '\t'
            << "primer"
            << '\t'
            << "backend"
            << '\t'
            << "max_seed_occurrences"
            << '\t'
            << "hit_count"
            << '\t'
            << "exhaustive_median_us"
            << '\t'
            << "adaptive_median_us"
            << '\t'
            << "adaptive_speedup"
            << '\t'
            << "all_checks"
            << '\n';

        for (
            std::size_t primer_index = 0;
            primer_index < primers.size();
            ++primer_index
        ) {
            const auto& primer =
                primers[primer_index];

            std::vector<double>
                exhaustive_times;

            std::vector<double>
                adaptive_times;

            exhaustive_times.reserve(
                repeats
            );

            adaptive_times.reserve(
                repeats
            );

            primerpair::
                SensitiveAdaptiveBackend
                    stable_backend =
                        primerpair::
                            SensitiveAdaptiveBackend::
                                Exhaustive;

            std::uint64_t
                stable_max_seed = 0;

            std::size_t
                stable_hit_count = 0;

            bool first = true;

            for (
                std::size_t repeat = 0;
                repeat < repeats;
                ++repeat
            ) {
                primerpair::
                    SensitivePrimerSearchResult
                        exhaustive_result;

                primerpair::
                    SensitiveAdaptiveSearchResult
                        adaptive_result;

                double exhaustive_us = 0.0;
                double adaptive_us = 0.0;

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

                auto run_adaptive =
                    [&]() {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        adaptive_result =
                            adaptive_engine.search(
                                primer,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        adaptive_us =
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

                if (
                    (
                        primer_index +
                        repeat
                    ) %
                    2 ==
                    0
                ) {
                    run_exhaustive();
                    run_adaptive();

                } else {

                    run_adaptive();
                    run_exhaustive();
                }

                if (
                    exhaustive_result.hits !=
                    adaptive_result
                        .search_result
                        .hits
                ) {
                    ++equality_violations;

                    throw std::logic_error(
                        "Adaptive result differs "
                        "from exhaustive reference."
                    );
                }

                if (first) {

                    stable_backend =
                        adaptive_result.backend;

                    stable_max_seed =
                        adaptive_result
                            .cost_estimate
                            .max_seed_occurrences;

                    stable_hit_count =
                        adaptive_result
                            .hit_count();

                    first = false;

                } else {

                    if (
                        adaptive_result.backend !=
                        stable_backend
                    ) {
                        throw std::logic_error(
                            "Adaptive backend is "
                            "not deterministic."
                        );
                    }

                    if (
                        adaptive_result
                            .cost_estimate
                            .max_seed_occurrences
                        !=
                        stable_max_seed
                    ) {
                        throw std::logic_error(
                            "Adaptive estimate is "
                            "not deterministic."
                        );
                    }
                }

                exhaustive_times.push_back(
                    exhaustive_us
                );

                adaptive_times.push_back(
                    adaptive_us
                );
            }

            const double ex_median =
                percentile(
                    exhaustive_times,
                    0.50
                );

            const double adaptive_median =
                percentile(
                    adaptive_times,
                    0.50
                );

            exhaustive_medians.push_back(
                ex_median
            );

            adaptive_medians.push_back(
                adaptive_median
            );

            if (
                stable_backend ==
                primerpair::
                    SensitiveAdaptiveBackend::
                        Candidate
            ) {
                ++candidate_routes;

            } else {

                ++exhaustive_routes;
            }

            const double speedup =
                adaptive_median > 0.0
                    ? ex_median /
                      adaptive_median
                    : 0.0;

            std::cout
                << primer_index
                << '\t'
                << primer
                << '\t'
                << (
                    stable_backend ==
                    primerpair::
                        SensitiveAdaptiveBackend::
                            Candidate
                        ? "CANDIDATE"
                        : "EXHAUSTIVE"
                )
                << '\t'
                << stable_max_seed
                << '\t'
                << stable_hit_count
                << '\t'
                << ex_median
                << '\t'
                << adaptive_median
                << '\t'
                << speedup
                << '\t'
                << "YES"
                << '\n';
        }

        const double exhaustive_mean =
            mean(
                exhaustive_medians
            );

        const double adaptive_mean =
            mean(
                adaptive_medians
            );

        const double exhaustive_median =
            percentile(
                exhaustive_medians,
                0.50
            );

        const double adaptive_median =
            percentile(
                adaptive_medians,
                0.50
            );

        const double exhaustive_p95 =
            percentile(
                exhaustive_medians,
                0.95
            );

        const double adaptive_p95 =
            percentile(
                adaptive_medians,
                0.95
            );

        const double mean_speedup =
            adaptive_mean > 0.0
                ? exhaustive_mean /
                  adaptive_mean
                : 0.0;

        std::cout
            << "# SUMMARY"
            << '\t'
            << "candidate_routes="
            << candidate_routes
            << '\t'
            << "exhaustive_routes="
            << exhaustive_routes
            << '\t'
            << "equality_violations="
            << equality_violations
            << '\n';

        std::cout
            << "# LATENCY"
            << '\t'
            << "exhaustive_mean_us="
            << exhaustive_mean
            << '\t'
            << "adaptive_mean_us="
            << adaptive_mean
            << '\t'
            << "mean_speedup="
            << mean_speedup
            << '\n';

        std::cout
            << "# LATENCY"
            << '\t'
            << "exhaustive_median_us="
            << exhaustive_median
            << '\t'
            << "adaptive_median_us="
            << adaptive_median
            << '\n';

        std::cout
            << "# LATENCY"
            << '\t'
            << "exhaustive_p95_us="
            << exhaustive_p95
            << '\t'
            << "adaptive_p95_us="
            << adaptive_p95
            << '\n';

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
