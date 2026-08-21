#include <algorithm>
#include <array>
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

#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/sensitive_adaptive_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace {

struct DesignedPair {
    std::string primer1;
    std::string primer2;

    std::uint64_t left_position{0};
    std::uint64_t right_position{0};
    std::uint64_t amplicon_length{0};
};


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
                    sequence.push_back(base);
                    break;

                case 'N':
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

    const std::size_t lower =
        static_cast<std::size_t>(
            position
        );

    const std::size_t upper =
        std::min(
            lower + 1,
            values.size() - 1
        );

    const double weight =
        position -
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


double positive_difference(
    const double total,
    const double component
) {
    return
        total > component
            ? total - component
            : 0.0;
}

}  // namespace


int main(
    int argc,
    char* argv[]
) {
    try {
        if (
            argc < 2 ||
            argc > 4
        ) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " [pair_count]"
                << " [repeats]\n";

            return 2;
        }

        constexpr std::size_t
            primer_length = 20;

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::uint64_t
            min_amplicon = 50;

        constexpr std::uint64_t
            max_amplicon = 3000;

        /*
         * Same deterministic sampler as the
         * STRICT-vs-SENSITIVE pair benchmark.
         */
        constexpr std::uint64_t
            sampling_stride = 130363ULL;

        constexpr std::uint64_t
            initial_offset = 50021ULL;

        const std::size_t target_pairs =
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
            target_pairs == 0 ||
            repeats == 0
        ) {
            throw std::invalid_argument(
                "pair_count and repeats must be > 0."
            );
        }

        const std::string reference =
            load_fasta(
                argv[1]
            );

        if (
            reference.size() <=
                max_amplicon +
                primer_length
        ) {
            throw std::runtime_error(
                "Reference is too short."
            );
        }

        /*
         * --------------------------------------------------
         * Same deterministic real-reference pair set.
         * --------------------------------------------------
         */

        std::vector<DesignedPair>
            pairs;

        pairs.reserve(
            target_pairs
        );

        std::unordered_set<std::string>
            seen_pairs;

        seen_pairs.reserve(
            target_pairs * 2
        );

        const std::uint64_t scan_span =
            static_cast<std::uint64_t>(
                reference.size()
            ) -
            max_amplicon -
            1;

        std::uint64_t cursor =
            initial_offset %
            scan_span;

        const std::size_t max_attempts =
            target_pairs *
            100000;

        for (
            std::size_t attempt = 0;
            attempt < max_attempts &&
            pairs.size() < target_pairs;
            ++attempt
        ) {
            cursor =
                (
                    cursor +
                    sampling_stride
                )
                %
                scan_span;

            const std::uint64_t
                amplicon_length =
                    min_amplicon +
                    (
                        (
                            pairs.size() * 97ULL +
                            attempt * 13ULL
                        )
                        %
                        (
                            max_amplicon -
                            min_amplicon +
                            1
                        )
                    );

            const std::uint64_t
                right_position =
                    cursor +
                    amplicon_length -
                    primer_length;

            const std::string_view
                left_site(
                    reference.data() +
                        cursor,
                    primer_length
                );

            const std::string_view
                right_site(
                    reference.data() +
                        right_position,
                    primer_length
                );

            if (
                !all_acgt(
                    left_site
                )
                ||
                !all_acgt(
                    right_site
                )
            ) {
                continue;
            }

            const std::string primer1(
                left_site
            );

            const std::string primer2 =
                primerpair::
                    reverse_complement(
                        right_site
                    );

            if (primer1 == primer2) {
                continue;
            }

            const std::string key =
                primer1 +
                "\t" +
                primer2;

            if (
                !seen_pairs
                    .insert(key)
                    .second
            ) {
                continue;
            }

            pairs.push_back(
                DesignedPair{
                    primer1,
                    primer2,
                    cursor,
                    right_position,
                    amplicon_length
                }
            );
        }

        if (
            pairs.size() !=
            target_pairs
        ) {
            throw std::runtime_error(
                "Could not collect requested "
                "number of primer pairs."
            );
        }

        /*
         * --------------------------------------------------
         * Engines
         * --------------------------------------------------
         */

        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::
            StrandAwarePrimerSearchEngine
                strict_single_engine(
                    index,
                    packed
                );

        const primerpair::
            SensitiveAdaptiveSearchEngine
                sensitive_single_engine(
                    index,
                    packed
                );

        const primerpair::
            PrimerPairSearchEngine
                pair_engine(
                    index,
                    packed
                );

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "# reference_bp\t"
            << reference.size()
            << '\n';

        std::cout
            << "# designed_pairs\t"
            << pairs.size()
            << '\n';

        std::cout
            << "# primer_length\t"
            << primer_length
            << '\n';

        std::cout
            << "# anchor_length\t"
            << anchor_length
            << '\n';

        std::cout
            << "# amplicon_range\t"
            << min_amplicon
            << '-'
            << max_amplicon
            << '\n';

        std::cout
            << "# repeats\t"
            << repeats
            << '\n';

        std::cout
            << "budget"
            << '\t'
            << "strict_two_single_mean_us"
            << '\t'
            << "strict_full_pair_mean_us"
            << '\t'
            << "strict_est_nonsearch_us"
            << '\t'
            << "strict_search_share_pct"
            << '\t'
            << "sensitive_two_single_mean_us"
            << '\t'
            << "sensitive_full_pair_mean_us"
            << '\t'
            << "sensitive_est_nonsearch_us"
            << '\t'
            << "sensitive_search_share_pct"
            << '\t'
            << "sensitive_vs_strict_single_ratio"
            << '\t'
            << "sensitive_vs_strict_pair_ratio"
            << '\t'
            << "candidate_routes"
            << '\t'
            << "exhaustive_routes"
            << '\t'
            << "estimator_used_routes"
            << '\t'
            << "strict_amplicons"
            << '\t'
            << "sensitive_amplicons"
            << '\t'
            << "hit_count_violations"
            << '\t'
            << "strict_single_median_us"
            << '\t'
            << "strict_pair_median_us"
            << '\t'
            << "sensitive_single_median_us"
            << '\t'
            << "sensitive_pair_median_us"
            << '\t'
            << "strict_pair_p95_us"
            << '\t'
            << "sensitive_pair_p95_us"
            << '\t'
            << "all_checks"
            << '\n';

        constexpr std::array<
            std::size_t,
            4
        > budgets{
            0,
            1,
            2,
            3
        };

        for (
            const std::size_t budget :
            budgets
        ) {
            std::vector<double>
                strict_single_medians;

            std::vector<double>
                strict_pair_medians;

            std::vector<double>
                sensitive_single_medians;

            std::vector<double>
                sensitive_pair_medians;

            strict_single_medians.reserve(
                pairs.size()
            );

            strict_pair_medians.reserve(
                pairs.size()
            );

            sensitive_single_medians.reserve(
                pairs.size()
            );

            sensitive_pair_medians.reserve(
                pairs.size()
            );

            std::uint64_t
                strict_amplicons = 0;

            std::uint64_t
                sensitive_amplicons = 0;

            std::size_t
                hit_count_violations = 0;

            std::size_t
                candidate_routes = 0;

            std::size_t
                exhaustive_routes = 0;

            std::size_t
                estimator_used_routes = 0;

            for (
                std::size_t pair_index = 0;
                pair_index < pairs.size();
                ++pair_index
            ) {
                const auto& pair =
                    pairs.at(
                        pair_index
                    );

                std::vector<double>
                    strict_single_times;

                std::vector<double>
                    strict_pair_times;

                std::vector<double>
                    sensitive_single_times;

                std::vector<double>
                    sensitive_pair_times;

                strict_single_times.reserve(
                    repeats
                );

                strict_pair_times.reserve(
                    repeats
                );

                sensitive_single_times.reserve(
                    repeats
                );

                sensitive_pair_times.reserve(
                    repeats
                );

                primerpair::
                    StrandAwarePrimerSearchResult
                        strict_p1;

                primerpair::
                    StrandAwarePrimerSearchResult
                        strict_p2;

                primerpair::
                    SensitiveAdaptiveSearchResult
                        sensitive_p1;

                primerpair::
                    SensitiveAdaptiveSearchResult
                        sensitive_p2;

                primerpair::
                    PrimerPairSearchResult
                        strict_pair;

                primerpair::
                    PrimerPairSearchResult
                        sensitive_pair;

                for (
                    std::size_t repeat = 0;
                    repeat < repeats;
                    ++repeat
                ) {
                    auto run_strict_single =
                        [&]() {
                            const auto start =
                                std::chrono::
                                    steady_clock::now();

                            strict_p1 =
                                strict_single_engine.search(
                                    pair.primer1,
                                    anchor_length,
                                    budget
                                );

                            strict_p2 =
                                strict_single_engine.search(
                                    pair.primer2,
                                    anchor_length,
                                    budget
                                );

                            const auto stop =
                                std::chrono::
                                    steady_clock::now();

                            strict_single_times.push_back(
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

                    auto run_sensitive_single =
                        [&]() {
                            const auto start =
                                std::chrono::
                                    steady_clock::now();

                            sensitive_p1 =
                                sensitive_single_engine.search(
                                    pair.primer1,
                                    budget
                                );

                            sensitive_p2 =
                                sensitive_single_engine.search(
                                    pair.primer2,
                                    budget
                                );

                            const auto stop =
                                std::chrono::
                                    steady_clock::now();

                            sensitive_single_times.push_back(
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

                    auto run_strict_pair =
                        [&]() {
                            const auto start =
                                std::chrono::
                                    steady_clock::now();

                            strict_pair =
                                pair_engine.search(
                                    pair.primer1,
                                    pair.primer2,
                                    primerpair::
                                        SearchProfile::
                                            Strict,
                                    anchor_length,
                                    budget,
                                    min_amplicon,
                                    max_amplicon
                                );

                            const auto stop =
                                std::chrono::
                                    steady_clock::now();

                            strict_pair_times.push_back(
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

                    auto run_sensitive_pair =
                        [&]() {
                            const auto start =
                                std::chrono::
                                    steady_clock::now();

                            sensitive_pair =
                                pair_engine.search(
                                    pair.primer1,
                                    pair.primer2,
                                    primerpair::
                                        SearchProfile::
                                            Sensitive,
                                    anchor_length,
                                    budget,
                                    min_amplicon,
                                    max_amplicon
                                );

                            const auto stop =
                                std::chrono::
                                    steady_clock::now();

                            sensitive_pair_times.push_back(
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
                     * Rotate execution order to reduce
                     * systematic cache/order bias.
                     */
                    switch (
                        (
                            pair_index +
                            repeat
                        )
                        %
                        4
                    ) {
                        case 0:
                            run_strict_single();
                            run_sensitive_single();
                            run_strict_pair();
                            run_sensitive_pair();
                            break;

                        case 1:
                            run_sensitive_single();
                            run_strict_pair();
                            run_sensitive_pair();
                            run_strict_single();
                            break;

                        case 2:
                            run_strict_pair();
                            run_sensitive_pair();
                            run_strict_single();
                            run_sensitive_single();
                            break;

                        default:
                            run_sensitive_pair();
                            run_strict_single();
                            run_sensitive_single();
                            run_strict_pair();
                            break;
                    }
                }

                const double
                    strict_single_median =
                        percentile(
                            strict_single_times,
                            0.50
                        );

                const double
                    strict_pair_median =
                        percentile(
                            strict_pair_times,
                            0.50
                        );

                const double
                    sensitive_single_median =
                        percentile(
                            sensitive_single_times,
                            0.50
                        );

                const double
                    sensitive_pair_median =
                        percentile(
                            sensitive_pair_times,
                            0.50
                        );

                strict_single_medians.push_back(
                    strict_single_median
                );

                strict_pair_medians.push_back(
                    strict_pair_median
                );

                sensitive_single_medians.push_back(
                    sensitive_single_median
                );

                sensitive_pair_medians.push_back(
                    sensitive_pair_median
                );

                /*
                 * Standalone and pair engines must see
                 * exactly the same single-primer hits.
                 */
                if (
                    strict_pair
                        .primer1_single_hit_count
                    !=
                    strict_p1.hit_count()
                    ||
                    strict_pair
                        .primer2_single_hit_count
                    !=
                    strict_p2.hit_count()
                ) {
                    ++hit_count_violations;
                }

                if (
                    sensitive_pair
                        .primer1_single_hit_count
                    !=
                    sensitive_p1.hit_count()
                    ||
                    sensitive_pair
                        .primer2_single_hit_count
                    !=
                    sensitive_p2.hit_count()
                ) {
                    ++hit_count_violations;
                }

                strict_amplicons +=
                    strict_pair
                        .amplicon_count();

                sensitive_amplicons +=
                    sensitive_pair
                        .amplicon_count();

                const auto count_backend =
                    [&](
                        const primerpair::
                            SensitiveAdaptiveSearchResult&
                                result
                    ) {
                        if (
                            result.estimator_used
                        ) {
                            ++estimator_used_routes;
                        }

                        if (
                            result.backend ==
                            primerpair::
                                SensitiveAdaptiveBackend::
                                    Candidate
                        ) {
                            ++candidate_routes;

                        } else {

                            ++exhaustive_routes;
                        }
                    };

                count_backend(
                    sensitive_p1
                );

                count_backend(
                    sensitive_p2
                );
            }

            const double
                strict_single_mean =
                    mean(
                        strict_single_medians
                    );

            const double
                strict_pair_mean =
                    mean(
                        strict_pair_medians
                    );

            const double
                sensitive_single_mean =
                    mean(
                        sensitive_single_medians
                    );

            const double
                sensitive_pair_mean =
                    mean(
                        sensitive_pair_medians
                    );

            const double
                strict_nonsearch =
                    positive_difference(
                        strict_pair_mean,
                        strict_single_mean
                    );

            const double
                sensitive_nonsearch =
                    positive_difference(
                        sensitive_pair_mean,
                        sensitive_single_mean
                    );

            const double
                strict_search_share =
                    strict_pair_mean > 0.0
                        ? 100.0 *
                            strict_single_mean /
                            strict_pair_mean
                        : 0.0;

            const double
                sensitive_search_share =
                    sensitive_pair_mean > 0.0
                        ? 100.0 *
                            sensitive_single_mean /
                            sensitive_pair_mean
                        : 0.0;

            const double
                sensitive_vs_strict_single =
                    strict_single_mean > 0.0
                        ? sensitive_single_mean /
                            strict_single_mean
                        : 0.0;

            const double
                sensitive_vs_strict_pair =
                    strict_pair_mean > 0.0
                        ? sensitive_pair_mean /
                            strict_pair_mean
                        : 0.0;

            const bool all_checks =
                hit_count_violations == 0;

            std::cout
                << budget
                << '\t'
                << strict_single_mean
                << '\t'
                << strict_pair_mean
                << '\t'
                << strict_nonsearch
                << '\t'
                << strict_search_share
                << '\t'
                << sensitive_single_mean
                << '\t'
                << sensitive_pair_mean
                << '\t'
                << sensitive_nonsearch
                << '\t'
                << sensitive_search_share
                << '\t'
                << sensitive_vs_strict_single
                << '\t'
                << sensitive_vs_strict_pair
                << '\t'
                << candidate_routes
                << '\t'
                << exhaustive_routes
                << '\t'
                << estimator_used_routes
                << '\t'
                << strict_amplicons
                << '\t'
                << sensitive_amplicons
                << '\t'
                << hit_count_violations
                << '\t'
                << percentile(
                       strict_single_medians,
                       0.50
                   )
                << '\t'
                << percentile(
                       strict_pair_medians,
                       0.50
                   )
                << '\t'
                << percentile(
                       sensitive_single_medians,
                       0.50
                   )
                << '\t'
                << percentile(
                       sensitive_pair_medians,
                       0.50
                   )
                << '\t'
                << percentile(
                       strict_pair_medians,
                       0.95
                   )
                << '\t'
                << percentile(
                       sensitive_pair_medians,
                       0.95
                   )
                << '\t'
                << (
                    all_checks
                        ? "YES"
                        : "NO"
                )
                << '\n';

            if (!all_checks) {
                throw std::logic_error(
                    "Standalone/pair single-hit "
                    "count mismatch."
                );
            }
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
