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

#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/sensitive_pair_constrained_search.hpp"

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
        [](const char base) {
            return
                base == 'A' ||
                base == 'C' ||
                base == 'G' ||
                base == 'T';
        }
    );
}

bool contains_target(
    const primerpair::PrimerPairSearchResult& result,
    const DesignedPair& pair
) {
    return std::any_of(
        result.amplicons.begin(),
        result.amplicons.end(),
        [&pair](const auto& hit) {
            return
                hit.left_primer ==
                    primerpair::
                        PrimerIdentity::
                            Primer1
                &&
                hit.right_primer ==
                    primerpair::
                        PrimerIdentity::
                            Primer2
                &&
                hit.left_position ==
                    pair.left_position
                &&
                hit.right_position ==
                    pair.right_position
                &&
                hit.amplicon_length ==
                    pair.amplicon_length;
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
                << " <fasta> [pair_count] [repeats]\n";

            return 2;
        }

        constexpr std::size_t primer_length = 20;
        constexpr std::size_t anchor_length = 12;
        constexpr std::size_t budget = 3;

        constexpr std::uint64_t min_amplicon = 50;
        constexpr std::uint64_t max_amplicon = 3000;

        /*
         * Same deterministic sampler used in the
         * earlier chr22 pair benchmarks.
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

        std::vector<DesignedPair> pairs;

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
                !all_acgt(left_site) ||
                !all_acgt(right_site)
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

        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        /*
         * Baseline:
         * both primers searched globally.
         */
        const primerpair::PrimerPairSearchEngine
            global_engine(
                index,
                packed
            );

        /*
         * H1 prototype:
         * one global anchor + local mate scan.
         */
        const primerpair::
            SensitivePairConstrainedSearchEngine
                constrained_engine(
                    index,
                    packed
                );

        std::vector<double>
            global_pair_medians;

        std::vector<double>
            constrained_pair_medians;

        global_pair_medians.reserve(
            pairs.size()
        );

        constrained_pair_medians.reserve(
            pairs.size()
        );


        std::vector<double>
            forced_p1_pair_medians;

        std::vector<double>
            forced_p2_pair_medians;

        std::vector<double>
            oracle_pair_medians;

        forced_p1_pair_medians.reserve(
            pairs.size()
        );

        forced_p2_pair_medians.reserve(
            pairs.size()
        );

        oracle_pair_medians.reserve(
            pairs.size()
        );

        struct OraclePairRow {
            std::size_t pair_index{0};

            double adaptive_us{0.0};
            double forced_p1_us{0.0};
            double forced_p2_us{0.0};
            double oracle_us{0.0};

            std::uint64_t
                p1_max_seed_occurrences{0};

            std::uint64_t
                p2_max_seed_occurrences{0};

            std::uint64_t
                p1_total_seed_occurrences{0};

            std::uint64_t
                p2_total_seed_occurrences{0};

            primerpair::PrimerIdentity
                adaptive_anchor{
                    primerpair::
                        PrimerIdentity::
                            Primer1
                };

            primerpair::PrimerIdentity
                oracle_anchor{
                    primerpair::
                        PrimerIdentity::
                            Primer1
                };
        };

        std::vector<OraclePairRow>
            oracle_rows;

        oracle_rows.reserve(
            pairs.size()
        );

        std::uint64_t global_amplicons = 0;
        std::uint64_t constrained_amplicons = 0;

        std::uint64_t scanned_mate_starts = 0;
        std::uint64_t global_mate_hits = 0;
        std::uint64_t local_mate_hits = 0;

        std::uint64_t
            merged_forward_windows = 0;

        std::uint64_t
            merged_reverse_windows = 0;

        std::size_t anchor_primer1 = 0;
        std::size_t anchor_primer2 = 0;

        std::size_t
            equality_violations = 0;

        std::size_t
            anchor_hit_count_violations = 0;

        std::size_t
            local_mate_count_violations = 0;

        std::size_t global_targets = 0;
        std::size_t constrained_targets = 0;


        std::size_t
            forced_p1_equality_violations = 0;

        std::size_t
            forced_p2_equality_violations = 0;

        std::size_t
            oracle_primer1_count = 0;

        std::size_t
            oracle_primer2_count = 0;

        std::size_t
            adaptive_oracle_matches = 0;

        long double timing_cost_us = 0.0L;
        long double timing_anchor_us = 0.0L;
        long double timing_windows_us = 0.0L;
        long double timing_seed_us = 0.0L;
        long double timing_filter_verify_us = 0.0L;
        long double timing_pair_us = 0.0L;
        long double timing_accounted_us = 0.0L;

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
                global_times;

            std::vector<double>
                constrained_times;


            std::vector<double>
                forced_p1_times;

            std::vector<double>
                forced_p2_times;

            global_times.reserve(
                repeats
            );

            constrained_times.reserve(
                repeats
            );


            forced_p1_times.reserve(
                repeats
            );

            forced_p2_times.reserve(
                repeats
            );

            primerpair::
                PrimerPairSearchResult
                    global_result;

            primerpair::
                SensitivePairConstrainedSearchResult
                    constrained_result;


            primerpair::
                SensitivePairConstrainedSearchResult
                    forced_p1_result;

            primerpair::
                SensitivePairConstrainedSearchResult
                    forced_p2_result;

            for (
                std::size_t repeat = 0;
                repeat < repeats;
                ++repeat
            ) {
                auto run_global =
                    [&]() {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        global_result =
                            global_engine.search(
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

                        global_times.push_back(
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

                auto run_constrained =
                    [&]() {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        constrained_result =
                            constrained_engine.search(
                                pair.primer1,
                                pair.primer2,
                                budget,
                                min_amplicon,
                                max_amplicon
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        constrained_times.push_back(
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

                auto run_forced_p1 =
                    [&]() {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        forced_p1_result =
                            constrained_engine.search(
                                pair.primer1,
                                pair.primer2,
                                budget,
                                min_amplicon,
                                max_amplicon,
                                primerpair::
                                    SensitivePairAnchorPolicy::
                                        ForcePrimer1
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        forced_p1_times.push_back(
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


                auto run_forced_p2 =
                    [&]() {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        forced_p2_result =
                            constrained_engine.search(
                                pair.primer1,
                                pair.primer2,
                                budget,
                                min_amplicon,
                                max_amplicon,
                                primerpair::
                                    SensitivePairAnchorPolicy::
                                        ForcePrimer2
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        forced_p2_times.push_back(
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
                 * Four-way rotation reduces systematic
                 * cache/order bias between:
                 *
                 *   global baseline
                 *   adaptive anchor
                 *   forced Primer1
                 *   forced Primer2
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
                        run_global();
                        run_constrained();
                        run_forced_p1();
                        run_forced_p2();
                        break;

                    case 1:
                        run_constrained();
                        run_forced_p1();
                        run_forced_p2();
                        run_global();
                        break;

                    case 2:
                        run_forced_p1();
                        run_forced_p2();
                        run_global();
                        run_constrained();
                        break;

                    default:
                        run_forced_p2();
                        run_global();
                        run_constrained();
                        run_forced_p1();
                        break;
                }
            }

            global_pair_medians.push_back(
                percentile(
                    global_times,
                    0.50
                )
            );

            constrained_pair_medians.push_back(
                percentile(
                    constrained_times,
                    0.50
                )
            );


            const double forced_p1_median =
                percentile(
                    forced_p1_times,
                    0.50
                );

            const double forced_p2_median =
                percentile(
                    forced_p2_times,
                    0.50
                );

            forced_p1_pair_medians.push_back(
                forced_p1_median
            );

            forced_p2_pair_medians.push_back(
                forced_p2_median
            );

            const double oracle_median =
                std::min(
                    forced_p1_median,
                    forced_p2_median
                );

            oracle_pair_medians.push_back(
                oracle_median
            );

            const primerpair::PrimerIdentity
                oracle_anchor =
                    forced_p1_median <=
                            forced_p2_median
                        ?
                        primerpair::
                            PrimerIdentity::
                                Primer1
                        :
                        primerpair::
                            PrimerIdentity::
                                Primer2;

            if (
                oracle_anchor ==
                primerpair::
                    PrimerIdentity::
                        Primer1
            ) {
                ++oracle_primer1_count;

            } else {

                ++oracle_primer2_count;
            }

            if (
                constrained_result
                    .anchor_primer ==
                oracle_anchor
            ) {
                ++adaptive_oracle_matches;
            }

            oracle_rows.push_back(
                OraclePairRow{
                    pair_index,

                    constrained_pair_medians.back(),
                    forced_p1_median,
                    forced_p2_median,
                    oracle_median,

                    constrained_result
                        .primer1_cost
                        .max_seed_occurrences,

                    constrained_result
                        .primer2_cost
                        .max_seed_occurrences,

                    constrained_result
                        .primer1_cost
                        .total_seed_occurrences,

                    constrained_result
                        .primer2_cost
                        .total_seed_occurrences,

                    constrained_result
                        .anchor_primer,

                    oracle_anchor
                }
            );

            /*
             * Critical lossless invariant.
             */
            if (
                constrained_result
                    .pair_result
                    .amplicons
                !=
                global_result
                    .amplicons
            ) {
                ++equality_violations;
            }

            if (
                forced_p1_result
                    .pair_result
                    .amplicons
                !=
                global_result
                    .amplicons
            ) {
                ++forced_p1_equality_violations;
            }

            if (
                forced_p2_result
                    .pair_result
                    .amplicons
                !=
                global_result
                    .amplicons
            ) {
                ++forced_p2_equality_violations;
            }


            if (
                contains_target(
                    global_result,
                    pair
                )
            ) {
                ++global_targets;
            }

            if (
                contains_target(
                    constrained_result
                        .pair_result,
                    pair
                )
            ) {
                ++constrained_targets;
            }

            global_amplicons +=
                global_result
                    .amplicon_count();

            constrained_amplicons +=
                constrained_result
                    .pair_result
                    .amplicon_count();

            scanned_mate_starts +=
                constrained_result
                    .scanned_mate_start_positions;

            local_mate_hits +=
                constrained_result
                    .mate_local_hit_count;

            timing_cost_us +=
                constrained_result
                    .timing
                    .cost_estimator_us;

            timing_anchor_us +=
                constrained_result
                    .timing
                    .anchor_search_us;

            timing_windows_us +=
                constrained_result
                    .timing
                    .window_build_merge_us;

            timing_seed_us +=
                constrained_result
                    .timing
                    .mate_seed_generation_us;

            timing_filter_verify_us +=
                constrained_result
                    .timing
                    .mate_filter_verify_us;

            timing_pair_us +=
                constrained_result
                    .timing
                    .pair_assembly_us;

            timing_accounted_us +=
                constrained_result
                    .timing
                    .accounted_total_us();

            merged_forward_windows +=
                constrained_result
                    .merged_forward_mate_windows;

            merged_reverse_windows +=
                constrained_result
                    .merged_reverse_mate_windows;

            std::uint64_t
                expected_anchor_hits = 0;

            std::uint64_t
                expected_mate_hits = 0;

            if (
                constrained_result
                    .anchor_primer ==
                primerpair::
                    PrimerIdentity::
                        Primer1
            ) {
                ++anchor_primer1;

                expected_anchor_hits =
                    global_result
                        .primer1_single_hit_count;

                expected_mate_hits =
                    global_result
                        .primer2_single_hit_count;

            } else {

                ++anchor_primer2;

                expected_anchor_hits =
                    global_result
                        .primer2_single_hit_count;

                expected_mate_hits =
                    global_result
                        .primer1_single_hit_count;
            }

            global_mate_hits +=
                expected_mate_hits;

            if (
                constrained_result
                    .anchor_global_hit_count
                !=
                expected_anchor_hits
            ) {
                ++anchor_hit_count_violations;
            }

            /*
             * Every local mate hit must also be a
             * valid genome-wide SENSITIVE mate hit.
             */
            if (
                constrained_result
                    .mate_local_hit_count
                >
                expected_mate_hits
            ) {
                ++local_mate_count_violations;
            }
        }

        const double global_mean =
            mean(
                global_pair_medians
            );

        const double constrained_mean =
            mean(
                constrained_pair_medians
            );

        const double speedup =
            constrained_mean > 0.0
                ? global_mean /
                    constrained_mean
                : 0.0;


        const double forced_p1_mean =
            mean(
                forced_p1_pair_medians
            );

        const double forced_p2_mean =
            mean(
                forced_p2_pair_medians
            );

        const double oracle_mean =
            mean(
                oracle_pair_medians
            );

        const double oracle_median =
            percentile(
                oracle_pair_medians,
                0.50
            );

        const double oracle_p95 =
            percentile(
                oracle_pair_medians,
                0.95
            );

        const double adaptive_p95 =
            percentile(
                constrained_pair_medians,
                0.95
            );

        const double
            oracle_speedup_vs_adaptive =
                oracle_mean > 0.0
                    ?
                    constrained_mean /
                        oracle_mean
                    :
                    0.0;

        const double
            oracle_speedup_vs_global =
                oracle_mean > 0.0
                    ?
                    global_mean /
                        oracle_mean
                    :
                    0.0;

        const double
            oracle_gain_vs_adaptive_pct =
                constrained_mean > 0.0
                    ?
                    100.0 *
                    (
                        1.0 -
                        oracle_mean /
                            constrained_mean
                    )
                    :
                    0.0;

        const double
            adaptive_oracle_accuracy_pct =
                pairs.empty()
                    ?
                    0.0
                    :
                    100.0 *
                    static_cast<double>(
                        adaptive_oracle_matches
                    )
                    /
                    static_cast<double>(
                        pairs.size()
                    );

        const long double
            possible_oriented_starts =
                static_cast<long double>(
                    pairs.size()
                )
                *
                2.0L
                *
                static_cast<long double>(
                    reference.size() -
                    primer_length +
                    1
                );

        const long double
            local_scan_pct =
                possible_oriented_starts > 0.0L
                    ?
                    100.0L *
                    static_cast<long double>(
                        scanned_mate_starts
                    )
                    /
                    possible_oriented_starts
                    :
                    0.0L;

        const long double
            local_mate_fraction =
                global_mate_hits > 0
                    ?
                    static_cast<long double>(
                        local_mate_hits
                    )
                    /
                    static_cast<long double>(
                        global_mate_hits
                    )
                    :
                    0.0L;

        const double
            mean_scanned_starts =
                static_cast<double>(
                    scanned_mate_starts
                )
                /
                static_cast<double>(
                    pairs.size()
                );

        const double
            mean_merged_windows =
                static_cast<double>(
                    merged_forward_windows +
                    merged_reverse_windows
                )
                /
                static_cast<double>(
                    pairs.size()
                );

        const bool all_checks =
            equality_violations == 0
            &&
            forced_p1_equality_violations == 0
            &&
            forced_p2_equality_violations == 0
            &&
            anchor_hit_count_violations == 0
            &&
            local_mate_count_violations == 0
            &&
            global_targets ==
                pairs.size()
            &&
            constrained_targets ==
                pairs.size();

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
            << "# budget\t"
            << budget
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
            << "global_amplicons\t"
            << global_amplicons
            << '\n';

        std::cout
            << "constrained_amplicons\t"
            << constrained_amplicons
            << '\n';

        std::cout
            << "amplicon_equality_violations\t"
            << equality_violations
            << '\n';

        std::cout
            << "global_targets\t"
            << global_targets
            << '/'
            << pairs.size()
            << '\n';

        std::cout
            << "constrained_targets\t"
            << constrained_targets
            << '/'
            << pairs.size()
            << '\n';

        std::cout
            << "anchor_primer1\t"
            << anchor_primer1
            << '\n';

        std::cout
            << "anchor_primer2\t"
            << anchor_primer2
            << '\n';

        std::cout
            << "anchor_hit_count_violations\t"
            << anchor_hit_count_violations
            << '\n';

        std::cout
            << "local_mate_count_violations\t"
            << local_mate_count_violations
            << '\n';

        std::cout
            << "global_mate_hits\t"
            << global_mate_hits
            << '\n';

        std::cout
            << "local_mate_hits\t"
            << local_mate_hits
            << '\n';

        std::cout
            << "local_mate_hit_fraction\t"
            << static_cast<double>(
                   local_mate_fraction
               )
            << '\n';

        std::cout
            << "scanned_mate_starts\t"
            << scanned_mate_starts
            << '\n';

        std::cout
            << "mean_scanned_starts_per_pair\t"
            << mean_scanned_starts
            << '\n';

        std::cout
            << "effective_local_scan_pct\t"
            << static_cast<double>(
                   local_scan_pct
               )
            << '\n';

        std::cout
            << "mean_merged_windows_per_pair\t"
            << mean_merged_windows
            << '\n';

        const long double timing_denominator =
            static_cast<long double>(
                pairs.size()
            );

        std::cout
            << "timing_cost_estimator_mean_us\t"
            << static_cast<double>(
                   timing_cost_us /
                   timing_denominator
               )
            << '\n';

        std::cout
            << "timing_anchor_search_mean_us\t"
            << static_cast<double>(
                   timing_anchor_us /
                   timing_denominator
               )
            << '\n';

        std::cout
            << "timing_window_build_merge_mean_us\t"
            << static_cast<double>(
                   timing_windows_us /
                   timing_denominator
               )
            << '\n';

        std::cout
            << "timing_mate_seed_generation_mean_us\t"
            << static_cast<double>(
                   timing_seed_us /
                   timing_denominator
               )
            << '\n';

        std::cout
            << "timing_mate_filter_verify_mean_us\t"
            << static_cast<double>(
                   timing_filter_verify_us /
                   timing_denominator
               )
            << '\n';

        std::cout
            << "timing_pair_assembly_mean_us\t"
            << static_cast<double>(
                   timing_pair_us /
                   timing_denominator
               )
            << '\n';

        std::cout
            << "timing_accounted_mean_us\t"
            << static_cast<double>(
                   timing_accounted_us /
                   timing_denominator
               )
            << '\n';

        std::cout
            << "forced_p1_equality_violations\t"
            << forced_p1_equality_violations
            << '\n';

        std::cout
            << "forced_p2_equality_violations\t"
            << forced_p2_equality_violations
            << '\n';

        std::cout
            << "forced_p1_mean_us\t"
            << forced_p1_mean
            << '\n';

        std::cout
            << "forced_p2_mean_us\t"
            << forced_p2_mean
            << '\n';

        std::cout
            << "oracle_mean_us\t"
            << oracle_mean
            << '\n';

        std::cout
            << "oracle_median_us\t"
            << oracle_median
            << '\n';

        std::cout
            << "oracle_p95_us\t"
            << oracle_p95
            << '\n';

        std::cout
            << "adaptive_p95_us\t"
            << adaptive_p95
            << '\n';

        std::cout
            << "oracle_speedup_vs_adaptive\t"
            << oracle_speedup_vs_adaptive
            << '\n';

        std::cout
            << "oracle_speedup_vs_global\t"
            << oracle_speedup_vs_global
            << '\n';

        std::cout
            << "oracle_gain_vs_adaptive_pct\t"
            << oracle_gain_vs_adaptive_pct
            << '\n';

        std::cout
            << "adaptive_oracle_matches\t"
            << adaptive_oracle_matches
            << '/'
            << pairs.size()
            << '\n';

        std::cout
            << "adaptive_oracle_accuracy_pct\t"
            << adaptive_oracle_accuracy_pct
            << '\n';

        std::cout
            << "oracle_primer1_count\t"
            << oracle_primer1_count
            << '\n';

        std::cout
            << "oracle_primer2_count\t"
            << oracle_primer2_count
            << '\n';

        std::cout
            << "global_mean_us\t"
            << global_mean
            << '\n';

        std::cout
            << "constrained_mean_us\t"
            << constrained_mean
            << '\n';

        std::cout
            << "speedup\t"
            << speedup
            << '\n';

        std::cout
            << "global_median_us\t"
            << percentile(
                   global_pair_medians,
                   0.50
               )
            << '\n';

        std::cout
            << "constrained_median_us\t"
            << percentile(
                   constrained_pair_medians,
                   0.50
               )
            << '\n';

        std::cout
            << "global_p95_us\t"
            << percentile(
                   global_pair_medians,
                   0.95
               )
            << '\n';

        std::cout
            << "constrained_p95_us\t"
            << percentile(
                   constrained_pair_medians,
                   0.95
               )
            << '\n';

        std::cout
            << "# oracle_pair_columns\t"
            << "pair_index"
            << "\tadaptive_us"
            << "\tforced_p1_us"
            << "\tforced_p2_us"
            << "\toracle_us"
            << "\tp1_max_seed"
            << "\tp2_max_seed"
            << "\tp1_total_seed"
            << "\tp2_total_seed"
            << "\tadaptive_anchor"
            << "\toracle_anchor"
            << '\n';

        for (
            const auto& row :
            oracle_rows
        ) {
            std::cout
                << "PAIR\t"
                << row.pair_index
                << '\t'
                << row.adaptive_us
                << '\t'
                << row.forced_p1_us
                << '\t'
                << row.forced_p2_us
                << '\t'
                << row.oracle_us
                << '\t'
                << row.p1_max_seed_occurrences
                << '\t'
                << row.p2_max_seed_occurrences
                << '\t'
                << row.p1_total_seed_occurrences
                << '\t'
                << row.p2_total_seed_occurrences
                << '\t'
                << (
                    row.adaptive_anchor ==
                    primerpair::
                        PrimerIdentity::
                            Primer1
                        ?
                        "P1"
                        :
                        "P2"
                )
                << '\t'
                << (
                    row.oracle_anchor ==
                    primerpair::
                        PrimerIdentity::
                            Primer1
                        ?
                        "P1"
                        :
                        "P2"
                )
                << '\n';
        }

        std::cout
            << "all_checks\t"
            << (
                all_checks
                    ? "YES"
                    : "NO"
            )
            << '\n';

        if (!all_checks) {
            throw std::logic_error(
                "Pair-constrained real-FASTA "
                "differential validation failed."
            );
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
