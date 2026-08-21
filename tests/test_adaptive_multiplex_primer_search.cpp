#include "primerpair/adaptive_multiplex_primer_search.hpp"

#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/multiplex_primer_search.hpp"
#include "primerpair/multiplex_primer_search_v2.hpp"
#include "primerpair/multiplex_strategy_router.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>


namespace {


void expect(
    const bool condition,
    const char* message
) {
    if (!condition) {
        throw std::runtime_error(
            message
        );
    }


    std::cout
        << "[PASS] "
        << message
        << '\n';
}


std::string make_reference(
    const std::size_t length
) {
    static constexpr char
        bases[] = {
            'A',
            'C',
            'G',
            'T'
        };


    std::string sequence;

    sequence.resize(
        length
    );


    std::uint64_t state =
        UINT64_C(
            0x9E3779B97F4A7C15
        );


    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state
            *
            UINT64_C(
                6364136223846793005
            )
            +
            UINT64_C(
                1442695040888963407
            );


        sequence.at(i) =
            bases[
                static_cast<
                    std::size_t
                >(
                    (
                        state >>
                        32
                    )
                    &
                    UINT64_C(3)
                )
            ];
    }


    return sequence;
}


bool cross_less(
    const primerpair::
        MultiplexCrossAmplicon& lhs,

    const primerpair::
        MultiplexCrossAmplicon& rhs
) noexcept {

    return
        std::tie(
            lhs.forward_pair_index,
            lhs.forward_slot,
            lhs.reverse_pair_index,
            lhs.reverse_slot,
            lhs.amplicon_start,
            lhs.amplicon_end_exclusive,
            lhs.amplicon_length,
            lhs.forward_mismatches,
            lhs.reverse_mismatches,
            lhs.forward_mismatch_mask,
            lhs.reverse_mismatch_mask
        )
        <
        std::tie(
            rhs.forward_pair_index,
            rhs.forward_slot,
            rhs.reverse_pair_index,
            rhs.reverse_slot,
            rhs.amplicon_start,
            rhs.amplicon_end_exclusive,
            rhs.amplicon_length,
            rhs.forward_mismatches,
            rhs.reverse_mismatches,
            rhs.forward_mismatch_mask,
            rhs.reverse_mismatch_mask
        );
}


void normalize(
    std::vector<
        primerpair::
            MultiplexCrossAmplicon
    >& values
) {
    std::sort(
        values.begin(),
        values.end(),
        cross_less
    );


    values.erase(
        std::unique(
            values.begin(),
            values.end()
        ),
        values.end()
    );
}


void expect_intended_equal(
    const std::vector<
        primerpair::
            PrimerPairSearchResult
    >& expected,

    const std::vector<
        primerpair::
            PrimerPairSearchResult
    >& observed
) {
    expect(
        expected.size() ==
            observed.size(),
        "Intended result counts equal"
    );


    for (
        std::size_t i = 0;
        i < expected.size();
        ++i
    ) {
        const auto& a =
            expected.at(i);

        const auto& b =
            observed.at(i);


        expect(
            a.amplicons ==
                b.amplicons,
            "Intended amplicons exact"
        );


        expect(
            a.primer1_single_hit_count ==
                b.primer1_single_hit_count,
            "Primer1 single-hit counts exact"
        );


        expect(
            a.primer2_single_hit_count ==
                b.primer2_single_hit_count,
            "Primer2 single-hit counts exact"
        );
    }
}


void expect_cross_equal(
    std::vector<
        primerpair::
            MultiplexCrossAmplicon
    > expected,

    std::vector<
        primerpair::
            MultiplexCrossAmplicon
    > observed
) {
    normalize(
        expected
    );

    normalize(
        observed
    );


    expect(
        expected ==
            observed,
        "Cross-amplicons exact"
    );
}


}  // namespace


int main() {
    try {
        using namespace primerpair;


        constexpr std::size_t
            anchor_length = 12;


        constexpr std::size_t
            primer_length = 24;


        constexpr std::uint64_t
            cross_min = 50;


        constexpr std::uint64_t
            cross_max = 1500;


        /*
         * Same general construction used by the existing
         * V1/V2 equivalence test.
         */
        const std::string reference =
            make_reference(
                150000
            );


        PackedReference packed(
            reference
        );


        BidirectionalFMIndex bifm(
            reference
        );


        IPBWTIndex ipbwt(
            reference,
            anchor_length
        );


        SearchDifficultyEstimator
            estimator(
                bifm
            );


        BatchedAnchorLookup
            anchor_lookup(
                ipbwt,
                estimator
            );


        AnchorCandidateSearcher
            verifier(
                bifm,
                packed
            );


        BatchedCandidateSearchEngine
            candidate_engine(
                ipbwt,
                anchor_lookup,
                verifier
            );


        SinglePrimerSearchEngine
            single_engine(
                bifm,
                packed
            );


        HybridBatchedPrimerSearchEngine
            forward_hybrid(
                candidate_engine,
                single_engine
            );


        StrandAwarePrimerSearchEngine
            legacy_strand(
                bifm,
                packed
            );


        HybridStrandAwarePrimerSearchEngine
            hybrid_strand(
                forward_hybrid,
                anchor_lookup,
                ipbwt,
                packed,
                legacy_strand
            );


        MultiplexPrimerSearchEngine
            production_v1(
                hybrid_strand
            );


        MultiplexPrimerSearchEngineV2
            global_v2(
                hybrid_strand
            );


        AdaptiveMultiplexPrimerSearchEngine
            adaptive(
                production_v1,
                global_v2
            );


        expect(
            adaptive
                .router()
                .crossover_pair_count()
            ==
            12,
            "Adaptive engine uses production threshold 12"
        );


        /*
         * Build a clustered 12-pair panel.
         *
         * This intentionally includes repeated primer
         * memberships, matching the stress pattern already
         * used in V1/V2 equivalence validation.
         */
        std::vector<
            std::string
        > forward_pool;


        std::vector<
            std::string
        > reverse_pool;


        for (
            std::size_t i = 0;
            i < 12;
            ++i
        ) {
            const std::size_t left =
                10000
                +
                i * 85;


            const std::size_t right =
                10400
                +
                i * 91;


            forward_pool.push_back(
                reference.substr(
                    left,
                    primer_length
                )
            );


            reverse_pool.push_back(
                reverse_complement(
                    reference.substr(
                        right,
                        primer_length
                    )
                )
            );
        }


        const std::vector<
            std::pair<
                std::size_t,
                std::size_t
            >
        > layout{
            {0, 0},
            {0, 1},
            {2, 0},
            {0, 0},
            {3, 3},
            {4, 0},
            {5, 5},
            {0, 6},
            {7, 7},
            {8, 1},
            {9, 9},
            {10, 10}
        };


        std::vector<
            MultiplexPrimerPairRequest
        > full_requests;


        full_requests.reserve(
            layout.size()
        );


        for (
            const auto [left, right] :
            layout
        ) {
            full_requests.push_back(
                MultiplexPrimerPairRequest{
                    forward_pool.at(left),
                    reverse_pool.at(right),
                    3,
                    50,
                    1500
                }
            );
        }


        expect(
            full_requests.size() == 12,
            "Synthetic full panel contains 12 pairs"
        );


        /*
         * ================================================
         * TEST A:
         *
         * 11 pairs -> Production V1
         * ================================================
         */
        std::vector<
            MultiplexPrimerPairRequest
        > requests_11(
            full_requests.begin(),
            full_requests.begin() + 11
        );


        const auto direct_v1 =
            production_v1.search(
                requests_11,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        const auto adaptive_11 =
            adaptive.search(
                requests_11,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        expect(
            adaptive_11
                .decision
                .strategy
            ==
            MultiplexExecutionStrategy::
                ProductionV1,
            "11 pairs route to Production V1"
        );


        expect(
            adaptive_11
                .used_production_v1(),
            "V1 helper reports selected backend"
        );


        expect(
            !adaptive_11
                .used_global_sweep_v2(),
            "V2 helper false for 11 pairs"
        );


        expect(
            !adaptive_11
                .global_cross_stats_available,
            "V1 result has no global sweep stats"
        );


        expect_intended_equal(
            direct_v1.intended_pairs,
            adaptive_11.intended_pairs
        );


        expect_cross_equal(
            direct_v1.cross_amplicons,
            adaptive_11.cross_amplicons
        );


        expect(
            direct_v1
                .stats
                .pair_requests
            ==
            adaptive_11
                .stats
                .pair_requests,
            "V1 pair-request stats exact"
        );


        expect(
            direct_v1
                .stats
                .unique_primer_queries
            ==
            adaptive_11
                .stats
                .unique_primer_queries,
            "V1 unique-primer stats exact"
        );


        expect(
            direct_v1
                .stats
                .cross_amplicon_records
            ==
            adaptive_11
                .stats
                .cross_amplicon_records,
            "V1 cross-record stats exact"
        );


        /*
         * ================================================
         * TEST B:
         *
         * 12 pairs -> Global Sweep V2
         * ================================================
         */
        const auto direct_v2 =
            global_v2.search(
                full_requests,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        const auto adaptive_12 =
            adaptive.search(
                full_requests,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        expect(
            adaptive_12
                .decision
                .strategy
            ==
            MultiplexExecutionStrategy::
                GlobalSweepV2,
            "12 pairs route to Global Sweep V2"
        );


        expect(
            adaptive_12
                .used_global_sweep_v2(),
            "V2 helper reports selected backend"
        );


        expect(
            !adaptive_12
                .used_production_v1(),
            "V1 helper false for 12 pairs"
        );


        expect(
            adaptive_12
                .global_cross_stats_available,
            "V2 result exposes global sweep stats"
        );


        expect_intended_equal(
            direct_v2.intended_pairs,
            adaptive_12.intended_pairs
        );


        expect_cross_equal(
            direct_v2.cross_amplicons,
            adaptive_12.cross_amplicons
        );


        expect(
            direct_v2
                .stats
                .pair_requests
            ==
            adaptive_12
                .stats
                .pair_requests,
            "V2 pair-request stats exact"
        );


        expect(
            direct_v2
                .stats
                .unique_primer_queries
            ==
            adaptive_12
                .stats
                .unique_primer_queries,
            "V2 unique-primer stats exact"
        );


        expect(
            direct_v2
                .stats
                .cross_amplicon_records
            ==
            adaptive_12
                .stats
                .cross_amplicon_records,
            "V2 cross-record stats exact"
        );


        expect(
            direct_v2
                .global_cross_stats
                .unique_primers
            ==
            adaptive_12
                .global_cross_stats
                .unique_primers,
            "Global unique-primer stats exact"
        );


        expect(
            direct_v2
                .global_cross_stats
                .forward_hits
            ==
            adaptive_12
                .global_cross_stats
                .forward_hits,
            "Global forward-hit stats exact"
        );


        expect(
            direct_v2
                .global_cross_stats
                .reverse_hits
            ==
            adaptive_12
                .global_cross_stats
                .reverse_hits,
            "Global reverse-hit stats exact"
        );


        expect(
            direct_v2
                .global_cross_stats
                .window_candidates
            ==
            adaptive_12
                .global_cross_stats
                .window_candidates,
            "Global window-candidate stats exact"
        );


        expect(
            direct_v2
                .global_cross_stats
                .unique_products
            ==
            adaptive_12
                .global_cross_stats
                .unique_products,
            "Global unique-product stats exact"
        );


        /*
         * ================================================
         * TEST C:
         *
         * Custom router override.
         *
         * threshold 5:
         * 5 pairs -> V2
         * ================================================
         */
        AdaptiveMultiplexPrimerSearchEngine
            custom_adaptive(
                production_v1,
                global_v2,
                MultiplexStrategyRouter{
                    5
                }
            );


        std::vector<
            MultiplexPrimerPairRequest
        > requests_5(
            full_requests.begin(),
            full_requests.begin() + 5
        );


        const auto custom_result =
            custom_adaptive.search(
                requests_5,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        expect(
            custom_result
                .decision
                .strategy
            ==
            MultiplexExecutionStrategy::
                GlobalSweepV2,
            "Custom threshold 5 routes five pairs to V2"
        );


        /*
         * ================================================
         * TEST D:
         *
         * Empty input.
         *
         * Router naturally sends zero-pair workload
         * to V1.
         * ================================================
         */
        const std::vector<
            MultiplexPrimerPairRequest
        > empty_requests;


        const auto empty_result =
            adaptive.search(
                empty_requests,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        expect(
            empty_result
                .decision
                .strategy
            ==
            MultiplexExecutionStrategy::
                ProductionV1,
            "Empty panel routes to V1"
        );


        expect(
            empty_result
                .intended_pairs
                .empty(),
            "Empty panel has no intended products"
        );


        expect(
            empty_result
                .cross_amplicons
                .empty(),
            "Empty panel has no cross products"
        );


        std::cout
            << "production_threshold\t"
            << adaptive
                .router()
                .crossover_pair_count()
            << '\n';


        std::cout
            << "v1_path_pairs\t"
            << requests_11.size()
            << '\n';


        std::cout
            << "v2_path_pairs\t"
            << full_requests.size()
            << '\n';


        std::cout
            << "ADAPTIVE_V1_EXACT_EQUIVALENCE\tYES\n";


        std::cout
            << "ADAPTIVE_V2_EXACT_EQUIVALENCE\tYES\n";


        std::cout
            << "ADAPTIVE_ROUTER_BOUNDARY\tYES\n";


        std::cout
            << "ADAPTIVE_CUSTOM_ROUTER\tYES\n";


        std::cout
            << "ALL_CHECKS\tYES\n";


        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';


        return 1;
    }
}
