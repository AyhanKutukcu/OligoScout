#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/multiplex_primer_search.hpp"
#include "primerpair/multiplex_primer_search_v2.hpp"
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
#include <vector>

namespace {


std::string make_reference(
    const std::size_t length
) {
    constexpr char bases[] = {
        'A', 'C', 'G', 'T'
    };

    std::uint64_t state =
        0x8CB92BA72F3D8DD7ULL;

    std::string reference;

    reference.reserve(length);


    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        reference.push_back(
            bases[
                static_cast<std::size_t>(
                    (state >> 32U) &
                    3ULL
                )
            ]
        );
    }


    return reference;
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
    >& products
) {
    std::sort(
        products.begin(),
        products.end(),
        cross_less
    );

    products.erase(
        std::unique(
            products.begin(),
            products.end()
        ),
        products.end()
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


        SearchDifficultyEstimator estimator(
            bifm
        );


        BatchedAnchorLookup anchor_lookup(
            ipbwt,
            estimator
        );


        AnchorCandidateSearcher verifier(
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
            v1(
                hybrid_strand
            );


        MultiplexPrimerSearchEngineV2
            v2(
                hybrid_strand
            );


        /*
         * Clustered primer pool.
         *
         * Cross products should be common.
         */
        std::vector<std::string>
            forward_pool;

        std::vector<std::string>
            reverse_pool;


        for (
            std::size_t i = 0;
            i < 12;
            ++i
        ) {
            const std::size_t left =
                10000 +
                i * 85;


            const std::size_t right =
                10400 +
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


        /*
         * Mixed panel:
         *
         * - shared forward primers
         * - shared reverse primers
         * - duplicated intended pair
         * - unique pairs
         */
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
        > requests;


        requests.reserve(
            layout.size()
        );


        for (
            const auto [left, right] :
            layout
        ) {
            requests.push_back(
                MultiplexPrimerPairRequest{
                    forward_pool.at(left),
                    reverse_pool.at(right),

                    3,

                    50,
                    1500
                }
            );
        }


        const auto result_v1 =
            v1.search(
                requests,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        const auto result_v2 =
            v2.search(
                requests,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        if (
            result_v1.intended_pairs.size() !=
            result_v2.intended_pairs.size()
        ) {
            throw std::runtime_error(
                "V1/V2 intended result-count mismatch."
            );
        }


        std::size_t intended_checks = 0;


        for (
            std::size_t i = 0;
            i < result_v1
                .intended_pairs
                .size();
            ++i
        ) {
            const auto& expected =
                result_v1
                    .intended_pairs
                    .at(i);


            const auto& observed =
                result_v2
                    .intended_pairs
                    .at(i);


            if (
                expected.amplicons !=
                observed.amplicons
            ) {
                std::cerr
                    << "V2_INTENDED_MISMATCH"
                    << '\t'
                    << "pair="
                    << i
                    << '\n';

                return 1;
            }


            if (
                expected
                    .primer1_single_hit_count !=
                observed
                    .primer1_single_hit_count
                ||
                expected
                    .primer2_single_hit_count !=
                observed
                    .primer2_single_hit_count
            ) {
                throw std::runtime_error(
                    "V1/V2 intended hit-count mismatch."
                );
            }


            ++intended_checks;
        }


        auto expected_cross =
            result_v1.cross_amplicons;


        auto observed_cross =
            result_v2.cross_amplicons;


        normalize(
            expected_cross
        );

        normalize(
            observed_cross
        );


        if (
            expected_cross !=
            observed_cross
        ) {
            std::cerr
                << "V2_CROSS_MISMATCH"
                << '\t'
                << "expected="
                << expected_cross.size()
                << '\t'
                << "observed="
                << observed_cross.size()
                << '\n';

            return 1;
        }


        if (
            result_v1
                .stats
                .unique_primer_queries !=
            result_v2
                .stats
                .unique_primer_queries
        ) {
            throw std::runtime_error(
                "V1/V2 unique-primer count mismatch."
            );
        }


        if (
            result_v2
                .global_cross_stats
                .window_candidates ==
            0
        ) {
            throw std::runtime_error(
                "V2 global sweep produced "
                "zero window candidates."
            );
        }


        std::cout
            << "pairs\t"
            << requests.size()
            << '\n';


        std::cout
            << "unique_primers\t"
            << result_v2
                .stats
                .unique_primer_queries
            << '\n';


        std::cout
            << "intended_checks\t"
            << intended_checks
            << '\n';


        std::cout
            << "cross_amplicons\t"
            << observed_cross.size()
            << '\n';


        std::cout
            << "global_forward_hits\t"
            << result_v2
                .global_cross_stats
                .forward_hits
            << '\n';


        std::cout
            << "global_reverse_hits\t"
            << result_v2
                .global_cross_stats
                .reverse_hits
            << '\n';


        std::cout
            << "global_window_candidates\t"
            << result_v2
                .global_cross_stats
                .window_candidates
            << '\n';


        std::cout
            << "global_unique_products\t"
            << result_v2
                .global_cross_stats
                .unique_products
            << '\n';


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
