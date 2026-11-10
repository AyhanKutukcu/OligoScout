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
#include <string_view>
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
        0xDB4F0B9175AE2165ULL;

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

        constexpr std::uint64_t
            cross_min = 50;

        constexpr std::uint64_t
            cross_max = 15000;


        std::string reference =
            make_reference(
                100000
            );


        /*
         * Two primer sequences deliberately placed
         * in BOTH genomic orientations.
         *
         * Therefore the same unique primer ID can
         * contribute both Forward and Reverse hits.
         */
        const std::string primer_a =
            "AGTCCGATGCTAACGTTGACCTGA";

        const std::string primer_b =
            "CTGATCGGACATGTCAGTACCGTA";


        const std::string rc_a =
            reverse_complement(
                primer_a
            );

        const std::string rc_b =
            reverse_complement(
                primer_b
            );


        /*
         * Primer A:
         * Forward hits at 10000, 30000
         * Reverse hits at 10400, 30400
         */
        reference.replace(
            10000,
            primer_a.size(),
            primer_a
        );

        reference.replace(
            10400,
            rc_a.size(),
            rc_a
        );

        reference.replace(
            30000,
            primer_a.size(),
            primer_a
        );

        reference.replace(
            30400,
            rc_a.size(),
            rc_a
        );


        /*
         * Primer B:
         * Forward hits at 15000, 35000
         * Reverse hits at 15450, 35450
         */
        reference.replace(
            15000,
            primer_b.size(),
            primer_b
        );

        reference.replace(
            15450,
            rc_b.size(),
            rc_b
        );

        reference.replace(
            35000,
            primer_b.size(),
            primer_b
        );

        reference.replace(
            35450,
            rc_b.size(),
            rc_b
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
         * Critical layout:
         *
         * same sequence appears:
         *
         * - in Primer1 slots
         * - in Primer2 slots
         * - on both sides of same intended pair
         * - in duplicated pair definitions
         *
         * Only TWO unique primer sequences exist.
         */
        const std::vector<
            MultiplexPrimerPairRequest
        > requests{
            {
                primer_a,
                primer_a,
                3,
                50,
                15000
            },

            {
                primer_a,
                primer_a,
                3,
                50,
                15000
            },

            {
                primer_a,
                primer_b,
                3,
                50,
                15000
            },

            {
                primer_b,
                primer_a,
                3,
                50,
                15000
            },

            {
                primer_a,
                primer_a,
                3,
                50,
                15000
            },

            {
                primer_b,
                primer_b,
                3,
                50,
                15000
            }
        };


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
                "Repeated-primer intended "
                "result-count mismatch."
            );
        }


        std::size_t intended_checks = 0;


        for (
            std::size_t i = 0;
            i < requests.size();
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
                    << "REPEATED_INTENDED_MISMATCH"
                    << '\t'
                    << "pair="
                    << i
                    << '\t'
                    << "expected="
                    << expected.amplicons.size()
                    << '\t'
                    << "observed="
                    << observed.amplicons.size()
                    << '\n';

                return 1;
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
                << "REPEATED_CROSS_MISMATCH"
                << '\t'
                << "expected="
                << expected_cross.size()
                << '\t'
                << "observed="
                << observed_cross.size()
                << '\n';

            return 1;
        }


        /*
         * There are 12 slots but only A and B.
         */
        if (
            result_v2
                .stats
                .total_primer_slots !=
            12
        ) {
            throw std::runtime_error(
                "Repeated-primer slot-count mismatch."
            );
        }


        if (
            result_v2
                .stats
                .unique_primer_queries !=
            2
        ) {
            throw std::runtime_error(
                "Repeated-primer deduplication failed."
            );
        }


        if (
            result_v2
                .stats
                .reused_primer_slots !=
            10
        ) {
            throw std::runtime_error(
                "Repeated-primer reuse-count mismatch."
            );
        }


        if (
            result_v2
                .global_cross_stats
                .forward_hits ==
            0
            ||
            result_v2
                .global_cross_stats
                .reverse_hits ==
            0
        ) {
            throw std::runtime_error(
                "Expected both forward and reverse "
                "global hits."
            );
        }


        if (
            result_v2
                .global_cross_stats
                .unique_products ==
            0
        ) {
            throw std::runtime_error(
                "Expected non-zero global products."
            );
        }


        std::cout
            << "pairs\t"
            << requests.size()
            << '\n';


        std::cout
            << "primer_slots\t"
            << result_v2
                .stats
                .total_primer_slots
            << '\n';


        std::cout
            << "unique_primers\t"
            << result_v2
                .stats
                .unique_primer_queries
            << '\n';


        std::cout
            << "reused_slots\t"
            << result_v2
                .stats
                .reused_primer_slots
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
