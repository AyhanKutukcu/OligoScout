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
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {


using Clock =
    std::chrono::steady_clock;


std::string make_reference(
    const std::size_t length
) {
    constexpr char bases[] = {
        'A', 'C', 'G', 'T'
    };

    std::uint64_t state =
        0xF1357AEA2E62A9C5ULL;

    std::string reference;

    reference.reserve(
        length
    );


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


double ns_between(
    const Clock::time_point begin,
    const Clock::time_point end
) {
    return
        static_cast<double>(
            std::chrono::duration_cast<
                std::chrono::nanoseconds
            >(
                end - begin
            ).count()
        );
}


}  // namespace


int main(
    int argc,
    char** argv
) {
    try {
        using namespace primerpair;


        std::size_t repetitions = 5;

        if (argc >= 2) {
            repetitions =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[1]
                    )
                );
        }


        if (repetitions == 0) {
            throw std::invalid_argument(
                "Repetitions must be > 0."
            );
        }


        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            primer_length = 24;

        constexpr std::uint64_t
            min_amplicon = 50;

        constexpr std::uint64_t
            max_amplicon = 1500;


        /*
         * Large enough for 128 separated pairs.
         */
        const std::string reference =
            make_reference(
                500000
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


        const std::vector<std::size_t>
            panel_sizes{
                8,
                16,
                32,
                64,
                128
            };


        std::uint64_t sink = 0;


        std::cout
            << "pairs"
            << '\t'
            << "primer_slots"
            << '\t'
            << "unique_primers"
            << '\t'
            << "logical_cross_requests"
            << '\t'
            << "global_forward_hits"
            << '\t'
            << "global_reverse_hits"
            << '\t'
            << "global_window_candidates"
            << '\t'
            << "global_unique_products"
            << '\t'
            << "cross_amplicons"
            << '\t'
            << "v1_ns"
            << '\t'
            << "v2_ns"
            << '\t'
            << "v1_ns_per_pair"
            << '\t'
            << "v2_ns_per_pair"
            << '\t'
            << "speedup"
            << '\n';


        for (
            const std::size_t pair_count :
            panel_sizes
        ) {
            std::vector<std::string>
                primer1_storage;

            std::vector<std::string>
                primer2_storage;


            primer1_storage.reserve(
                pair_count
            );

            primer2_storage.reserve(
                pair_count
            );


            /*
             * All primer sequences are expected to be
             * unique.
             *
             * Pair starts are separated by 500 nt.
             *
             * Within each pair:
             *
             * forward ----350 nt---- reverse
             *
             * max amplicon 1500 means only local
             * genomic neighbours can form products.
             */
            for (
                std::size_t i = 0;
                i < pair_count;
                ++i
            ) {
                const std::size_t left =
                    20000 +
                    i * 500;


                const std::size_t right =
                    left +
                    350 +
                    (
                        i %
                        31
                    );


                primer1_storage.push_back(
                    reference.substr(
                        left,
                        primer_length
                    )
                );


                primer2_storage.push_back(
                    reverse_complement(
                        reference.substr(
                            right,
                            primer_length
                        )
                    )
                );
            }


            std::vector<
                MultiplexPrimerPairRequest
            > requests;


            requests.reserve(
                pair_count
            );


            for (
                std::size_t i = 0;
                i < pair_count;
                ++i
            ) {
                requests.push_back(
                    MultiplexPrimerPairRequest{
                        primer1_storage.at(i),
                        primer2_storage.at(i),

                        3,

                        min_amplicon,
                        max_amplicon
                    }
                );
            }


            /*
             * ----------------------------------------
             * Correctness before timing.
             * ----------------------------------------
             */
            const auto check_v1 =
                v1.search(
                    requests,
                    anchor_length,
                    true,
                    min_amplicon,
                    max_amplicon
                );


            const auto check_v2 =
                v2.search(
                    requests,
                    anchor_length,
                    true,
                    min_amplicon,
                    max_amplicon
                );


            if (
                check_v1.intended_pairs.size() !=
                check_v2.intended_pairs.size()
            ) {
                throw std::runtime_error(
                    "Scaling intended result-count mismatch."
                );
            }


            for (
                std::size_t i = 0;
                i < pair_count;
                ++i
            ) {
                if (
                    check_v1
                        .intended_pairs
                        .at(i)
                        .amplicons
                    !=
                    check_v2
                        .intended_pairs
                        .at(i)
                        .amplicons
                ) {
                    throw std::runtime_error(
                        "Scaling intended equivalence failure."
                    );
                }
            }


            auto cross_v1 =
                check_v1.cross_amplicons;


            auto cross_v2 =
                check_v2.cross_amplicons;


            normalize(
                cross_v1
            );

            normalize(
                cross_v2
            );


            if (
                cross_v1 !=
                cross_v2
            ) {
                std::cerr
                    << "SCALING_CROSS_MISMATCH"
                    << '\t'
                    << "pairs="
                    << pair_count
                    << '\t'
                    << "v1="
                    << cross_v1.size()
                    << '\t'
                    << "v2="
                    << cross_v2.size()
                    << '\n';

                return 1;
            }


            /*
             * This benchmark is meant to be
             * all-unique.
             */
            if (
                check_v2
                    .stats
                    .unique_primer_queries
                !=
                pair_count * 2
            ) {
                std::cerr
                    << "pairs="
                    << pair_count
                    << '\t'
                    << "slots="
                    << pair_count * 2
                    << '\t'
                    << "unique="
                    << check_v2
                        .stats
                        .unique_primer_queries
                    << '\n';

                throw std::runtime_error(
                    "Scaling panel unexpectedly "
                    "contains duplicate primers."
                );
            }


            /*
             * Warm-up.
             */
            {
                const auto warm_v1 =
                    v1.search(
                        requests,
                        anchor_length,
                        true,
                        min_amplicon,
                        max_amplicon
                    );

                const auto warm_v2 =
                    v2.search(
                        requests,
                        anchor_length,
                        true,
                        min_amplicon,
                        max_amplicon
                    );

                sink +=
                    warm_v1
                        .cross_amplicons
                        .size();

                sink +=
                    warm_v2
                        .cross_amplicons
                        .size();
            }


            const auto v1_begin =
                Clock::now();


            for (
                std::size_t rep = 0;
                rep < repetitions;
                ++rep
            ) {
                const auto result =
                    v1.search(
                        requests,
                        anchor_length,
                        true,
                        min_amplicon,
                        max_amplicon
                    );

                sink +=
                    result
                        .cross_amplicons
                        .size();
            }


            const auto v1_end =
                Clock::now();


            const auto v2_begin =
                Clock::now();


            for (
                std::size_t rep = 0;
                rep < repetitions;
                ++rep
            ) {
                const auto result =
                    v2.search(
                        requests,
                        anchor_length,
                        true,
                        min_amplicon,
                        max_amplicon
                    );

                sink +=
                    result
                        .cross_amplicons
                        .size();
            }


            const auto v2_end =
                Clock::now();


            const double v1_ns =
                ns_between(
                    v1_begin,
                    v1_end
                )
                /
                repetitions;


            const double v2_ns =
                ns_between(
                    v2_begin,
                    v2_end
                )
                /
                repetitions;


            const std::size_t
                logical_cross_requests =
                    (
                        pair_count *
                        (
                            pair_count - 1
                        )
                        /
                        2
                    )
                    *
                    4;


            std::cout
                << pair_count
                << '\t'
                << pair_count * 2
                << '\t'
                << check_v2
                    .stats
                    .unique_primer_queries
                << '\t'
                << logical_cross_requests
                << '\t'
                << check_v2
                    .global_cross_stats
                    .forward_hits
                << '\t'
                << check_v2
                    .global_cross_stats
                    .reverse_hits
                << '\t'
                << check_v2
                    .global_cross_stats
                    .window_candidates
                << '\t'
                << check_v2
                    .global_cross_stats
                    .unique_products
                << '\t'
                << cross_v2.size()
                << '\t'
                << v1_ns
                << '\t'
                << v2_ns
                << '\t'
                << (
                    v1_ns /
                    pair_count
                )
                << '\t'
                << (
                    v2_ns /
                    pair_count
                )
                << '\t'
                << (
                    v1_ns /
                    v2_ns
                )
                << '\n';
        }


        std::cout
            << "VERIFY_ALL_EQUIVALENT\tYES\n";


        std::cerr
            << "sink\t"
            << sink
            << '\n';


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
