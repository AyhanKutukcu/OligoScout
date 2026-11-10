#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/multiplex_primer_search.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {


std::string make_reference(
    const std::size_t length
) {
    constexpr char bases[] = {
        'A', 'C', 'G', 'T'
    };

    std::uint64_t state =
        0xBF58476D1CE4E5B9ULL;

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


}  // namespace


int main() {
    try {
        using namespace primerpair;

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            primer_length = 24;

        constexpr std::size_t
            pair_count = 60;


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
            multiplex(
                hybrid_strand
            );

        PrimerPairSearchEngine
            legacy_pair(
                bifm,
                packed
            );


        /*
         * Shared primer appears in every pair.
         */
        const std::string shared_primer =
            reference.substr(
                1000,
                primer_length
            );


        std::vector<std::string>
            right_storage;

        right_storage.reserve(
            pair_count
        );


        for (
            std::size_t i = 0;
            i < pair_count;
            ++i
        ) {
            const std::size_t position =
                1300 +
                i * 1700;

            right_storage.push_back(
                reverse_complement(
                    reference.substr(
                        position,
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
                    shared_primer,
                    right_storage.at(i),
                    3,
                    50,
                    3000
                }
            );
        }


        const auto observed =
            multiplex.search(
                requests,
                anchor_length,
                false
            );


        if (
            observed.intended_pairs.size() !=
            requests.size()
        ) {
            throw std::runtime_error(
                "Multiplex intended result-count mismatch."
            );
        }


        std::size_t checks = 0;


        for (
            std::size_t i = 0;
            i < requests.size();
            ++i
        ) {
            const auto expected =
                legacy_pair.search(
                    requests.at(i).primer1,
                    requests.at(i).primer2,

                    anchor_length,
                    3,
                    50,
                    3000
                );


            if (
                expected.amplicons !=
                observed
                    .intended_pairs
                    .at(i)
                    .amplicons
            ) {
                std::cerr
                    << "MULTIPLEX_INTENDED_MISMATCH"
                    << '\t'
                    << "pair="
                    << i
                    << '\n';

                return 1;
            }


            ++checks;
        }


        if (
            observed
                .stats
                .total_primer_slots
            !=
            pair_count * 2
        ) {
            throw std::runtime_error(
                "Total primer-slot count mismatch."
            );
        }


        /*
         * shared primer + 60 distinct rights
         */
        if (
            observed
                .stats
                .unique_primer_queries
            !=
            pair_count + 1
        ) {
            std::cerr
                << "unique="
                << observed
                    .stats
                    .unique_primer_queries
                << '\n';

            throw std::runtime_error(
                "Unique-primer cache did not deduplicate correctly."
            );
        }


        if (
            observed
                .stats
                .reused_primer_slots
            !=
            pair_count - 1
        ) {
            throw std::runtime_error(
                "Primer-search reuse count mismatch."
            );
        }


        std::cout
            << "pairs\t"
            << pair_count
            << '\n';

        std::cout
            << "primer_slots\t"
            << observed.stats.total_primer_slots
            << '\n';

        std::cout
            << "unique_primer_queries\t"
            << observed.stats.unique_primer_queries
            << '\n';

        std::cout
            << "reused_primer_slots\t"
            << observed.stats.reused_primer_slots
            << '\n';

        std::cout
            << "checks\t"
            << checks
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
