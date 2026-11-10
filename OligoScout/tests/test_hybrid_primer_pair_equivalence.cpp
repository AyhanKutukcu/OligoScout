#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/hybrid_primer_pair_search.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
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
        'A',
        'C',
        'G',
        'T'
    };

    std::uint64_t state =
        0x94D049BB133111EBULL;

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


    /*
     * Repeat-rich bölgeler:
     * routing testinin hem candidate hem branching
     * yollarını görme olasılığını artırır.
     */
    const std::string motif =
        "ACGTACGTACGTACGT"
        "AAAAAAAAAAAAAAAA"
        "CCCCCCCCCCCCCCCC"
        "GGGGGGGGGGGGGGGG"
        "TTTTTTTTTTTTTTTT";


    for (
        std::size_t start = 5000;
        start + 512 <
            reference.size();
        start += 13000
    ) {
        for (
            std::size_t offset = 0;
            offset < 512;
            ++offset
        ) {
            reference[
                start + offset
            ] =
                motif[
                    offset %
                    motif.size()
                ];
        }
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


        const std::string reference =
            make_reference(
                200000
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


        HybridPrimerPairSearchEngine
            hybrid_pair(
                hybrid_strand
            );


        PrimerPairSearchEngine
            legacy_pair(
                bifm,
                packed
            );


        std::vector<std::string>
            primer1_storage;

        std::vector<std::string>
            primer2_storage;

        primer1_storage.reserve(
            500
        );

        primer2_storage.reserve(
            500
        );


        for (
            std::size_t i = 0;
            i < 400;
            ++i
        ) {
            const std::size_t p1 =
                (
                    i *
                    3571 +
                    127
                )
                %
                (
                    reference.size() -
                    2000
                );


            const std::size_t distance =
                80 +
                (
                    i %
                    850
                );


            const std::size_t p2 =
                p1 +
                distance;


            primer1_storage.push_back(
                reference.substr(
                    p1,
                    primer_length
                )
            );


            primer2_storage.push_back(
                reverse_complement(
                    reference.substr(
                        p2,
                        primer_length
                    )
                )
            );
        }


        std::vector<
            HybridPrimerPairRequest
        > requests;

        requests.reserve(
            primer1_storage.size()
        );


        for (
            std::size_t i = 0;
            i < primer1_storage.size();
            ++i
        ) {
            requests.push_back(
                HybridPrimerPairRequest{
                    primer1_storage.at(i),
                    primer2_storage.at(i),
                    3,
                    50,
                    1000
                }
            );
        }


        const auto observed =
            hybrid_pair.search(
                requests,
                anchor_length
            );


        if (
            observed.size() !=
            requests.size()
        ) {
            throw std::runtime_error(
                "Hybrid pair result-count mismatch."
            );
        }


        std::size_t checks = 0;

        std::size_t
            forward_candidate = 0;

        std::size_t
            reverse_candidate = 0;


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

                    requests.at(i)
                        .max_mismatches,

                    requests.at(i)
                        .min_amplicon_length,

                    requests.at(i)
                        .max_amplicon_length
                );


            const auto& actual =
                observed.at(i)
                    .pair_result;


            if (
                expected.amplicons !=
                actual.amplicons
            ) {
                std::cerr
                    << "HYBRID_PAIR_MISMATCH"
                    << '\t'
                    << "request="
                    << i
                    << '\t'
                    << "expected="
                    << expected
                        .amplicons
                        .size()
                    << '\t'
                    << "actual="
                    << actual
                        .amplicons
                        .size()
                    << '\n';

                return 1;
            }


            if (
                expected
                    .primer1_single_hit_count
                !=
                actual
                    .primer1_single_hit_count
                ||
                expected
                    .primer2_single_hit_count
                !=
                actual
                    .primer2_single_hit_count
            ) {
                throw std::runtime_error(
                    "Hybrid pair single-hit "
                    "count mismatch."
                );
            }


            if (
                observed.at(i)
                    .primer1_search
                    .forward_candidate_backend
            ) {
                ++forward_candidate;
            }


            if (
                observed.at(i)
                    .primer2_search
                    .reverse_candidate_backend
            ) {
                ++reverse_candidate;
            }


            ++checks;
        }


        std::cout
            << "requests\t"
            << requests.size()
            << '\n';


        std::cout
            << "checks\t"
            << checks
            << '\n';


        std::cout
            << "primer1_forward_candidate\t"
            << forward_candidate
            << '\n';


        std::cout
            << "primer2_reverse_candidate\t"
            << reverse_candidate
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
