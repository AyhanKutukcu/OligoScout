#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
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
#include <utility>
#include <vector>

namespace {


std::string make_reference(
    const std::size_t length
) {
    constexpr char alphabet[] = {
        'A',
        'C',
        'G',
        'T'
    };

    std::uint64_t state =
        0xD6E8FEB86659FD93ULL;

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
            alphabet[
                static_cast<std::size_t>(
                    (state >> 32U) &
                    3ULL
                )
            ]
        );
    }


    /*
     * Repeat-rich islands.
     */
    const std::string motif =
        "ACGTACGTACGTACGT"
        "AAAAAAAAAAAAAAAA"
        "CCCCCCCCCCCCCCCC"
        "GGGGGGGGGGGGGGGG"
        "TTTTTTTTTTTTTTTT"
        "ACACACACGTGTGTGT";

    for (
        std::size_t start = 4096;
        start + 768 < reference.size();
        start += 12000
    ) {
        for (
            std::size_t offset = 0;
            offset < 768;
            ++offset
        ) {
            reference[
                start +
                offset
            ] =
                motif[
                    offset %
                    motif.size()
                ];
        }
    }

    return reference;
}


char mutate(
    const char base
) {
    switch (base) {
        case 'A': return 'C';
        case 'C': return 'G';
        case 'G': return 'T';
        case 'T': return 'A';

        default:
            throw std::runtime_error(
                "Unexpected primer base."
            );
    }
}


void normalize(
    std::vector<
        primerpair::OrientedPrimerSearchHit
    >& hits
) {
    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const auto& lhs,
            const auto& rhs
        ) {
            if (
                lhs.position !=
                rhs.position
            ) {
                return
                    lhs.position <
                    rhs.position;
            }

            if (
                lhs.orientation !=
                rhs.orientation
            ) {
                return
                    static_cast<int>(
                        lhs.orientation
                    ) <
                    static_cast<int>(
                        rhs.orientation
                    );
            }

            if (
                lhs.mismatches !=
                rhs.mismatches
            ) {
                return
                    lhs.mismatches <
                    rhs.mismatches;
            }

            return
                lhs.mismatch_mask <
                rhs.mismatch_mask;
        }
    );

    hits.erase(
        std::unique(
            hits.begin(),
            hits.end()
        ),
        hits.end()
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
            hybrid(
                forward_hybrid,
                anchor_lookup,
                ipbwt,
                packed,
                legacy_strand
            );


        std::vector<std::string>
            owned_primers;

        std::vector<std::size_t>
            mismatch_budgets;

        owned_primers.reserve(
            1600
        );

        mismatch_budgets.reserve(
            1600
        );


        for (
            std::size_t i = 0;
            i < 300;
            ++i
        ) {
            const std::size_t position =
                (
                    i *
                    3571 +
                    211
                )
                %
                (
                    reference.size() -
                    primer_length
                );

            const std::string original =
                reference.substr(
                    position,
                    primer_length
                );


            /*
             * k=3:
             * usually candidate route.
             */
            owned_primers.push_back(
                original
            );

            mismatch_budgets.push_back(
                3
            );


            /*
             * Allowed mismatch outside forward
             * 3-prime suffix anchor.
             */
            std::string mutated_one =
                original;

            mutated_one.at(1) =
                mutate(
                    mutated_one.at(1)
                );

            owned_primers.push_back(
                std::move(
                    mutated_one
                )
            );

            mismatch_budgets.push_back(
                3
            );


            /*
             * k=1 always DirectBranching.
             */
            owned_primers.push_back(
                original
            );

            mismatch_budgets.push_back(
                1
            );


            /*
             * k=0 DirectBranching exact path.
             */
            owned_primers.push_back(
                original
            );

            mismatch_budgets.push_back(
                0
            );
        }


        /*
         * Explicit repeat-rich requests.
         */
        owned_primers.push_back(
            "TTTTTTTTTTTT"
            "ACGTACGTACGT"
        );

        mismatch_budgets.push_back(
            3
        );

        owned_primers.push_back(
            "AAAAAAAAAAAA"
            "ACACACACACAC"
        );

        mismatch_budgets.push_back(
            3
        );

        owned_primers.push_back(
            "CCCCCCCCCCCC"
            "GTGTGTGTGTGT"
        );

        mismatch_budgets.push_back(
            1
        );


        std::vector<
            HybridStrandAwarePrimerRequest
        > requests;

        requests.reserve(
            owned_primers.size()
        );

        for (
            std::size_t i = 0;
            i < owned_primers.size();
            ++i
        ) {
            requests.push_back(
                HybridStrandAwarePrimerRequest{
                    owned_primers.at(i),
                    mismatch_budgets.at(i)
                }
            );
        }


        const auto observed =
            hybrid.search(
                requests,
                anchor_length
            );

        if (
            observed.size() !=
            requests.size()
        ) {
            throw std::runtime_error(
                "Hybrid strand result-count mismatch."
            );
        }


        std::size_t checks = 0;

        std::size_t
            forward_candidate_routes = 0;

        std::size_t
            forward_branching_routes = 0;

        std::size_t
            reverse_candidate_routes = 0;

        std::size_t
            reverse_branching_routes = 0;


        for (
            std::size_t i = 0;
            i < requests.size();
            ++i
        ) {
            auto expected =
                legacy_strand.search(
                    requests.at(i).primer,
                    anchor_length,
                    requests.at(i)
                        .max_mismatches
                );

            auto expected_hits =
                expected.hits;

            auto observed_hits =
                observed.at(i).hits;

            normalize(
                expected_hits
            );

            normalize(
                observed_hits
            );


            if (
                expected
                    .forward_decision
                    .recommended_strategy
                !=
                observed.at(i)
                    .forward_strategy
            ) {
                throw std::runtime_error(
                    "Forward route mismatch."
                );
            }


            if (
                expected
                    .reverse_decision
                    .recommended_strategy
                !=
                observed.at(i)
                    .reverse_strategy
            ) {
                throw std::runtime_error(
                    "Reverse route mismatch."
                );
            }


            if (
                expected_hits !=
                observed_hits
            ) {
                std::cerr
                    << "STRAND_HIT_MISMATCH"
                    << '\t'
                    << "request="
                    << i
                    << '\t'
                    << "expected="
                    << expected_hits.size()
                    << '\t'
                    << "observed="
                    << observed_hits.size()
                    << '\n';

                throw std::runtime_error(
                    "Hybrid strand hit-set mismatch."
                );
            }


            if (
                observed.at(i)
                    .forward_candidate_backend
            ) {
                ++forward_candidate_routes;
            } else {
                ++forward_branching_routes;
            }


            if (
                observed.at(i)
                    .reverse_candidate_backend
            ) {
                ++reverse_candidate_routes;
            } else {
                ++reverse_branching_routes;
            }


            ++checks;
        }


        if (
            forward_candidate_routes == 0 ||
            forward_branching_routes == 0 ||
            reverse_candidate_routes == 0 ||
            reverse_branching_routes == 0
        ) {
            throw std::runtime_error(
                "Both-strand mixed routing was not exercised."
            );
        }


        std::cout
            << "requests\t"
            << requests.size()
            << '\n';

        std::cout
            << "forward_candidate_routes\t"
            << forward_candidate_routes
            << '\n';

        std::cout
            << "forward_branching_routes\t"
            << forward_branching_routes
            << '\n';

        std::cout
            << "reverse_candidate_routes\t"
            << reverse_candidate_routes
            << '\n';

        std::cout
            << "reverse_branching_routes\t"
            << reverse_branching_routes
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
