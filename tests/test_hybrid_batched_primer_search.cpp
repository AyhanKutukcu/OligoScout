#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

std::string
make_reference(
    const std::size_t length
) {
    constexpr char alphabet[] = {
        'A', 'C', 'G', 'T'
    };


    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;


    std::string reference;
    reference.reserve(length);


    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL
            +
            1442695040888963407ULL;


        reference.push_back(
            alphabet[
                static_cast<std::size_t>(
                    (state >> 32U) & 3ULL
                )
            ]
        );
    }


    const std::string motif =
        "ACGTACGTACGTACGT"
        "AAAAAAAACCCCCCCC"
        "GGGGGGGGTTTTTTTT"
        "ACACACACGTGTGTGT";


    for (
        std::size_t start = 4096;
        start + 512 < reference.size();
        start += 8192
    ) {
        for (
            std::size_t offset = 0;
            offset < 512;
            ++offset
        ) {
            reference[start + offset] =
                motif[
                    offset %
                    motif.size()
                ];
        }
    }


    return reference;
}


char
mutated_base(
    const char base
) {
    return
        base == 'A'
            ? 'T'
            : 'A';
}


void
normalize_hits(
    std::vector<primerpair::PrimerSearchHit>& hits
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

            return
                lhs.mismatches <
                rhs.mismatches;
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


int
main() {
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
            legacy_engine(
                bifm,
                packed
            );


        HybridBatchedPrimerSearchEngine
            hybrid_engine(
                candidate_engine,
                legacy_engine
            );


        std::vector<std::string>
            owned_primers;

        owned_primers.reserve(
            2048
        );


        for (
            std::size_t i = 0;
            i < 400;
            ++i
        ) {
            const std::size_t position =
                (
                    i *
                    3571
                    +
                    123
                )
                %
                (
                    reference.size()
                    -
                    primer_length
                );


            const std::string original =
                reference.substr(
                    position,
                    primer_length
                );


            owned_primers.push_back(
                original
            );


            for (
                std::size_t introduced = 1;
                introduced <= 3;
                ++introduced
            ) {
                std::string primer =
                    original;


                for (
                    std::size_t m = 0;
                    m < introduced;
                    ++m
                ) {
                    const std::size_t p =
                        1 +
                        m * 3;


                    primer[p] =
                        mutated_base(
                            primer[p]
                        );
                }


                owned_primers.push_back(
                    std::move(
                        primer
                    )
                );
            }
        }


        /*
         * Explicit repeat-rich primers so the
         * panel contains DirectBranching routes.
         */
        owned_primers.push_back(
            "TTTTTTTTTTTT"
            "ACGTACGTACGT"
        );

        owned_primers.push_back(
            "AAAAAAAAAAAA"
            "ACACACACACAC"
        );

        owned_primers.push_back(
            "CCCCCCCCCCCC"
            "GTGTGTGTGTGT"
        );


        std::vector<
            HybridBatchedPrimerRequest
        > requests;

        requests.reserve(
            owned_primers.size()
        );


        for (
            const auto& primer :
            owned_primers
        ) {
            requests.push_back(
                HybridBatchedPrimerRequest{
                    primer,
                    3
                }
            );
        }


        const auto observed =
            hybrid_engine.search(
                requests,
                anchor_length
            );


        if (
            observed.size() !=
            requests.size()
        ) {
            throw std::runtime_error(
                "Hybrid result-count mismatch."
            );
        }


        std::size_t checks = 0;
        std::size_t candidate_routes = 0;
        std::size_t branching_routes = 0;


        for (
            std::size_t i = 0;
            i < requests.size();
            ++i
        ) {
            auto expected =
                legacy_engine.search(
                    requests[i].primer,
                    anchor_length,
                    requests[i]
                        .max_mismatches
                );


            auto expected_hits =
                expected.hits;

            auto observed_hits =
                observed[i].hits;


            normalize_hits(
                expected_hits
            );

            normalize_hits(
                observed_hits
            );


            if (
                expected.decision
                    .recommended_strategy
                !=
                observed[i].strategy
            ) {
                throw std::runtime_error(
                    "Hybrid strategy mismatch."
                );
            }


            if (
                expected_hits !=
                observed_hits
            ) {
                std::cerr
                    << "HIT_MISMATCH"
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
                    "Hybrid final hit-set mismatch."
                );
            }


            if (
                observed[i]
                    .used_candidate_backend
            ) {
                ++candidate_routes;


                if (
                    observed[i].strategy !=
                    SearchStrategy::
                        AnchorCandidateVerification
                ) {
                    throw std::runtime_error(
                        "Candidate backend used "
                        "for wrong strategy."
                    );
                }

            } else {
                ++branching_routes;


                if (
                    observed[i].strategy !=
                    SearchStrategy::
                        DirectBranching
                ) {
                    throw std::runtime_error(
                        "Legacy backend used "
                        "for wrong strategy."
                    );
                }
            }


            ++checks;
        }


        if (
            candidate_routes == 0 ||
            branching_routes == 0
        ) {
            throw std::runtime_error(
                "Mixed routing was not exercised."
            );
        }


        std::cout
            << "requests\t"
            << requests.size()
            << '\n';


        std::cout
            << "candidate_routes\t"
            << candidate_routes
            << '\n';


        std::cout
            << "branching_routes\t"
            << branching_routes
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
