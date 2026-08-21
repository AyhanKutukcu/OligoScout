#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/search_strategy.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string
make_reference(
    const std::size_t length
) {
    constexpr char alphabet[] = {
        'A',
        'C',
        'G',
        'T'
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
    std::vector<primerpair::AnchorCandidateHit>& hits
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
}


void
compare_candidate_results(
    primerpair::AnchorCandidateSearchResult expected,
    primerpair::AnchorCandidateSearchResult observed
) {
    if (
        expected.primer_length !=
            observed.primer_length
        ||
        expected.anchor_length !=
            observed.anchor_length
        ||
        expected.max_mismatches !=
            observed.max_mismatches
        ||
        expected.anchor_occurrences !=
            observed.anchor_occurrences
        ||
        expected.candidates_verified !=
            observed.candidates_verified
    ) {
        throw std::runtime_error(
            "Candidate metadata mismatch."
        );
    }

    normalize_hits(
        expected.hits
    );

    normalize_hits(
        observed.hits
    );

    if (
        expected.hits.size() !=
        observed.hits.size()
    ) {
        throw std::runtime_error(
            "Candidate hit-count mismatch."
        );
    }

    for (
        std::size_t i = 0;
        i < expected.hits.size();
        ++i
    ) {
        if (
            expected.hits[i].position !=
                observed.hits[i].position
            ||
            expected.hits[i].mismatches !=
                observed.hits[i].mismatches
        ) {
            throw std::runtime_error(
                "Candidate hit mismatch."
            );
        }
    }
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
            batch_engine(
                ipbwt,
                anchor_lookup,
                verifier
            );


        std::vector<std::string>
            owned_primers;

        owned_primers.reserve(
            2048
        );


        /*
         * Build a mixed panel containing:
         *
         * - exact primers
         * - 1/2/3 prefix mismatch primers
         * - repetitive-anchor primers
         *
         * Exact 3' anchor is preserved in the
         * mutation cases.
         */
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
                        m *
                        3;


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
         * Explicit repeat-rich requests to force
         * non-candidate routing cases as well.
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


        std::vector<BatchedCandidateRequest>
            requests;

        requests.reserve(
            owned_primers.size()
        );


        for (
            const auto& primer :
            owned_primers
        ) {
            requests.push_back(
                BatchedCandidateRequest{
                    primer,
                    3
                }
            );
        }


        const auto observed =
            batch_engine.search(
                requests,
                anchor_length
            );


        if (
            observed.size() !=
            requests.size()
        ) {
            throw std::runtime_error(
                "Batched candidate result-count mismatch."
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
            const auto legacy_decision =
                SearchStrategyRouter(
                    bifm
                ).decide(
                    requests[i].primer,
                    anchor_length,
                    requests[i].max_mismatches
                );


            if (
                observed[i].decision.strategy !=
                legacy_decision.recommended_strategy
            ) {
                throw std::runtime_error(
                    "Routing strategy mismatch."
                );
            }


            if (
                observed[i].decision.occurrences !=
                legacy_decision
                    .difficulty_profile
                    .anchor_occurrences
            ) {
                throw std::runtime_error(
                    "Anchor occurrence mismatch."
                );
            }


            if (
                observed[i].decision.difficulty !=
                legacy_decision
                    .difficulty_profile
                    .difficulty
            ) {
                throw std::runtime_error(
                    "Difficulty mismatch."
                );
            }


            if (
                legacy_decision.recommended_strategy ==
                SearchStrategy::
                    AnchorCandidateVerification
            ) {
                ++candidate_routes;


                if (
                    !observed[i].candidate_executed
                ) {
                    throw std::runtime_error(
                        "Candidate route was not executed."
                    );
                }


                const auto expected =
                    verifier.search(
                        requests[i].primer,
                        anchor_length,
                        requests[i].max_mismatches
                    );


                compare_candidate_results(
                    expected,
                    observed[i].candidate_result
                );

            } else {
                ++branching_routes;


                if (
                    observed[i].candidate_executed
                ) {
                    throw std::runtime_error(
                        "DirectBranching request "
                        "incorrectly executed candidate backend."
                    );
                }
            }


            ++checks;
        }


        if (
            candidate_routes == 0
        ) {
            throw std::runtime_error(
                "Test produced no candidate-routed requests."
            );
        }


        if (
            branching_routes == 0
        ) {
            throw std::runtime_error(
                "Test produced no direct-branching requests."
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
