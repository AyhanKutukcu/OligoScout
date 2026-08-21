#include "primerpair/hybrid_strand_aware_primer_search.hpp"

#include <algorithm>
#include <bit>
#include <cctype>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace primerpair {

namespace {


std::uint64_t build_mismatch_mask(
    const PackedReference& reference,
    const std::uint64_t start,
    const std::string_view query,
    const bool reverse_to_original
) {
    if (query.size() > 64) {
        throw std::invalid_argument(
            "Mismatch mask supports primers up to 64 nt."
        );
    }

    if (start > reference.size()) {
        throw std::out_of_range(
            "Mismatch-mask start outside reference."
        );
    }

    if (
        query.size() >
        reference.size() - start
    ) {
        throw std::out_of_range(
            "Mismatch-mask query exceeds reference."
        );
    }

    std::uint64_t mask = 0;

    for (
        std::size_t query_position = 0;
        query_position < query.size();
        ++query_position
    ) {
        const char query_base =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        query.at(
                            query_position
                        )
                    )
                )
            );

        const char reference_base =
            reference.base_at(
                start +
                static_cast<std::uint64_t>(
                    query_position
                )
            );

        if (
            query_base ==
            reference_base
        ) {
            continue;
        }

        const std::size_t original_position =
            reverse_to_original
                ? query.size() -
                    1 -
                    query_position
                : query_position;

        mask |=
            (
                std::uint64_t{1}
                <<
                original_position
            );
    }

    return mask;
}


void validate_mismatch_count(
    const std::uint64_t mask,
    const std::size_t expected
) {
    const std::size_t observed =
        static_cast<std::size_t>(
            std::popcount(
                mask
            )
        );

    if (observed != expected) {
        throw std::logic_error(
            "Mismatch count and mismatch mask disagree."
        );
    }
}


void normalize_hits(
    std::vector<
        OrientedPrimerSearchHit
    >& hits
) {
    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const OrientedPrimerSearchHit& lhs,
            const OrientedPrimerSearchHit& rhs
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


HybridStrandAwarePrimerSearchEngine::
HybridStrandAwarePrimerSearchEngine(
    const HybridBatchedPrimerSearchEngine&
        forward_hybrid,
    const BatchedAnchorLookup&
        anchor_lookup,
    const IPBWTIndex&
        ipbwt,
    const PackedReference&
        reference,
    const StrandAwarePrimerSearchEngine&
        legacy_strand_engine
) noexcept
    : forward_hybrid_(
          forward_hybrid
      ),
      anchor_lookup_(
          anchor_lookup
      ),
      ipbwt_(
          ipbwt
      ),
      reference_(
          reference
      ),
      legacy_strand_engine_(
          legacy_strand_engine
      ) {
}


std::vector<OrientedPrimerSearchHit>
HybridStrandAwarePrimerSearchEngine::
verify_reverse_candidates(
    const std::string_view reverse_query,
    std::vector<std::uint64_t> positions,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    if (
        anchor_length == 0 ||
        anchor_length >
            reverse_query.size()
    ) {
        throw std::invalid_argument(
            "Invalid reverse anchor length."
        );
    }

    /*
     * Prefix anchor is exact.
     *
     * Only the suffix of reverse_query is allowed
     * to contain mismatches.
     */
    const std::size_t suffix_length =
        reverse_query.size() -
        anchor_length;

    const std::string_view suffix(
        reverse_query.data() +
            anchor_length,
        suffix_length
    );

    std::sort(
        positions.begin(),
        positions.end()
    );

    positions.erase(
        std::unique(
            positions.begin(),
            positions.end()
        ),
        positions.end()
    );

    std::vector<
        OrientedPrimerSearchHit
    > hits;

    hits.reserve(
        positions.size()
    );

    for (
        const std::uint64_t candidate_start :
        positions
    ) {
        if (
            candidate_start >
            reference_.size()
        ) {
            continue;
        }

        if (
            reverse_query.size() >
            reference_.size() -
                candidate_start
        ) {
            continue;
        }

        std::size_t mismatches = 0;

        if (suffix_length != 0) {
            const std::uint64_t
                suffix_start =
                    candidate_start +
                    static_cast<std::uint64_t>(
                        anchor_length
                    );

            mismatches =
                reference_
                    .bounded_hamming_distance(
                        suffix_start,
                        suffix,
                        max_mismatches
                    );

            if (
                mismatches >
                max_mismatches
            ) {
                continue;
            }
        }

        const std::uint64_t mismatch_mask =
            build_mismatch_mask(
                reference_,
                candidate_start,
                reverse_query,
                true
            );

        validate_mismatch_count(
            mismatch_mask,
            mismatches
        );

        hits.push_back(
            OrientedPrimerSearchHit{
                candidate_start,
                mismatches,
                PrimerOrientation::Reverse,
                mismatch_mask
            }
        );
    }

    normalize_hits(
        hits
    );

    return hits;
}


std::vector<HybridStrandAwarePrimerResult>
HybridStrandAwarePrimerSearchEngine::search(
    const std::vector<
        HybridStrandAwarePrimerRequest
    >& requests,
    const std::size_t anchor_length
) const {
    if (requests.empty()) {
        return {};
    }

    /*
     * ------------------------------------------------
     * Forward panel.
     * ------------------------------------------------
     */

    std::vector<
        HybridBatchedPrimerRequest
    > forward_requests;

    forward_requests.reserve(
        requests.size()
    );

    for (
        const auto& request :
        requests
    ) {
        forward_requests.push_back(
            HybridBatchedPrimerRequest{
                request.primer,
                request.max_mismatches
            }
        );
    }

    const auto forward =
        forward_hybrid_.search(
            forward_requests,
            anchor_length
        );

    if (
        forward.size() !=
        requests.size()
    ) {
        throw std::logic_error(
            "Forward hybrid result-count mismatch."
        );
    }


    /*
     * ------------------------------------------------
     * Reverse queries.
     *
     * reserve() is important because anchor_requests
     * store string_views into these strings.
     * ------------------------------------------------
     */

    std::vector<std::string>
        reverse_queries;

    reverse_queries.reserve(
        requests.size()
    );

    for (
        const auto& request :
        requests
    ) {
        reverse_queries.push_back(
            reverse_complement(
                request.primer
            )
        );
    }


    /*
     * ------------------------------------------------
     * Batched biological 3-prime PREFIX anchors.
     * ------------------------------------------------
     */

    std::vector<
        BatchedAnchorRequest
    > reverse_anchor_requests;

    reverse_anchor_requests.reserve(
        requests.size()
    );

    for (
        std::size_t i = 0;
        i < requests.size();
        ++i
    ) {
        reverse_anchor_requests.push_back(
            BatchedAnchorRequest{
                reverse_queries.at(i),
                AnchorPlacement::Prefix,
                requests.at(i)
                    .max_mismatches
            }
        );
    }

    const auto reverse_decisions =
        anchor_lookup_.lookup(
            reverse_anchor_requests,
            anchor_length
        );

    if (
        reverse_decisions.size() !=
        requests.size()
    ) {
        throw std::logic_error(
            "Reverse anchor decision-count mismatch."
        );
    }


    /*
     * ------------------------------------------------
     * Merge forward + reverse.
     * ------------------------------------------------
     */

    std::vector<
        HybridStrandAwarePrimerResult
    > results;

    results.reserve(
        requests.size()
    );

    for (
        std::size_t i = 0;
        i < requests.size();
        ++i
    ) {
        HybridStrandAwarePrimerResult
            result;

        result.primer_length =
            requests.at(i)
                .primer
                .size();

        result.forward_strategy =
            forward.at(i)
                .strategy;

        result.forward_candidate_backend =
            forward.at(i)
                .used_candidate_backend;

        result.reverse_strategy =
            reverse_decisions.at(i)
                .strategy;


        /*
         * Forward hits -> oriented hits.
         */
        result.hits.reserve(
            forward.at(i)
                .hits
                .size()
        );

        for (
            const auto& hit :
            forward.at(i).hits
        ) {
            const std::uint64_t mismatch_mask =
                build_mismatch_mask(
                    reference_,
                    hit.position,
                    requests.at(i).primer,
                    false
                );

            validate_mismatch_count(
                mismatch_mask,
                hit.mismatches
            );

            result.hits.push_back(
                OrientedPrimerSearchHit{
                    hit.position,
                    hit.mismatches,
                    PrimerOrientation::Forward,
                    mismatch_mask
                }
            );
        }


        /*
         * Reverse route.
         */
        std::vector<
            OrientedPrimerSearchHit
        > reverse_hits;

        switch (
            reverse_decisions.at(i)
                .strategy
        ) {
            case SearchStrategy::
                AnchorCandidateVerification: {

                auto positions =
                    ipbwt_.locate(
                        reverse_decisions
                            .at(i)
                            .interval
                    );

                reverse_hits =
                    verify_reverse_candidates(
                        reverse_queries.at(i),
                        std::move(
                            positions
                        ),
                        anchor_length,
                        requests.at(i)
                            .max_mismatches
                    );

                result.reverse_candidate_backend =
                    true;

                break;
            }


            case SearchStrategy::
                DirectBranching: {

                /*
                 * Minority fallback.
                 *
                 * Uses the already validated
                 * biological reverse-prefix +
                 * extend_right implementation.
                 */
                const auto legacy =
                    legacy_strand_engine_.search(
                        requests.at(i).primer,
                        anchor_length,
                        requests.at(i)
                            .max_mismatches
                    );

                if (
                    legacy
                        .reverse_decision
                        .recommended_strategy
                    !=
                    result.reverse_strategy
                ) {
                    throw std::logic_error(
                        "Hybrid/legacy reverse routing disagreement."
                    );
                }

                for (
                    const auto& hit :
                    legacy.hits
                ) {
                    if (
                        hit.orientation ==
                        PrimerOrientation::Reverse
                    ) {
                        reverse_hits.push_back(
                            hit
                        );
                    }
                }

                result.reverse_candidate_backend =
                    false;

                break;
            }


            case SearchStrategy::
                SplitSeedCandidate:

                throw std::logic_error(
                    "SplitSeedCandidate backend "
                    "is not implemented."
                );
        }


        result.hits.insert(
            result.hits.end(),
            reverse_hits.begin(),
            reverse_hits.end()
        );

        normalize_hits(
            result.hits
        );

        results.push_back(
            std::move(
                result
            )
        );
    }

    return results;
}

}  // namespace primerpair
