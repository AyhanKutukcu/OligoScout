#include "primerpair/hybrid_batched_primer_search.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace primerpair {

namespace {

void
normalize_hits(
    std::vector<PrimerSearchHit>& hits
) {
    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const PrimerSearchHit& lhs,
            const PrimerSearchHit& rhs
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


HybridBatchedPrimerSearchEngine::
HybridBatchedPrimerSearchEngine(
    const BatchedCandidateSearchEngine&
        candidate_engine,
    const SinglePrimerSearchEngine&
        legacy_engine
) noexcept
    : candidate_engine_(
          candidate_engine
      ),
      legacy_engine_(
          legacy_engine
      ) {
}


std::vector<HybridBatchedPrimerResult>
HybridBatchedPrimerSearchEngine::search(
    const std::vector<
        HybridBatchedPrimerRequest
    >& requests,
    const std::size_t anchor_length
) const {
    std::vector<BatchedCandidateRequest>
        candidate_requests;

    candidate_requests.reserve(
        requests.size()
    );


    for (
        const auto& request :
        requests
    ) {
        candidate_requests.push_back(
            BatchedCandidateRequest{
                request.primer,
                request.max_mismatches
            }
        );
    }


    /*
     * IMPORTANT:
     *
     * This performs one batched IP-BWT anchor
     * lookup/routing pass for the entire panel.
     */
    const auto batched =
        candidate_engine_.search(
            candidate_requests,
            anchor_length
        );


    if (
        batched.size() !=
        requests.size()
    ) {
        throw std::logic_error(
            "Hybrid batched result-count mismatch."
        );
    }


    std::vector<HybridBatchedPrimerResult>
        results;

    results.reserve(
        requests.size()
    );


    for (
        std::size_t i = 0;
        i < requests.size();
        ++i
    ) {
        HybridBatchedPrimerResult result;

        result.strategy =
            batched[i]
                .decision
                .strategy;


        if (
            batched[i].candidate_executed
        ) {
            result.used_candidate_backend =
                true;


            result.hits.reserve(
                batched[i]
                    .candidate_result
                    .hits
                    .size()
            );


            for (
                const auto& hit :
                batched[i]
                    .candidate_result
                    .hits
            ) {
                result.hits.push_back(
                    PrimerSearchHit{
                        hit.position,
                        hit.mismatches
                    }
                );
            }


            normalize_hits(
                result.hits
            );

        } else {
            /*
             * First production hybrid version:
             *
             * DirectBranching requests fall back
             * to the already validated legacy
             * SinglePrimerSearchEngine.
             *
             * This intentionally recomputes the
             * BiFM anchor for this minority path.
             * We optimize that only if profiling
             * later shows it matters.
             */
            const auto legacy =
                legacy_engine_.search(
                    requests[i].primer,
                    anchor_length,
                    requests[i]
                        .max_mismatches
                );


            if (
                legacy
                    .decision
                    .recommended_strategy
                !=
                result.strategy
            ) {
                throw std::logic_error(
                    "Hybrid/legacy routing "
                    "disagreement."
                );
            }


            if (
                result.strategy !=
                SearchStrategy::
                    DirectBranching
            ) {
                throw std::logic_error(
                    "Non-candidate request did "
                    "not route to DirectBranching."
                );
            }


            result.hits =
                legacy.hits;


            result.used_candidate_backend =
                false;
        }


        results.push_back(
            std::move(
                result
            )
        );
    }


    return results;
}

}  // namespace primerpair
