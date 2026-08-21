#include "primerpair/batched_candidate_search.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace primerpair {

BatchedCandidateSearchEngine::
BatchedCandidateSearchEngine(
    const IPBWTIndex& ipbwt,
    const BatchedAnchorLookup& anchor_lookup,
    const AnchorCandidateSearcher& verifier
) noexcept
    : ipbwt_(ipbwt),
      anchor_lookup_(anchor_lookup),
      verifier_(verifier) {
}


std::vector<BatchedCandidateResult>
BatchedCandidateSearchEngine::search(
    const std::vector<BatchedCandidateRequest>& requests,
    const std::size_t anchor_length
) const {
    std::vector<BatchedAnchorRequest>
        anchor_requests;

    anchor_requests.reserve(
        requests.size()
    );


    for (
        const auto& request :
        requests
    ) {
        anchor_requests.push_back(
            BatchedAnchorRequest{
                request.primer,
                AnchorPlacement::Suffix,
                request.max_mismatches
            }
        );
    }


    const auto decisions =
        anchor_lookup_.lookup(
            anchor_requests,
            anchor_length
        );


    if (
        decisions.size() !=
        requests.size()
    ) {
        throw std::logic_error(
            "Batched anchor decision count mismatch."
        );
    }


    std::vector<BatchedCandidateResult>
        results;

    results.reserve(
        requests.size()
    );


    for (
        std::size_t i = 0;
        i < requests.size();
        ++i
    ) {
        BatchedCandidateResult result;

        result.decision =
            decisions[i];


        if (
            decisions[i].strategy ==
            SearchStrategy::
                AnchorCandidateVerification
        ) {
            auto positions =
                ipbwt_.locate(
                    decisions[i].interval
                );


            result.candidate_result =
                verifier_
                    .verify_from_anchor_positions(
                        requests[i].primer,
                        std::move(
                            positions
                        ),
                        decisions[i].occurrences,
                        anchor_length,
                        requests[i].max_mismatches
                    );


            result.candidate_executed =
                true;
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
