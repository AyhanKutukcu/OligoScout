#pragma once

#include <cstddef>
#include <string_view>
#include <vector>

#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/ip_bwt_index.hpp"

namespace primerpair {

struct BatchedCandidateRequest {
    std::string_view primer{};
    std::size_t max_mismatches{3};
};

struct BatchedCandidateResult {
    BatchedAnchorDecision decision{};

    /*
     * True only when routing selected
     * AnchorCandidateVerification.
     *
     * DirectBranching requests are intentionally
     * left for the existing BiFM backend.
     */
    bool candidate_executed{false};

    AnchorCandidateSearchResult candidate_result{};
};

class BatchedCandidateSearchEngine {
public:
    BatchedCandidateSearchEngine(
        const IPBWTIndex& ipbwt,
        const BatchedAnchorLookup& anchor_lookup,
        const AnchorCandidateSearcher& verifier
    ) noexcept;

    [[nodiscard]]
    std::vector<BatchedCandidateResult>
    search(
        const std::vector<BatchedCandidateRequest>& requests,
        std::size_t anchor_length = 12
    ) const;

private:
    const IPBWTIndex& ipbwt_;
    const BatchedAnchorLookup& anchor_lookup_;
    const AnchorCandidateSearcher& verifier_;
};

}  // namespace primerpair
