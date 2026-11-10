#pragma once

#include <cstddef>
#include <string_view>
#include <vector>

#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/search_strategy.hpp"
#include "primerpair/single_primer_search.hpp"

namespace primerpair {

struct HybridBatchedPrimerRequest {
    std::string_view primer{};
    std::size_t max_mismatches{3};
};

struct HybridBatchedPrimerResult {
    SearchStrategy strategy{
        SearchStrategy::DirectBranching
    };

    /*
     * true:
     *     IP-BWT batched candidate backend
     *
     * false:
     *     legacy BiFM DirectBranching backend
     */
    bool used_candidate_backend{false};

    std::vector<PrimerSearchHit> hits{};
};

class HybridBatchedPrimerSearchEngine {
public:
    HybridBatchedPrimerSearchEngine(
        const BatchedCandidateSearchEngine&
            candidate_engine,
        const SinglePrimerSearchEngine&
            legacy_engine
    ) noexcept;

    [[nodiscard]]
    std::vector<HybridBatchedPrimerResult>
    search(
        const std::vector<
            HybridBatchedPrimerRequest
        >& requests,
        std::size_t anchor_length = 12
    ) const;

private:
    const BatchedCandidateSearchEngine&
        candidate_engine_;

    const SinglePrimerSearchEngine&
        legacy_engine_;
};

}  // namespace primerpair
