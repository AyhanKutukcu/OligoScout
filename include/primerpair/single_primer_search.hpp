#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_strategy.hpp"

namespace primerpair {

struct PrimerSearchHit {
    std::uint64_t position{0};
    std::size_t mismatches{0};

    bool operator==(
        const PrimerSearchHit&
    ) const = default;
};

struct SinglePrimerSearchResult {
    SearchStrategyDecision decision{};

    std::vector<PrimerSearchHit> hits{};

    [[nodiscard]]
    bool empty() const noexcept {
        return hits.empty();
    }

    [[nodiscard]]
    std::size_t hit_count() const noexcept {
        return hits.size();
    }
};

class SinglePrimerSearchEngine {
public:
    SinglePrimerSearchEngine(
        const BidirectionalFMIndex& index,
        const PackedReference& reference,
        SearchDifficultyThresholds thresholds = {}
    );

    [[nodiscard]]
    SinglePrimerSearchResult search(
        std::string_view primer,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 3
    ) const;

    [[nodiscard]]
    const SearchStrategyRouter&
    router() const noexcept {
        return router_;
    }

private:
    const BidirectionalFMIndex& index_;
    const PackedReference& reference_;

    SearchStrategyRouter router_;

    /*
     * Executes STRICT 5-prime mismatch branching
     * starting directly from the anchor interval
     * already produced by SearchDifficultyEstimator.
     */
    [[nodiscard]]
    std::vector<PrimerSearchHit>
    execute_branching(
        std::string_view primer,
        std::size_t anchor_length,
        std::size_t max_mismatches,
        const BidirectionalInterval& anchor_state
    ) const;

    /*
     * Executes anchor-locate + PackedReference
     * verification using the same precomputed
     * anchor interval.
     */
    [[nodiscard]]
    std::vector<PrimerSearchHit>
    execute_candidate(
        std::string_view primer,
        std::size_t anchor_length,
        std::size_t max_mismatches,
        const BidirectionalInterval& anchor_state
    ) const;
};

}  // namespace primerpair
