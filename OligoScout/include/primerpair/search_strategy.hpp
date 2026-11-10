#pragma once

#include <cstddef>
#include <string_view>

#include "primerpair/search_difficulty.hpp"

namespace primerpair {

enum class SearchStrategy {
    DirectBranching,
    AnchorCandidateVerification,

    /*
     * Reserved for future repeat-rich backend.
     * Not selected by production policy v2.
     */
    SplitSeedCandidate
};

struct SearchStrategyDecision {
    SearchDifficultyProfile difficulty_profile{};

    SearchStrategy recommended_strategy{
        SearchStrategy::DirectBranching
    };

    std::size_t max_mismatches{0};
};

class SearchStrategyRouter {
public:
    explicit SearchStrategyRouter(
        const BidirectionalFMIndex& index,
        SearchDifficultyThresholds thresholds = {}
    );

    [[nodiscard]]
    SearchStrategyDecision decide(
        std::string_view primer,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 3
    ) const;

    /*
     * Budget-aware routing policy.
     */
    [[nodiscard]]
    static SearchStrategy choose(
        SearchDifficulty difficulty,
        std::size_t max_mismatches
    );

    /*
     * Compatibility overload:
     * assumes MVP maximum k=3.
     */
    [[nodiscard]]
    static SearchStrategy choose(
        SearchDifficulty difficulty
    ) noexcept;

    [[nodiscard]]
    const SearchDifficultyEstimator&
    estimator() const noexcept {
        return estimator_;
    }

private:
    SearchDifficultyEstimator estimator_;
};

[[nodiscard]]
const char* to_string(
    SearchStrategy strategy
) noexcept;

}  // namespace primerpair
