#include "primerpair/search_strategy.hpp"

#include <stdexcept>

namespace primerpair {

SearchStrategyRouter::
SearchStrategyRouter(
    const BidirectionalFMIndex& index,
    SearchDifficultyThresholds thresholds
)
    : estimator_(
          index,
          thresholds
      ) {
}

SearchStrategy
SearchStrategyRouter::choose(
    const SearchDifficulty difficulty,
    const std::size_t max_mismatches
) {
    if (max_mismatches > 3) {
        throw std::invalid_argument(
            "SearchStrategyRouter currently supports "
            "at most 3 mismatches."
        );
    }

    /*
     * --------------------------------------------------
     * Conservative routing policy v3
     *
     * Derived from the 9-repeat GRCh38.p14 chr22
     * end-to-end stability benchmark.
     *
     * Candidate routing is enabled only when:
     *
     * median speedup >= 1.10
     * AND
     * routed wins >= 7/9
     * AND
     * output equivalence is preserved.
     *
     * Stable candidate wins:
     *
     * EASY     k=3
     * MODERATE k=3
     *
     * All k=0..2 remain DirectBranching.
     * HARD and REPEAT_RICH remain DirectBranching.
     * --------------------------------------------------
     */

    if (max_mismatches < 3) {
        return
            SearchStrategy::DirectBranching;
    }

    switch (difficulty) {

        case SearchDifficulty::Easy:
        case SearchDifficulty::Moderate:

            return
                SearchStrategy::
                    AnchorCandidateVerification;

        case SearchDifficulty::Hard:
        case SearchDifficulty::RepeatRich:

            return
                SearchStrategy::DirectBranching;
    }

    return
        SearchStrategy::DirectBranching;
}


SearchStrategy
SearchStrategyRouter::choose(
    const SearchDifficulty difficulty
) noexcept {
    /*
     * Compatibility path assumes k=3.
     */
    switch (difficulty) {

        case SearchDifficulty::Easy:
        case SearchDifficulty::Moderate:
            return
                SearchStrategy::
                    AnchorCandidateVerification;

        case SearchDifficulty::Hard:
        case SearchDifficulty::RepeatRich:
            return
                SearchStrategy::DirectBranching;
    }

    return
        SearchStrategy::DirectBranching;
}

SearchStrategyDecision
SearchStrategyRouter::decide(
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    const SearchDifficultyProfile profile =
        estimator_.estimate(
            primer,
            anchor_length
        );

    return SearchStrategyDecision{
        profile,
        choose(
            profile.difficulty,
            max_mismatches
        ),
        max_mismatches
    };
}

const char* to_string(
    const SearchStrategy strategy
) noexcept {
    switch (strategy) {

        case SearchStrategy::DirectBranching:
            return "DIRECT_BRANCHING";

        case SearchStrategy::
            AnchorCandidateVerification:
            return "ANCHOR_CANDIDATE_VERIFICATION";

        case SearchStrategy::SplitSeedCandidate:
            return "SPLIT_SEED_CANDIDATE";
    }

    return "UNKNOWN";
}

}  // namespace primerpair
