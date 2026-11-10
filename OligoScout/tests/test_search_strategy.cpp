#include "primerpair/search_strategy.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void expect(
    const bool condition,
    const std::string& name
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " + name
        );
    }

    std::cout
        << "[PASS] "
        << name
        << '\n';
}

}  // namespace

int main() {
    try {
        using primerpair::SearchDifficulty;
        using primerpair::SearchStrategy;
        using primerpair::SearchStrategyRouter;

        /*
         * EASY
         */
        expect(
            SearchStrategyRouter::choose(
                SearchDifficulty::Easy,
                0
            ) ==
                SearchStrategy::DirectBranching,
            "EASY k0 -> direct"
        );

        expect(
            SearchStrategyRouter::choose(
                SearchDifficulty::Easy,
                1
            ) ==
                SearchStrategy::DirectBranching,
            "EASY k1 -> direct"
        );

        expect(
            SearchStrategyRouter::choose(
                SearchDifficulty::Easy,
                2
            ) ==
                SearchStrategy::DirectBranching,
            "EASY k2 -> direct"
        );

        expect(
            SearchStrategyRouter::choose(
                SearchDifficulty::Easy,
                3
            ) ==
                SearchStrategy::
                    AnchorCandidateVerification,
            "EASY k3 -> candidate"
        );

        /*
         * MODERATE
         */
        expect(
            SearchStrategyRouter::choose(
                SearchDifficulty::Moderate,
                2
            ) ==
                SearchStrategy::DirectBranching,
            "MODERATE k2 -> direct"
        );

        expect(
            SearchStrategyRouter::choose(
                SearchDifficulty::Moderate,
                3
            ) ==
                SearchStrategy::
                    AnchorCandidateVerification,
            "MODERATE k3 -> candidate"
        );

        /*
         * HARD
         */
        for (std::size_t k = 0;
             k <= 3;
             ++k) {

            expect(
                SearchStrategyRouter::choose(
                    SearchDifficulty::Hard,
                    k
                ) ==
                    SearchStrategy::DirectBranching,
                "HARD k" +
                    std::to_string(k) +
                    " -> direct"
            );
        }

        /*
         * REPEAT_RICH
         */
        for (std::size_t k = 0;
             k <= 3;
             ++k) {

            expect(
                SearchStrategyRouter::choose(
                    SearchDifficulty::RepeatRich,
                    k
                ) ==
                    SearchStrategy::DirectBranching,
                "REPEAT_RICH k" +
                    std::to_string(k) +
                    " -> direct"
            );
        }

        expect(
            std::string(
                primerpair::to_string(
                    SearchStrategy::
                        AnchorCandidateVerification
                )
            ) ==
                "ANCHOR_CANDIDATE_VERIFICATION",
            "Candidate strategy name"
        );

        expect(
            std::string(
                primerpair::to_string(
                    SearchStrategy::
                        SplitSeedCandidate
                )
            ) ==
                "SPLIT_SEED_CANDIDATE",
            "Experimental split-seed name"
        );

        /*
         * Integration.
         */
        const std::string reference =
            "ACGTACGTACGTACGT"
            "TTTTCCCCAAAAGGGG"
            "ACGTACGTACGTACGT"
            "GATCGATCGATCGATC";

        const primerpair::BidirectionalFMIndex
            index(reference);

        const SearchStrategyRouter
            router(index);

        const std::string primer =
            "GGGGACGTACGTACGT";

        const auto decision =
            router.decide(
                primer,
                12,
                3
            );

        expect(
            decision.max_mismatches == 3,
            "Decision stores mismatch budget"
        );

        expect(
            decision.recommended_strategy ==
                SearchStrategyRouter::choose(
                    decision
                        .difficulty_profile
                        .difficulty,
                    3
                ),
            "Decision follows budget-aware policy"
        );

        bool excessive_budget_rejected =
            false;

        try {
            static_cast<void>(
                SearchStrategyRouter::choose(
                    SearchDifficulty::Easy,
                    4
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            excessive_budget_rejected = true;
        }

        expect(
            excessive_budget_rejected,
            "Mismatch budget >3 rejected"
        );

        std::cout
            << "All search strategy v3 tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
