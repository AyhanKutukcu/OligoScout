#include "primerpair/search_difficulty.hpp"

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

        const std::string reference =
            "ACGTACGTACGTACGT"
            "TTTTCCCCAAAAGGGG"
            "ACGTACGTACGTACGT"
            "GATCGATCGATCGATC";

        const primerpair::BidirectionalFMIndex
            index(reference);

        const primerpair::SearchDifficultyEstimator
            estimator(index);

        expect(
            estimator.classify(0) ==
                SearchDifficulty::Easy,
            "Zero occurrence is EASY"
        );

        expect(
            estimator.classify(10) ==
                SearchDifficulty::Easy,
            "EASY upper boundary"
        );

        expect(
            estimator.classify(11) ==
                SearchDifficulty::Moderate,
            "MODERATE lower boundary"
        );

        expect(
            estimator.classify(100) ==
                SearchDifficulty::Moderate,
            "MODERATE upper boundary"
        );

        expect(
            estimator.classify(101) ==
                SearchDifficulty::Hard,
            "HARD lower boundary"
        );

        expect(
            estimator.classify(1000) ==
                SearchDifficulty::Hard,
            "HARD upper boundary"
        );

        expect(
            estimator.classify(1001) ==
                SearchDifficulty::RepeatRich,
            "REPEAT_RICH lower boundary"
        );

        expect(
            std::string(
                primerpair::to_string(
                    SearchDifficulty::RepeatRich
                )
            ) ==
                "REPEAT_RICH",
            "Difficulty name"
        );

        const std::string primer =
            "GGGGACGTACGTACGT";

        const auto profile =
            estimator.estimate(
                primer,
                12
            );

        expect(
            profile.primer_length ==
                primer.size(),
            "Profile primer length"
        );

        expect(
            profile.anchor_length == 12,
            "Profile anchor length"
        );

        expect(
            profile.anchor_occurrences > 0,
            "Profile anchor occurrence count"
        );

        expect(
            profile.anchor_state.size() ==
                profile.anchor_occurrences,
            "Profile preserves anchor interval"
        );

        expect(
            profile.difficulty ==
                estimator.classify(
                    profile.anchor_occurrences
                ),
            "Profile classification"
        );

        const auto lowercase =
            estimator.estimate(
                "ggggacgtacgtacgt",
                12
            );

        expect(
            lowercase.anchor_occurrences ==
                profile.anchor_occurrences,
            "Lowercase normalization"
        );

        bool zero_anchor_rejected = false;

        try {
            static_cast<void>(
                estimator.estimate(
                    primer,
                    0
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            zero_anchor_rejected = true;
        }

        expect(
            zero_anchor_rejected,
            "Zero anchor rejection"
        );

        bool invalid_thresholds_rejected =
            false;

        try {
            const primerpair::
                SearchDifficultyEstimator bad(
                    index,
                    {
                        100,
                        10,
                        1000
                    }
                );

            static_cast<void>(bad);

        } catch (
            const std::invalid_argument&
        ) {
            invalid_thresholds_rejected =
                true;
        }

        expect(
            invalid_thresholds_rejected,
            "Invalid thresholds rejection"
        );

        std::cout
            << "All search difficulty tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
