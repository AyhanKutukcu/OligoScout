#include "primerpair/search_difficulty.hpp"

#include <cctype>
#include <stdexcept>
#include <string>

namespace primerpair {

namespace {

std::string normalize_primer(
    const std::string_view primer
) {
    if (primer.empty()) {
        throw std::invalid_argument(
            "Primer cannot be empty."
        );
    }

    std::string normalized;
    normalized.reserve(
        primer.size()
    );

    for (const char raw : primer) {

        const char base =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        raw
                    )
                )
            );

        switch (base) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                normalized.push_back(
                    base
                );
                break;

            default:
                throw std::invalid_argument(
                    "Search difficulty primer "
                    "must contain only A/C/G/T."
                );
        }
    }

    return normalized;
}

void validate_thresholds(
    const SearchDifficultyThresholds& thresholds
) {
    if (
        thresholds.easy_max >=
        thresholds.moderate_max
    ) {
        throw std::invalid_argument(
            "easy_max must be smaller "
            "than moderate_max."
        );
    }

    if (
        thresholds.moderate_max >=
        thresholds.hard_max
    ) {
        throw std::invalid_argument(
            "moderate_max must be smaller "
            "than hard_max."
        );
    }
}

}  // namespace

SearchDifficultyEstimator::
SearchDifficultyEstimator(
    const BidirectionalFMIndex& index,
    SearchDifficultyThresholds thresholds
)
    : index_(index),
      thresholds_(thresholds) {

    validate_thresholds(
        thresholds_
    );
}

SearchDifficulty
SearchDifficultyEstimator::classify(
    const std::uint64_t anchor_occurrences
) const noexcept {

    if (
        anchor_occurrences <=
        thresholds_.easy_max
    ) {
        return SearchDifficulty::Easy;
    }

    if (
        anchor_occurrences <=
        thresholds_.moderate_max
    ) {
        return SearchDifficulty::Moderate;
    }

    if (
        anchor_occurrences <=
        thresholds_.hard_max
    ) {
        return SearchDifficulty::Hard;
    }

    return SearchDifficulty::RepeatRich;
}

SearchDifficultyProfile
SearchDifficultyEstimator::estimate(
    const std::string_view primer,
    const std::size_t anchor_length
) const {
    if (anchor_length == 0) {
        throw std::invalid_argument(
            "Anchor length must be > 0."
        );
    }

    const std::string normalized =
        normalize_primer(
            primer
        );

    if (
        anchor_length >
        normalized.size()
    ) {
        throw std::invalid_argument(
            "Anchor length cannot exceed "
            "primer length."
        );
    }

    const std::size_t anchor_begin =
        normalized.size() -
        anchor_length;

    const std::string_view anchor(
        normalized.data() +
            anchor_begin,
        anchor_length
    );

    const BidirectionalInterval anchor_state =
        index_.search(
            anchor
        );

    const std::uint64_t occurrences =
        anchor_state.size();

    return SearchDifficultyProfile{
        classify(
            occurrences
        ),
        anchor_state,
        occurrences,
        normalized.size(),
        anchor_length
    };
}

const char* to_string(
    const SearchDifficulty difficulty
) noexcept {
    switch (difficulty) {

        case SearchDifficulty::Easy:
            return "EASY";

        case SearchDifficulty::Moderate:
            return "MODERATE";

        case SearchDifficulty::Hard:
            return "HARD";

        case SearchDifficulty::RepeatRich:
            return "REPEAT_RICH";
    }

    return "UNKNOWN";
}

}  // namespace primerpair
