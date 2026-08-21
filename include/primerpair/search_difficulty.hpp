#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "primerpair/bidirectional_fm_index.hpp"

namespace primerpair {

enum class SearchDifficulty {
    Easy,
    Moderate,
    Hard,
    RepeatRich
};

struct SearchDifficultyThresholds {
    std::uint64_t easy_max{10};
    std::uint64_t moderate_max{100};
    std::uint64_t hard_max{1000};
};

struct SearchDifficultyProfile {
    SearchDifficulty difficulty{
        SearchDifficulty::Easy
    };

    /*
     * Exact 3-prime anchor interval computed during
     * difficulty estimation.
     *
     * This is reused by the selected backend so the
     * same FM-index anchor search is not repeated.
     */
    BidirectionalInterval anchor_state{};

    std::uint64_t anchor_occurrences{0};

    std::size_t primer_length{0};
    std::size_t anchor_length{0};
};

class SearchDifficultyEstimator {
public:
    explicit SearchDifficultyEstimator(
        const BidirectionalFMIndex& index,
        SearchDifficultyThresholds thresholds = {}
    );

    [[nodiscard]]
    SearchDifficultyProfile estimate(
        std::string_view primer,
        std::size_t anchor_length = 12
    ) const;

    [[nodiscard]]
    SearchDifficulty classify(
        std::uint64_t anchor_occurrences
    ) const noexcept;

    [[nodiscard]]
    const SearchDifficultyThresholds&
    thresholds() const noexcept {
        return thresholds_;
    }

private:
    const BidirectionalFMIndex& index_;
    SearchDifficultyThresholds thresholds_;
};

[[nodiscard]]
const char* to_string(
    SearchDifficulty difficulty
) noexcept;

}  // namespace primerpair
