#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>

#include "primerpair/primer_pair_search.hpp"

namespace primerpair {

struct PrimerMismatchFeatures {
    std::size_t primer_length{0};

    std::uint64_t mismatch_mask{0};

    std::size_t mismatch_count{0};

    /*
     * Number of mismatches in biologically useful
     * primer windows.
     */
    std::size_t first_5_count{0};

    std::size_t last_8_count{0};
    std::size_t last_5_count{0};
    std::size_t last_3_count{0};

    bool terminal_3prime_mismatch{false};

    /*
     * Distance from the primer's 3-prime terminal base.
     *
     * 0 = terminal 3-prime base
     * 1 = penultimate base
     * ...
     *
     * nullopt = exact match
     */
    std::optional<std::size_t>
        nearest_mismatch_to_3prime{};

    bool operator==(
        const PrimerMismatchFeatures&
    ) const = default;
};

struct PrimerPairMismatchFeatures {
    PrimerMismatchFeatures left{};
    PrimerMismatchFeatures right{};

    std::size_t total_mismatches{0};

    bool operator==(
        const PrimerPairMismatchFeatures&
    ) const = default;
};

[[nodiscard]]
std::size_t count_mismatches_in_5prime_window(
    std::uint64_t mismatch_mask,
    std::size_t primer_length,
    std::size_t window_length
);

[[nodiscard]]
std::size_t count_mismatches_in_3prime_window(
    std::uint64_t mismatch_mask,
    std::size_t primer_length,
    std::size_t window_length
);

[[nodiscard]]
bool mismatch_mask_overlaps_3prime_anchor(
    std::uint64_t mismatch_mask,
    std::size_t primer_length,
    std::size_t anchor_length
);

[[nodiscard]]
PrimerMismatchFeatures extract_mismatch_features(
    std::uint64_t mismatch_mask,
    std::size_t primer_length
);

[[nodiscard]]
PrimerPairMismatchFeatures
extract_pair_mismatch_features(
    const PrimerPairHit& hit,
    std::size_t left_primer_length,
    std::size_t right_primer_length
);

}  // namespace primerpair
