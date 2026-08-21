#include "primerpair/mismatch_features.hpp"

#include <cstddef>
#include <cstdint>
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
        /*
         * 20-mer:
         *
         * mismatches at:
         * 0  = 5-prime terminal
         * 4
         * 17 = inside last 3
         * 19 = 3-prime terminal
         */
        const std::uint64_t mask =
            (
                std::uint64_t{1} << 0
            ) |
            (
                std::uint64_t{1} << 4
            ) |
            (
                std::uint64_t{1} << 17
            ) |
            (
                std::uint64_t{1} << 19
            );

        const auto features =
            primerpair::
                extract_mismatch_features(
                    mask,
                    20
                );

        expect(
            features.mismatch_count == 4,
            "Total mismatch count"
        );

        expect(
            features.first_5_count == 2,
            "First 5-prime window count"
        );

        expect(
            features.last_8_count == 2,
            "Last 8 count"
        );

        expect(
            features.last_5_count == 2,
            "Last 5 count"
        );

        expect(
            features.last_3_count == 2,
            "Last 3 count"
        );

        expect(
            features.terminal_3prime_mismatch,
            "Terminal 3-prime mismatch"
        );

        expect(
            features.nearest_mismatch_to_3prime
                .has_value(),
            "Nearest 3-prime mismatch exists"
        );

        expect(
            features.nearest_mismatch_to_3prime
                .value() == 0,
            "Nearest mismatch distance zero"
        );

        /*
         * Current STRICT profile:
         *
         * 20-mer + exact 12-nt 3-prime anchor means
         * the closest legal mismatch position is 7.
         *
         * Position 7 is exactly 12 nt away from the
         * terminal 3-prime base.
         */
        const std::uint64_t strict_mask =
            std::uint64_t{1} << 7;

        const auto strict =
            primerpair::
                extract_mismatch_features(
                    strict_mask,
                    20
                );

        expect(
            strict.nearest_mismatch_to_3prime
                .value() == 12,
            "Strict nearest mismatch is 12 nt from 3-prime end"
        );

        expect(
            !primerpair::
                mismatch_mask_overlaps_3prime_anchor(
                    strict_mask,
                    20,
                    12
                ),
            "Strict mismatch outside exact 12-nt anchor"
        );

        expect(
            primerpair::
                mismatch_mask_overlaps_3prime_anchor(
                    std::uint64_t{1} << 8,
                    20,
                    12
                ),
            "Position 8 belongs to exact 12-nt anchor"
        );

        const auto exact =
            primerpair::
                extract_mismatch_features(
                    0,
                    20
                );

        expect(
            exact.mismatch_count == 0,
            "Exact hit has zero mismatches"
        );

        expect(
            !exact.nearest_mismatch_to_3prime
                .has_value(),
            "Exact hit has no nearest mismatch"
        );

        /*
         * Pair-level extraction.
         */
        primerpair::PrimerPairHit pair_hit;

        pair_hit.left_mismatches = 1;
        pair_hit.right_mismatches = 1;

        pair_hit.left_mismatch_mask =
            std::uint64_t{1} << 0;

        pair_hit.right_mismatch_mask =
            std::uint64_t{1} << 7;

        const auto pair_features =
            primerpair::
                extract_pair_mismatch_features(
                    pair_hit,
                    20,
                    20
                );

        expect(
            pair_features.total_mismatches == 2,
            "Pair total mismatch count"
        );

        expect(
            pair_features.left.first_5_count == 1,
            "Pair left 5-prime feature"
        );

        expect(
            pair_features.right
                .nearest_mismatch_to_3prime
                .value() == 12,
            "Pair right nearest 3-prime distance"
        );

        /*
         * Mask outside primer length must fail.
         */
        bool invalid_mask_rejected = false;

        try {
            static_cast<void>(
                primerpair::
                    extract_mismatch_features(
                        std::uint64_t{1} << 20,
                        20
                    )
            );

        } catch (
            const std::invalid_argument&
        ) {
            invalid_mask_rejected = true;
        }

        expect(
            invalid_mask_rejected,
            "Out-of-range mask bit rejected"
        );

        /*
         * Count/mask disagreement must fail.
         */
        bool inconsistent_pair_rejected =
            false;

        try {
            primerpair::PrimerPairHit bad;

            bad.left_mismatches = 2;
            bad.left_mismatch_mask = 1;

            static_cast<void>(
                primerpair::
                    extract_pair_mismatch_features(
                        bad,
                        20,
                        20
                    )
            );

        } catch (
            const std::logic_error&
        ) {
            inconsistent_pair_rejected =
                true;
        }

        expect(
            inconsistent_pair_rejected,
            "Pair count-mask disagreement rejected"
        );

        std::cout
            << "All mismatch feature tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
