#include "primerpair/mismatch_features.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <stdexcept>

namespace primerpair {

namespace {

void validate_primer_length(
    const std::size_t primer_length
) {
    if (
        primer_length == 0 ||
        primer_length > 64
    ) {
        throw std::invalid_argument(
            "Primer length for mismatch-mask analysis "
            "must be in the range 1..64."
        );
    }
}

std::uint64_t valid_mask_for_length(
    const std::size_t primer_length
) {
    validate_primer_length(
        primer_length
    );

    if (primer_length == 64) {
        return
            std::numeric_limits<
                std::uint64_t
            >::max();
    }

    return
        (
            std::uint64_t{1}
            <<
            primer_length
        ) -
        1;
}

void validate_mask(
    const std::uint64_t mismatch_mask,
    const std::size_t primer_length
) {
    const std::uint64_t valid =
        valid_mask_for_length(
            primer_length
        );

    if (
        (
            mismatch_mask &
            ~valid
        ) != 0
    ) {
        throw std::invalid_argument(
            "Mismatch mask contains bits outside "
            "the primer length."
        );
    }
}

std::uint64_t low_bits(
    const std::size_t count
) {
    if (count == 0) {
        return 0;
    }

    if (count >= 64) {
        return
            std::numeric_limits<
                std::uint64_t
            >::max();
    }

    return
        (
            std::uint64_t{1}
            <<
            count
        ) -
        1;
}

}  // namespace

std::size_t count_mismatches_in_5prime_window(
    const std::uint64_t mismatch_mask,
    const std::size_t primer_length,
    const std::size_t window_length
) {
    validate_mask(
        mismatch_mask,
        primer_length
    );

    const std::size_t width =
        std::min(
            primer_length,
            window_length
        );

    const std::uint64_t window_mask =
        low_bits(
            width
        );

    return
        static_cast<std::size_t>(
            std::popcount(
                mismatch_mask &
                window_mask
            )
        );
}

std::size_t count_mismatches_in_3prime_window(
    const std::uint64_t mismatch_mask,
    const std::size_t primer_length,
    const std::size_t window_length
) {
    validate_mask(
        mismatch_mask,
        primer_length
    );

    const std::size_t width =
        std::min(
            primer_length,
            window_length
        );

    if (width == 0) {
        return 0;
    }

    const std::size_t start =
        primer_length -
        width;

    std::uint64_t window_mask = 0;

    if (width == 64) {

        window_mask =
            std::numeric_limits<
                std::uint64_t
            >::max();

    } else {

        window_mask =
            low_bits(
                width
            )
            <<
            start;
    }

    return
        static_cast<std::size_t>(
            std::popcount(
                mismatch_mask &
                window_mask
            )
        );
}

bool mismatch_mask_overlaps_3prime_anchor(
    const std::uint64_t mismatch_mask,
    const std::size_t primer_length,
    const std::size_t anchor_length
) {
    if (
        anchor_length == 0 ||
        anchor_length >
            primer_length
    ) {
        throw std::invalid_argument(
            "Invalid 3-prime anchor length."
        );
    }

    return
        count_mismatches_in_3prime_window(
            mismatch_mask,
            primer_length,
            anchor_length
        ) != 0;
}

PrimerMismatchFeatures extract_mismatch_features(
    const std::uint64_t mismatch_mask,
    const std::size_t primer_length
) {
    validate_mask(
        mismatch_mask,
        primer_length
    );

    PrimerMismatchFeatures features;

    features.primer_length =
        primer_length;

    features.mismatch_mask =
        mismatch_mask;

    features.mismatch_count =
        static_cast<std::size_t>(
            std::popcount(
                mismatch_mask
            )
        );

    features.first_5_count =
        count_mismatches_in_5prime_window(
            mismatch_mask,
            primer_length,
            5
        );

    features.last_8_count =
        count_mismatches_in_3prime_window(
            mismatch_mask,
            primer_length,
            8
        );

    features.last_5_count =
        count_mismatches_in_3prime_window(
            mismatch_mask,
            primer_length,
            5
        );

    features.last_3_count =
        count_mismatches_in_3prime_window(
            mismatch_mask,
            primer_length,
            3
        );

    features.terminal_3prime_mismatch =
        (
            mismatch_mask &
            (
                std::uint64_t{1}
                <<
                (
                    primer_length -
                    1
                )
            )
        ) != 0;

    if (mismatch_mask != 0) {

        /*
         * Scan from the 3-prime-most primer position
         * toward the 5-prime end.
         */
        for (
            std::size_t offset = 0;
            offset < primer_length;
            ++offset
        ) {
            const std::size_t position =
                primer_length -
                1 -
                offset;

            if (
                (
                    mismatch_mask &
                    (
                        std::uint64_t{1}
                        <<
                        position
                    )
                ) != 0
            ) {
                features
                    .nearest_mismatch_to_3prime =
                        offset;

                break;
            }
        }
    }

    return features;
}

PrimerPairMismatchFeatures
extract_pair_mismatch_features(
    const PrimerPairHit& hit,
    const std::size_t left_primer_length,
    const std::size_t right_primer_length
) {
    PrimerPairMismatchFeatures output;

    output.left =
        extract_mismatch_features(
            hit.left_mismatch_mask,
            left_primer_length
        );

    output.right =
        extract_mismatch_features(
            hit.right_mismatch_mask,
            right_primer_length
        );

    /*
     * Treat disagreement as an invariant violation.
     */
    if (
        output.left.mismatch_count !=
            hit.left_mismatches ||
        output.right.mismatch_count !=
            hit.right_mismatches
    ) {
        throw std::logic_error(
            "PrimerPairHit mismatch counts disagree "
            "with mismatch masks."
        );
    }

    output.total_mismatches =
        output.left.mismatch_count +
        output.right.mismatch_count;

    return output;
}

}  // namespace primerpair
