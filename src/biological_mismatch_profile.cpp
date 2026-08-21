// BIOLOGICAL_MISMATCH_PROFILE_V1
#include "primerpair/biological_mismatch_profile.hpp"

#include <algorithm>
#include <cstdint>

#include "primerpair/mismatch_features.hpp"

namespace primerpair {

namespace {

double positional_burden(
    const std::uint64_t mismatch_mask,
    const std::size_t primer_length
) {
    long double weighted_sum = 0.0L;

    for (
        std::size_t position = 0;
        position < primer_length;
        ++position
    ) {
        const std::uint64_t bit =
            std::uint64_t{1}
            <<
            position;

        if (
            (
                mismatch_mask &
                bit
            ) == 0
        ) {
            continue;
        }

        /*
         * 5-prime position 0 has weight 1.
         * Terminal 3-prime position has
         * weight primer_length.
         */
        weighted_sum +=
            static_cast<long double>(
                position + 1
            );
    }

    const long double length =
        static_cast<long double>(
            primer_length
        );

    const long double maximum_weight =
        (
            length *
            (
                length + 1.0L
            )
        )
        /
        2.0L;

    return
        static_cast<double>(
            weighted_sum /
            maximum_weight
        );
}

}  // namespace


BiologicalMismatchProfile
build_biological_mismatch_profile(
    const std::uint64_t mismatch_mask,
    const std::size_t primer_length
) {
    /*
     * Reuse the already validated mismatch-feature
     * implementation for:
     *
     * - primer-length validation
     * - mismatch-mask validation
     * - 3-prime coordinate semantics
     */
    const PrimerMismatchFeatures base =
        extract_mismatch_features(
            mismatch_mask,
            primer_length
        );

    BiologicalMismatchProfile output;

    output.primer_length =
        primer_length;

    output.mismatch_mask =
        mismatch_mask;

    output.mismatch_count =
        base.mismatch_count;

    output.mismatch_fraction =
        static_cast<double>(
            output.mismatch_count
        )
        /
        static_cast<double>(
            primer_length
        );

    output.last_1_count =
        count_mismatches_in_3prime_window(
            mismatch_mask,
            primer_length,
            1
        );

    output.last_2_count =
        count_mismatches_in_3prime_window(
            mismatch_mask,
            primer_length,
            2
        );

    output.last_3_count =
        base.last_3_count;

    output.last_5_count =
        base.last_5_count;

    output.last_8_count =
        base.last_8_count;

    output.last_12_count =
        count_mismatches_in_3prime_window(
            mismatch_mask,
            primer_length,
            12
        );

    output.terminal_3prime_mismatch =
        base.terminal_3prime_mismatch;

    output.nearest_mismatch_to_3prime =
        base.nearest_mismatch_to_3prime;

    if (
        output.nearest_mismatch_to_3prime
            .has_value()
    ) {
        output.exact_3prime_run_length =
            output
                .nearest_mismatch_to_3prime
                .value();
    } else {
        output.exact_3prime_run_length =
            primer_length;
    }

    output.normalized_3prime_positional_burden =
        positional_burden(
            mismatch_mask,
            primer_length
        );

    return output;
}


PrimerPairBiologicalMismatchProfile
build_pair_biological_mismatch_profile(
    const PrimerPairHit& hit,
    const std::size_t left_primer_length,
    const std::size_t right_primer_length
) {
    /*
     * This call also validates that mismatch counts
     * stored in PrimerPairHit agree with its masks.
     */
    const PrimerPairMismatchFeatures
        mismatch_features =
            extract_pair_mismatch_features(
                hit,
                left_primer_length,
                right_primer_length
            );

    PrimerPairBiologicalMismatchProfile output;

    output.left =
        build_biological_mismatch_profile(
            hit.left_mismatch_mask,
            left_primer_length
        );

    output.right =
        build_biological_mismatch_profile(
            hit.right_mismatch_mask,
            right_primer_length
        );

    output.total_mismatches =
        mismatch_features.total_mismatches;

    output
        .mean_normalized_3prime_positional_burden =
            (
                output
                    .left
                    .normalized_3prime_positional_burden
                +
                output
                    .right
                    .normalized_3prime_positional_burden
            )
            /
            2.0;

    output
        .max_normalized_3prime_positional_burden =
            std::max(
                output
                    .left
                    .normalized_3prime_positional_burden,
                output
                    .right
                    .normalized_3prime_positional_burden
            );

    return output;
}

}  // namespace primerpair
