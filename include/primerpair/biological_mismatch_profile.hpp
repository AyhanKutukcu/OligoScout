// BIOLOGICAL_MISMATCH_PROFILE_V1
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>

#include "primerpair/primer_pair_search.hpp"

namespace primerpair {

/*
 * Positional mismatch feature representation.
 *
 * IMPORTANT:
 * This structure is not itself a thermodynamic model
 * and is not interpreted as PCR probability.
 *
 * It converts already validated mismatch masks into
 * biologically useful, strand-normalized features.
 *
 * Primer coordinates always use the original
 * primer 5' -> 3' coordinate system.
 */
struct BiologicalMismatchProfile {
    std::size_t primer_length{0};

    std::uint64_t mismatch_mask{0};

    std::size_t mismatch_count{0};

    /*
     * mismatch_count / primer_length.
     *
     * Range: [0, 1].
     */
    double mismatch_fraction{0.0};

    /*
     * Mismatch counts in progressively larger
     * windows measured from the biological
     * 3-prime terminus.
     */
    std::size_t last_1_count{0};
    std::size_t last_2_count{0};
    std::size_t last_3_count{0};
    std::size_t last_5_count{0};
    std::size_t last_8_count{0};
    std::size_t last_12_count{0};

    bool terminal_3prime_mismatch{false};

    /*
     * 0 = mismatch at terminal 3-prime base.
     * 1 = mismatch at penultimate base.
     * ...
     *
     * nullopt = exact primer/reference match.
     */
    std::optional<std::size_t>
        nearest_mismatch_to_3prime{};

    /*
     * Number of consecutive perfectly matched
     * nucleotides extending from the biological
     * 3-prime end toward 5-prime.
     *
     * exact primer:
     *     exact_3prime_run_length == primer_length
     *
     * terminal 3-prime mismatch:
     *     exact_3prime_run_length == 0
     */
    std::size_t exact_3prime_run_length{0};

    /*
     * Normalized positional mismatch burden.
     *
     * Primer positions are weighted monotonically
     * from 5-prime to 3-prime:
     *
     *     position 0       -> weight 1
     *     ...
     *     terminal 3-prime -> weight primer_length
     *
     * Sum of mismatch-position weights is divided by:
     *
     *     1 + 2 + ... + primer_length
     *
     * Therefore:
     *
     *     exact match            -> 0
     *     every base mismatched  -> 1
     *
     * For the same number of mismatches, moving a
     * mismatch toward the 3-prime terminus strictly
     * increases this feature.
     *
     * This is a positional feature, NOT a calibrated
     * thermodynamic penalty.
     */
    double normalized_3prime_positional_burden{0.0};

    bool operator==(
        const BiologicalMismatchProfile&
    ) const = default;
};


struct PrimerPairBiologicalMismatchProfile {
    BiologicalMismatchProfile left{};
    BiologicalMismatchProfile right{};

    std::size_t total_mismatches{0};

    double mean_normalized_3prime_positional_burden{
        0.0
    };

    double max_normalized_3prime_positional_burden{
        0.0
    };

    bool operator==(
        const PrimerPairBiologicalMismatchProfile&
    ) const = default;
};


[[nodiscard]]
BiologicalMismatchProfile
build_biological_mismatch_profile(
    std::uint64_t mismatch_mask,
    std::size_t primer_length
);


[[nodiscard]]
PrimerPairBiologicalMismatchProfile
build_pair_biological_mismatch_profile(
    const PrimerPairHit& hit,
    std::size_t left_primer_length,
    std::size_t right_primer_length
);

}  // namespace primerpair
