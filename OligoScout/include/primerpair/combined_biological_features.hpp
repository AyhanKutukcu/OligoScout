#pragma once

/*
 * COMBINED_BIOLOGICAL_FEATURE_VECTOR_V1
 *
 * OligoScout
 *
 * This layer combines already validated biological
 * and thermodynamic measurements.
 *
 * It does NOT:
 *
 * - alter search hit generation
 * - reject candidate hits
 * - calculate PCR probability
 * - calculate a calibrated risk score
 *
 * The output is an interpretable feature vector
 * intended for downstream calibration / ranking.
 */

#include <cstddef>
#include <cstdint>
#include <optional>

#include "primerpair/biological_mismatch_profile.hpp"
#include "primerpair/binding_site_thermodynamics.hpp"
#include "primerpair/primer_binding_site.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/thermodynamic_backend.hpp"

namespace primerpair {


struct PrimerCombinedBiologicalFeatures {

    PrimerIdentity primer{
        PrimerIdentity::Primer1
    };

    bool reverse_strand{false};

    std::uint64_t genomic_start{0};

    std::uint64_t
        genomic_end_exclusive{0};

    std::size_t primer_length{0};


    /*
     * Sequence mismatch features.
     */
    std::size_t mismatch_count{0};

    double mismatch_fraction{0.0};

    std::size_t last_1_count{0};
    std::size_t last_2_count{0};
    std::size_t last_3_count{0};
    std::size_t last_5_count{0};
    std::size_t last_8_count{0};
    std::size_t last_12_count{0};

    std::optional<std::size_t>
        nearest_mismatch_to_3prime{};

    std::size_t
        exact_3prime_run_length{0};

    double
        normalized_3prime_positional_burden{
            0.0
        };


    /*
     * Binding-site thermodynamic features.
     *
     * delta_tm =
     * perfect_match_tm -
     * observed_binding_tm
     */
    double
        perfect_match_tm_celsius{0.0};

    double
        observed_binding_tm_celsius{0.0};

    double
        delta_tm_celsius{0.0};


    /*
     * Primer intrinsic thermodynamic features.
     */
    double
        oligo_tm_celsius{0.0};

    double
        hairpin_tm_celsius{0.0};

    double
        homodimer_any_tm_celsius{0.0};

    double
        homodimer_end1_tm_celsius{0.0};

    double
        homodimer_end2_tm_celsius{0.0};


    bool operator==(
        const PrimerCombinedBiologicalFeatures&
    ) const = default;
};


struct PrimerPairCombinedBiologicalFeatures {

    PrimerCombinedBiologicalFeatures
        left{};

    PrimerCombinedBiologicalFeatures
        right{};


    /*
     * PCR product geometry.
     */
    std::uint64_t amplicon_start{0};

    std::uint64_t
        amplicon_end_exclusive{0};

    std::uint64_t
        amplicon_length{0};


    /*
     * Pair mismatch summaries.
     */
    std::size_t
        total_mismatches{0};

    double
        mean_mismatch_fraction{0.0};

    double
        max_mismatch_fraction{0.0};

    double
        mean_normalized_3prime_positional_burden{
            0.0
        };

    double
        max_normalized_3prime_positional_burden{
            0.0
        };


    /*
     * Pair binding-site thermodynamic summaries.
     */
    double
        mean_delta_tm_celsius{0.0};

    double
        max_delta_tm_celsius{0.0};


    /*
     * Primer-primer interaction features.
     *
     * Sequence order is genomic LEFT primer
     * followed by genomic RIGHT primer.
     */
    double
        heterodimer_any_tm_celsius{0.0};

    double
        heterodimer_end1_tm_celsius{0.0};

    double
        heterodimer_end2_tm_celsius{0.0};


    bool operator==(
        const PrimerPairCombinedBiologicalFeatures&
    ) const = default;
};


/*
 * Build an interpretable biological feature vector
 * for one already-discovered PCR product.
 *
 * Search correctness is intentionally outside this
 * layer.
 */
[[nodiscard]]
PrimerPairCombinedBiologicalFeatures
build_combined_biological_features(
    const PrimerPairHit& hit,
    const PrimerPairBindingSites& sites,
    const Primer3ThermodynamicBackend& backend,
    const ThermodynamicConditions& conditions = {}
);


}  // namespace primerpair
