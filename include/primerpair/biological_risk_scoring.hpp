#pragma once

/*
 * BIOLOGICAL_RISK_SCORING_V1
 *
 * OligoScout
 *
 * Transparent, deterministic and configurable
 * ranking layer over the already validated
 * Combined Biological Feature Vector.
 *
 * IMPORTANT:
 *
 * This is NOT:
 * - PCR amplification probability
 * - experimentally calibrated risk
 * - clinical probability
 *
 * Higher score means:
 *
 * "more compatible with off-target amplification
 * under the explicitly supplied heuristic scoring
 * configuration."
 *
 * No search hit is removed or modified here.
 */

#include <cstddef>

#include "primerpair/combined_biological_features.hpp"

namespace primerpair {


struct PrimerRiskNormalizedComponents {

    /*
     * 1 = exact sequence identity
     * 0 = fully mismatched
     */
    double sequence_similarity{0.0};


    /*
     * Derived from the validated positional burden.
     *
     * 1 = no mismatch burden
     * 0 = maximal normalized burden
     */
    double positional_3prime_integrity{0.0};


    /*
     * exact biological 3' run / primer length
     */
    double exact_3prime_fraction{0.0};


    /*
     * Fraction of exact bases in the indicated
     * biological terminal window.
     */
    double terminal_3_integrity{0.0};
    double terminal_5_integrity{0.0};
    double terminal_8_integrity{0.0};
    double terminal_12_integrity{0.0};


    /*
     * Thermodynamic retention derived from delta-Tm.
     *
     * retention =
     * 1 - clamp(max(delta_tm, 0) / scale, 0, 1)
     *
     * Negative delta-Tm is conservatively capped at 1.
     */
    double thermodynamic_retention{0.0};


    bool operator==(
        const PrimerRiskNormalizedComponents&
    ) const = default;
};


struct BiologicalRiskScoringWeights {

    double sequence_similarity{0.0};

    double positional_3prime_integrity{0.0};

    double exact_3prime_fraction{0.0};

    double terminal_3_integrity{0.0};
    double terminal_5_integrity{0.0};
    double terminal_8_integrity{0.0};
    double terminal_12_integrity{0.0};

    double thermodynamic_retention{0.0};


    bool operator==(
        const BiologicalRiskScoringWeights&
    ) const = default;
};


struct BiologicalRiskScoringConfig {

    BiologicalRiskScoringWeights
        weights{};


    /*
     * Positive delta-Tm corresponding to complete
     * thermodynamic-retention loss in this
     * normalization contract.
     *
     * This is a SCORING PARAMETER.
     *
     * It is not claimed to be a universal
     * biological threshold.
     */
    double delta_tm_scale_celsius{0.0};


    /*
     * Pair ranking:
     *
     * pair =
     * weakest_fraction * min(left,right)
     * +
     * (1 - weakest_fraction) * mean(left,right)
     *
     * Range 0..1.
     */
    double weakest_primer_fraction{0.0};


    bool operator==(
        const BiologicalRiskScoringConfig&
    ) const = default;
};


struct PrimerBiologicalRiskScore {

    PrimerRiskNormalizedComponents
        components{};

    /*
     * Weighted normalized score in [0,1].
     *
     * Explicitly uncalibrated.
     */
    double uncalibrated_ranking_score{0.0};


    bool operator==(
        const PrimerBiologicalRiskScore&
    ) const = default;
};


struct PrimerPairBiologicalRiskScore {

    PrimerBiologicalRiskScore left{};
    PrimerBiologicalRiskScore right{};

    double primer_mean_score{0.0};
    double primer_min_score{0.0};

    /*
     * Final pair ranking value for this layer.
     *
     * NOT a probability.
     */
    double uncalibrated_ranking_score{0.0};


    bool operator==(
        const PrimerPairBiologicalRiskScore&
    ) const = default;
};


[[nodiscard]]
PrimerRiskNormalizedComponents
normalize_biological_risk_components(
    const PrimerCombinedBiologicalFeatures& features,
    const BiologicalRiskScoringConfig& config
);


[[nodiscard]]
PrimerBiologicalRiskScore
score_primer_biological_risk(
    const PrimerCombinedBiologicalFeatures& features,
    const BiologicalRiskScoringConfig& config
);


[[nodiscard]]
PrimerPairBiologicalRiskScore
score_primer_pair_biological_risk(
    const PrimerPairCombinedBiologicalFeatures& features,
    const BiologicalRiskScoringConfig& config
);


}  // namespace primerpair
