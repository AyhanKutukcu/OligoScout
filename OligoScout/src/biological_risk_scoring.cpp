/*
 * BIOLOGICAL_RISK_SCORING_V1
 */

#include "primerpair/biological_risk_scoring.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace primerpair {

namespace {


double clamp_unit(
    const double value
)
{
    return
        std::clamp(
            value,
            0.0,
            1.0
        );
}


void require_finite(
    const double value,
    const char* message
)
{
    if (!std::isfinite(value)) {
        throw std::invalid_argument(
            message
        );
    }
}


double weight_sum(
    const BiologicalRiskScoringWeights& weights
)
{
    const std::array<double, 8> values = {
        weights.sequence_similarity,
        weights.positional_3prime_integrity,
        weights.exact_3prime_fraction,
        weights.terminal_3_integrity,
        weights.terminal_5_integrity,
        weights.terminal_8_integrity,
        weights.terminal_12_integrity,
        weights.thermodynamic_retention
    };


    double sum = 0.0;


    for (const double value : values) {

        require_finite(
            value,
            "Biological risk weight must be finite."
        );


        if (value < 0.0) {
            throw std::invalid_argument(
                "Biological risk weight must "
                "be non-negative."
            );
        }


        sum += value;
    }


    if (
        !std::isfinite(sum) ||
        sum <= 0.0
    ) {
        throw std::invalid_argument(
            "At least one biological risk weight "
            "must be positive."
        );
    }


    return sum;
}


void validate_config(
    const BiologicalRiskScoringConfig& config
)
{
    static_cast<void>(
        weight_sum(
            config.weights
        )
    );


    require_finite(
        config.delta_tm_scale_celsius,
        "delta-Tm scale must be finite."
    );


    if (
        config.delta_tm_scale_celsius <=
        0.0
    ) {
        throw std::invalid_argument(
            "delta-Tm scale must be > 0."
        );
    }


    require_finite(
        config.weakest_primer_fraction,
        "Weakest-primer fraction must be finite."
    );


    if (
        config.weakest_primer_fraction <
            0.0 ||
        config.weakest_primer_fraction >
            1.0
    ) {
        throw std::invalid_argument(
            "Weakest-primer fraction must "
            "be in the range 0..1."
        );
    }
}


void validate_features(
    const PrimerCombinedBiologicalFeatures&
        features
)
{
    if (features.primer_length == 0) {
        throw std::invalid_argument(
            "Primer length must be > 0."
        );
    }


    require_finite(
        features.mismatch_fraction,
        "Mismatch fraction must be finite."
    );


    if (
        features.mismatch_fraction < 0.0 ||
        features.mismatch_fraction > 1.0
    ) {
        throw std::invalid_argument(
            "Mismatch fraction must be "
            "in the range 0..1."
        );
    }


    require_finite(
        features
            .normalized_3prime_positional_burden,
        "3-prime positional burden must be finite."
    );


    if (
        features
            .normalized_3prime_positional_burden <
            0.0 ||
        features
            .normalized_3prime_positional_burden >
            1.0
    ) {
        throw std::invalid_argument(
            "Normalized 3-prime positional burden "
            "must be in the range 0..1."
        );
    }


    if (
        features.exact_3prime_run_length >
        features.primer_length
    ) {
        throw std::invalid_argument(
            "Exact 3-prime run cannot exceed "
            "primer length."
        );
    }


    require_finite(
        features.delta_tm_celsius,
        "delta-Tm must be finite."
    );
}


double terminal_integrity(
    const std::size_t mismatch_count,
    const std::size_t primer_length,
    const std::size_t window_length
)
{
    const std::size_t width =
        std::min(
            primer_length,
            window_length
        );


    if (width == 0) {
        throw std::logic_error(
            "Terminal-window width is zero."
        );
    }


    if (
        mismatch_count >
        width
    ) {
        throw std::invalid_argument(
            "Terminal-window mismatch count "
            "exceeds window width."
        );
    }


    return
        1.0 -
        (
            static_cast<double>(
                mismatch_count
            ) /
            static_cast<double>(
                width
            )
        );
}


double weighted_component_score(
    const PrimerRiskNormalizedComponents&
        components,

    const BiologicalRiskScoringWeights&
        weights
)
{
    const double denominator =
        weight_sum(
            weights
        );


    const double numerator =
        components.sequence_similarity *
            weights.sequence_similarity
        +
        components.positional_3prime_integrity *
            weights.positional_3prime_integrity
        +
        components.exact_3prime_fraction *
            weights.exact_3prime_fraction
        +
        components.terminal_3_integrity *
            weights.terminal_3_integrity
        +
        components.terminal_5_integrity *
            weights.terminal_5_integrity
        +
        components.terminal_8_integrity *
            weights.terminal_8_integrity
        +
        components.terminal_12_integrity *
            weights.terminal_12_integrity
        +
        components.thermodynamic_retention *
            weights.thermodynamic_retention;


    return
        clamp_unit(
            numerator /
            denominator
        );
}


}  // namespace


PrimerRiskNormalizedComponents
normalize_biological_risk_components(
    const PrimerCombinedBiologicalFeatures& features,
    const BiologicalRiskScoringConfig& config
)
{
    validate_config(
        config
    );

    validate_features(
        features
    );


    PrimerRiskNormalizedComponents output;


    output.sequence_similarity =
        clamp_unit(
            1.0 -
            features.mismatch_fraction
        );


    output.positional_3prime_integrity =
        clamp_unit(
            1.0 -
            features
                .normalized_3prime_positional_burden
        );


    output.exact_3prime_fraction =
        clamp_unit(
            static_cast<double>(
                features.exact_3prime_run_length
            ) /
            static_cast<double>(
                features.primer_length
            )
        );


    output.terminal_3_integrity =
        terminal_integrity(
            features.last_3_count,
            features.primer_length,
            3
        );


    output.terminal_5_integrity =
        terminal_integrity(
            features.last_5_count,
            features.primer_length,
            5
        );


    output.terminal_8_integrity =
        terminal_integrity(
            features.last_8_count,
            features.primer_length,
            8
        );


    output.terminal_12_integrity =
        terminal_integrity(
            features.last_12_count,
            features.primer_length,
            12
        );


    /*
     * Positive delta-Tm means the observed genomic
     * site is weaker than the perfect-match duplex.
     *
     * Negative values are conservatively capped at
     * full retention instead of creating >1 scores.
     */
    const double nonnegative_delta_tm =
        std::max(
            features.delta_tm_celsius,
            0.0
        );


    output.thermodynamic_retention =
        clamp_unit(
            1.0 -
            (
                nonnegative_delta_tm /
                config.delta_tm_scale_celsius
            )
        );


    return output;
}


PrimerBiologicalRiskScore
score_primer_biological_risk(
    const PrimerCombinedBiologicalFeatures& features,
    const BiologicalRiskScoringConfig& config
)
{
    PrimerBiologicalRiskScore output;


    output.components =
        normalize_biological_risk_components(
            features,
            config
        );


    output.uncalibrated_ranking_score =
        weighted_component_score(
            output.components,
            config.weights
        );


    return output;
}


PrimerPairBiologicalRiskScore
score_primer_pair_biological_risk(
    const PrimerPairCombinedBiologicalFeatures& features,
    const BiologicalRiskScoringConfig& config
)
{
    validate_config(
        config
    );


    PrimerPairBiologicalRiskScore output;


    output.left =
        score_primer_biological_risk(
            features.left,
            config
        );


    output.right =
        score_primer_biological_risk(
            features.right,
            config
        );


    output.primer_mean_score =
        (
            output.left
                .uncalibrated_ranking_score +
            output.right
                .uncalibrated_ranking_score
        ) /
        2.0;


    output.primer_min_score =
        std::min(
            output.left
                .uncalibrated_ranking_score,
            output.right
                .uncalibrated_ranking_score
        );


    output.uncalibrated_ranking_score =
        clamp_unit(
            (
                config.weakest_primer_fraction *
                output.primer_min_score
            )
            +
            (
                (
                    1.0 -
                    config.weakest_primer_fraction
                )
                *
                output.primer_mean_score
            )
        );


    return output;
}


}  // namespace primerpair
