/*
 * BIOLOGICAL_RISK_SCORING_V1
 * Test #73
 */

#include "primerpair/biological_risk_scoring.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {


void expect(
    const bool condition,
    const std::string_view message
)
{
    if (!condition) {
        throw std::runtime_error(
            std::string(message)
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}


bool almost_equal(
    const double lhs,
    const double rhs,
    const double tolerance = 1.0e-12
)
{
    return
        std::abs(
            lhs -
            rhs
        ) <=
        tolerance;
}


primerpair::PrimerCombinedBiologicalFeatures
make_exact_primer()
{
    using namespace primerpair;

    PrimerCombinedBiologicalFeatures output;

    output.primer =
        PrimerIdentity::Primer1;

    output.reverse_strand =
        false;

    output.genomic_start =
        1000;

    output.genomic_end_exclusive =
        1020;

    output.primer_length =
        20;

    output.mismatch_count =
        0;

    output.mismatch_fraction =
        0.0;

    output.last_1_count =
        0;

    output.last_2_count =
        0;

    output.last_3_count =
        0;

    output.last_5_count =
        0;

    output.last_8_count =
        0;

    output.last_12_count =
        0;

    output
        .nearest_mismatch_to_3prime =
            std::nullopt;

    output.exact_3prime_run_length =
        20;

    output
        .normalized_3prime_positional_burden =
            0.0;

    output.perfect_match_tm_celsius =
        62.0;

    output.observed_binding_tm_celsius =
        62.0;

    output.delta_tm_celsius =
        0.0;

    output.oligo_tm_celsius =
        62.0;

    output.hairpin_tm_celsius =
        10.0;

    output.homodimer_any_tm_celsius =
        5.0;

    output.homodimer_end1_tm_celsius =
        4.0;

    output.homodimer_end2_tm_celsius =
        4.0;

    return output;
}


primerpair::BiologicalRiskScoringConfig
make_reference_test_config()
{
    using namespace primerpair;

    BiologicalRiskScoringConfig config;

    /*
     * Equal weights are used ONLY as a deterministic
     * Test #73 fixture.
     *
     * They are not claimed to be production or
     * biological calibration coefficients.
     */
    config.weights.sequence_similarity =
        1.0;

    config.weights
        .positional_3prime_integrity =
            1.0;

    config.weights.exact_3prime_fraction =
        1.0;

    config.weights.terminal_3_integrity =
        1.0;

    config.weights.terminal_5_integrity =
        1.0;

    config.weights.terminal_8_integrity =
        1.0;

    config.weights.terminal_12_integrity =
        1.0;

    config.weights.thermodynamic_retention =
        1.0;

    config.delta_tm_scale_celsius =
        10.0;

    config.weakest_primer_fraction =
        0.75;

    return config;
}


}  // namespace


int main()
{
    try {

        using namespace primerpair;


        const BiologicalRiskScoringConfig config =
            make_reference_test_config();


        const PrimerCombinedBiologicalFeatures
            exact =
                make_exact_primer();


        const PrimerBiologicalRiskScore
            exact_score =
                score_primer_biological_risk(
                    exact,
                    config
                );


        expect(
            almost_equal(
                exact_score
                    .components
                    .sequence_similarity,
                1.0
            ),
            "Exact sequence similarity equals one"
        );


        expect(
            almost_equal(
                exact_score
                    .components
                    .positional_3prime_integrity,
                1.0
            ),
            "Exact 3-prime positional integrity equals one"
        );


        expect(
            almost_equal(
                exact_score
                    .components
                    .exact_3prime_fraction,
                1.0
            ),
            "Full exact 3-prime run fraction equals one"
        );


        expect(
            almost_equal(
                exact_score
                    .components
                    .terminal_3_integrity,
                1.0
            ) &&
            almost_equal(
                exact_score
                    .components
                    .terminal_5_integrity,
                1.0
            ) &&
            almost_equal(
                exact_score
                    .components
                    .terminal_8_integrity,
                1.0
            ) &&
            almost_equal(
                exact_score
                    .components
                    .terminal_12_integrity,
                1.0
            ),
            "Exact terminal-window integrity equals one"
        );


        expect(
            almost_equal(
                exact_score
                    .components
                    .thermodynamic_retention,
                1.0
            ),
            "Zero delta-Tm retention equals one"
        );


        expect(
            almost_equal(
                exact_score
                    .uncalibrated_ranking_score,
                1.0
            ),
            "Exact fixture ranking score equals one"
        );


        /*
         * A deliberately degraded synthetic fixture.
         *
         * This is used to validate scoring direction,
         * not to claim an empirical PCR probability.
         */
        PrimerCombinedBiologicalFeatures degraded =
            exact;


        degraded.mismatch_count =
            2;

        degraded.mismatch_fraction =
            0.10;

        degraded.last_3_count =
            1;

        degraded.last_5_count =
            1;

        degraded.last_8_count =
            1;

        degraded.last_12_count =
            1;

        degraded
            .nearest_mismatch_to_3prime =
                2;

        degraded.exact_3prime_run_length =
            2;

        degraded
            .normalized_3prime_positional_burden =
                0.35;

        degraded.observed_binding_tm_celsius =
            57.0;

        degraded.delta_tm_celsius =
            5.0;


        const PrimerBiologicalRiskScore
            degraded_score =
                score_primer_biological_risk(
                    degraded,
                    config
                );


        expect(
            degraded_score
                .uncalibrated_ranking_score <
            exact_score
                .uncalibrated_ranking_score,
            "Degraded fixture ranks below exact fixture"
        );


        expect(
            almost_equal(
                degraded_score
                    .components
                    .sequence_similarity,
                0.90
            ),
            "Mismatch fraction normalization exact"
        );


        expect(
            almost_equal(
                degraded_score
                    .components
                    .positional_3prime_integrity,
                0.65
            ),
            "3-prime positional normalization exact"
        );


        expect(
            almost_equal(
                degraded_score
                    .components
                    .exact_3prime_fraction,
                0.10
            ),
            "Exact 3-prime fraction normalization exact"
        );


        expect(
            almost_equal(
                degraded_score
                    .components
                    .terminal_3_integrity,
                2.0 / 3.0
            ),
            "Terminal-3 integrity exact"
        );


        expect(
            almost_equal(
                degraded_score
                    .components
                    .terminal_5_integrity,
                0.8
            ),
            "Terminal-5 integrity exact"
        );


        expect(
            almost_equal(
                degraded_score
                    .components
                    .terminal_8_integrity,
                0.875
            ),
            "Terminal-8 integrity exact"
        );


        expect(
            almost_equal(
                degraded_score
                    .components
                    .terminal_12_integrity,
                11.0 / 12.0
            ),
            "Terminal-12 integrity exact"
        );


        expect(
            almost_equal(
                degraded_score
                    .components
                    .thermodynamic_retention,
                0.5
            ),
            "delta-Tm normalization exact"
        );


        /*
         * Negative delta-Tm must not create
         * a score greater than one.
         */
        PrimerCombinedBiologicalFeatures
            stronger_than_perfect =
                exact;


        stronger_than_perfect
            .delta_tm_celsius =
                -3.0;


        const auto stronger_score =
            score_primer_biological_risk(
                stronger_than_perfect,
                config
            );


        expect(
            almost_equal(
                stronger_score
                    .components
                    .thermodynamic_retention,
                1.0
            ),
            "Negative delta-Tm capped at full retention"
        );


        /*
         * Validate individual configurable weight.
         */
        BiologicalRiskScoringConfig
            delta_only = config;


        delta_only.weights =
            BiologicalRiskScoringWeights{};


        delta_only
            .weights
            .thermodynamic_retention =
                1.0;


        const auto delta_only_score =
            score_primer_biological_risk(
                degraded,
                delta_only
            );


        expect(
            almost_equal(
                delta_only_score
                    .uncalibrated_ranking_score,
                0.5
            ),
            "Explicit thermodynamic-only weighting exact"
        );


        BiologicalRiskScoringConfig
            terminal3_only = config;


        terminal3_only.weights =
            BiologicalRiskScoringWeights{};


        terminal3_only
            .weights
            .terminal_3_integrity =
                1.0;


        const auto terminal3_score =
            score_primer_biological_risk(
                degraded,
                terminal3_only
            );


        expect(
            almost_equal(
                terminal3_score
                    .uncalibrated_ranking_score,
                2.0 / 3.0
            ),
            "Explicit terminal-3-only weighting exact"
        );


        /*
         * Pair aggregation.
         */
        PrimerPairCombinedBiologicalFeatures
            pair_features;


        pair_features.left =
            exact;

        pair_features.right =
            degraded;


        const auto pair_score =
            score_primer_pair_biological_risk(
                pair_features,
                config
            );


        expect(
            almost_equal(
                pair_score
                    .left
                    .uncalibrated_ranking_score,
                exact_score
                    .uncalibrated_ranking_score
            ),
            "Pair left score exact"
        );


        expect(
            almost_equal(
                pair_score
                    .right
                    .uncalibrated_ranking_score,
                degraded_score
                    .uncalibrated_ranking_score
            ),
            "Pair right score exact"
        );


        expect(
            almost_equal(
                pair_score.primer_min_score,
                degraded_score
                    .uncalibrated_ranking_score
            ),
            "Pair minimum primer score exact"
        );


        expect(
            pair_score.primer_min_score <=
                pair_score
                    .uncalibrated_ranking_score &&
            pair_score
                    .uncalibrated_ranking_score <=
                pair_score.primer_mean_score,
            "Pair score bounded by minimum and mean"
        );


        BiologicalRiskScoringConfig
            weakest_only = config;


        weakest_only.weakest_primer_fraction =
            1.0;


        const auto weakest_pair =
            score_primer_pair_biological_risk(
                pair_features,
                weakest_only
            );


        expect(
            almost_equal(
                weakest_pair
                    .uncalibrated_ranking_score,
                weakest_pair
                    .primer_min_score
            ),
            "Weakest-primer fraction one selects minimum"
        );


        BiologicalRiskScoringConfig
            mean_only = config;


        mean_only.weakest_primer_fraction =
            0.0;


        const auto mean_pair =
            score_primer_pair_biological_risk(
                pair_features,
                mean_only
            );


        expect(
            almost_equal(
                mean_pair
                    .uncalibrated_ranking_score,
                mean_pair.primer_mean_score
            ),
            "Weakest-primer fraction zero selects mean"
        );


        const auto repeated =
            score_primer_pair_biological_risk(
                pair_features,
                config
            );


        expect(
            repeated == pair_score,
            "Risk scoring deterministic"
        );


        /*
         * Invalid configuration checks.
         */
        bool zero_weights_rejected =
            false;


        try {

            BiologicalRiskScoringConfig invalid;

            invalid.delta_tm_scale_celsius =
                10.0;

            invalid.weakest_primer_fraction =
                0.5;


            static_cast<void>(
                score_primer_biological_risk(
                    exact,
                    invalid
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            zero_weights_rejected =
                true;
        }


        expect(
            zero_weights_rejected,
            "All-zero scoring weights rejected"
        );


        bool negative_weight_rejected =
            false;


        try {

            auto invalid =
                config;

            invalid.weights
                .sequence_similarity =
                    -1.0;


            static_cast<void>(
                score_primer_biological_risk(
                    exact,
                    invalid
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            negative_weight_rejected =
                true;
        }


        expect(
            negative_weight_rejected,
            "Negative scoring weight rejected"
        );


        bool delta_scale_rejected =
            false;


        try {

            auto invalid =
                config;

            invalid.delta_tm_scale_celsius =
                0.0;


            static_cast<void>(
                score_primer_biological_risk(
                    exact,
                    invalid
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            delta_scale_rejected =
                true;
        }


        expect(
            delta_scale_rejected,
            "Non-positive delta-Tm scale rejected"
        );


        bool pair_fraction_rejected =
            false;


        try {

            auto invalid =
                config;

            invalid.weakest_primer_fraction =
                1.5;


            static_cast<void>(
                score_primer_pair_biological_risk(
                    pair_features,
                    invalid
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            pair_fraction_rejected =
                true;
        }


        expect(
            pair_fraction_rejected,
            "Invalid pair aggregation fraction rejected"
        );


        bool feature_range_rejected =
            false;


        try {

            auto invalid_feature =
                exact;

            invalid_feature
                .mismatch_fraction =
                    1.1;


            static_cast<void>(
                score_primer_biological_risk(
                    invalid_feature,
                    config
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            feature_range_rejected =
                true;
        }


        expect(
            feature_range_rejected,
            "Invalid normalized feature range rejected"
        );


        std::cout
            << "exact_score\t"
            << exact_score
                   .uncalibrated_ranking_score
            << '\n';


        std::cout
            << "degraded_score\t"
            << degraded_score
                   .uncalibrated_ranking_score
            << '\n';


        std::cout
            << "pair_min_score\t"
            << pair_score.primer_min_score
            << '\n';


        std::cout
            << "pair_mean_score\t"
            << pair_score.primer_mean_score
            << '\n';


        std::cout
            << "pair_uncalibrated_score\t"
            << pair_score
                   .uncalibrated_ranking_score
            << '\n';


        std::cout
            << "BIOLOGICAL_RISK_NORMALIZATION_VALID\tYES\n";

        std::cout
            << "BIOLOGICAL_RISK_CONFIGURABLE_WEIGHTS\tYES\n";

        std::cout
            << "BIOLOGICAL_RISK_PAIR_AGGREGATION_VALID\tYES\n";

        std::cout
            << "BIOLOGICAL_RISK_DETERMINISTIC\tYES\n";

        std::cout
            << "BIOLOGICAL_RISK_SCORE_RANGE_0_1\tYES\n";

        std::cout
            << "BIOLOGICAL_RISK_UNCALIBRATED\tYES\n";

        std::cout
            << "BIOLOGICAL_RISK_PCR_PROBABILITY\tNO\n";

        std::cout
            << "BIOLOGICAL_RISK_EMPIRICAL_CALIBRATION\tNO\n";

        std::cout
            << "BIOLOGICAL_RISK_SEARCH_FILTERING\tNO\n";

        std::cout
            << "BIOLOGICAL_RISK_SCORING_V1_COMPLETE\tYES\n";

        std::cout
            << "ALL_CHECKS\tYES\n";


        return 0;

    } catch (
        const std::exception& e
    ) {

        std::cerr
            << "TEST_FAILURE\t"
            << e.what()
            << '\n';

        return 1;
    }
}
