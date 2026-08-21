/*
 * STRUCTURED_REPORT_MODEL_V1
 * Test #75
 */

#include "primerpair/report_model.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <optional>
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
make_left_features()
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


    output.nearest_mismatch_to_3prime =
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


primerpair::PrimerCombinedBiologicalFeatures
make_right_features()
{
    using namespace primerpair;


    PrimerCombinedBiologicalFeatures output;


    output.primer =
        PrimerIdentity::Primer2;


    output.reverse_strand =
        true;


    output.genomic_start =
        1200;


    output.genomic_end_exclusive =
        1220;


    output.primer_length =
        20;


    /*
     * One mismatch at original-primer
     * 5-prime position 0.
     *
     * Strict terminal-12 region remains exact.
     */
    output.mismatch_count =
        1;


    output.mismatch_fraction =
        0.05;


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


    output.nearest_mismatch_to_3prime =
        19;


    output.exact_3prime_run_length =
        19;


    output
        .normalized_3prime_positional_burden =
            1.0 / 210.0;


    output.perfect_match_tm_celsius =
        62.0;


    output.observed_binding_tm_celsius =
        61.3;


    output.delta_tm_celsius =
        0.7;


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
make_test_scoring_config()
{
    using namespace primerpair;


    BiologicalRiskScoringConfig config;


    /*
     * Unit-test fixture only.
     *
     * Not production calibration.
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


        PrimerPairHit hit;


        hit.left_primer =
            PrimerIdentity::Primer1;


        hit.right_primer =
            PrimerIdentity::Primer2;


        hit.left_position =
            1000;


        hit.right_position =
            1200;


        hit.left_mismatches =
            0;


        hit.right_mismatches =
            1;


        hit.left_mismatch_mask =
            0;


        hit.right_mismatch_mask =
            std::uint64_t{1};


        hit.amplicon_start =
            1000;


        hit.amplicon_end_exclusive =
            1220;


        hit.amplicon_length =
            220;


        PrimerPairCombinedBiologicalFeatures
            features;


        features.left =
            make_left_features();


        features.right =
            make_right_features();


        features.amplicon_start =
            hit.amplicon_start;


        features.amplicon_end_exclusive =
            hit.amplicon_end_exclusive;


        features.amplicon_length =
            hit.amplicon_length;


        const BiologicalRiskScoringConfig
            scoring_config =
                make_test_scoring_config();


        const PrimerPairBiologicalRiskScore
            risk =
                score_primer_pair_biological_risk(
                    features,
                    scoring_config
                );


        std::string dynamic_chromosome =
            "chr22";


        std::string dynamic_chemistry_name =
            "custom-lab-profile";


        PcrChemistryProfile chemistry =
            custom_pcr_chemistry_profile();


        chemistry.name =
            dynamic_chemistry_name;


        const PrimerPairReportRecord record =
            build_primer_pair_report_record(
                dynamic_chromosome,
                hit,
                features,
                risk,
                chemistry
            );


        expect(
            record.schema_version ==
                kPrimerPairReportSchemaVersion,
            "Report schema version preserved"
        );


        expect(
            record.schema_version == 1,
            "Report schema V1 exact"
        );


        expect(
            record.chromosome ==
                "chr22",
            "Chromosome preserved"
        );


        expect(
            record.hit.amplicon_start ==
                hit.amplicon_start &&
            record.hit.amplicon_end_exclusive ==
                hit.amplicon_end_exclusive &&
            record.hit.amplicon_length ==
                hit.amplicon_length,
            "Amplicon geometry preserved losslessly"
        );


        expect(
            record.hit.left_position ==
                hit.left_position &&
            record.hit.right_position ==
                hit.right_position,
            "Primer genomic positions preserved"
        );


        expect(
            record.hit.left_mismatch_mask ==
                hit.left_mismatch_mask &&
            record.hit.right_mismatch_mask ==
                hit.right_mismatch_mask,
            "Mismatch masks preserved"
        );


        expect(
            record
                .biological_features
                .left
                .mismatch_count == 0 &&
            record
                .biological_features
                .right
                .mismatch_count == 1,
            "Biological mismatch counts preserved"
        );


        expect(
            almost_equal(
                record
                    .biological_features
                    .right
                    .normalized_3prime_positional_burden,
                1.0 / 210.0
            ),
            "3-prime positional burden preserved"
        );


        expect(
            record
                .biological_features
                .right
                .exact_3prime_run_length ==
                19,
            "Exact 3-prime run preserved"
        );


        expect(
            almost_equal(
                record
                    .biological_features
                    .right
                    .delta_tm_celsius,
                0.7
            ),
            "Binding-site delta-Tm preserved"
        );


        expect(
            almost_equal(
                record
                    .biological_risk
                    .uncalibrated_ranking_score,
                risk
                    .uncalibrated_ranking_score
            ),
            "Uncalibrated ranking score preserved"
        );


        expect(
            record.score_semantics ==
                ReportScoreSemantics::
                    UncalibratedRanking,
            "Ranking semantics explicitly uncalibrated"
        );


        expect(
            record.chemistry.kind ==
                PcrChemistryKind::Custom,
            "Chemistry identity preserved"
        );


        expect(
            record.chemistry.name ==
                "custom-lab-profile",
            "Chemistry name preserved"
        );


        expect(
            !record.chemistry.source_backed,
            "Custom chemistry source status preserved"
        );


        expect(
            !record.calibration
                .ranking_score_calibrated,
            "Ranking calibration explicitly false"
        );


        expect(
            !record.calibration
                .empirical_calibration_applied,
            "Empirical calibration explicitly false"
        );


        expect(
            !record.calibration
                .pcr_probability_available,
            "PCR probability explicitly unavailable"
        );


        expect(
            record.source_search_hit_preserved,
            "Source search hit explicitly preserved"
        );


        /*
         * Verify report owns strings independently
         * of source string/string_view lifetimes.
         */
        dynamic_chromosome =
            "changed";


        dynamic_chemistry_name =
            "changed-profile";


        expect(
            record.chromosome ==
                "chr22",
            "Report chromosome uses owning storage"
        );


        expect(
            record.chemistry.name ==
                "custom-lab-profile",
            "Report chemistry name uses owning storage"
        );


        /*
         * Built-in chemistry metadata.
         */
        const auto taq =
            standard_taq_chemistry_profile();


        const auto taq_metadata =
            make_report_chemistry_metadata(
                taq
            );


        expect(
            taq_metadata.kind ==
                PcrChemistryKind::
                    StandardTaq,
            "Built-in Taq chemistry copied"
        );


        expect(
            taq_metadata.source_backed,
            "Built-in source-backed status copied"
        );


        /*
         * Empty chromosome must be rejected.
         */
        bool empty_chromosome_rejected =
            false;


        try {

            static_cast<void>(
                build_primer_pair_report_record(
                    "",
                    hit,
                    features,
                    risk,
                    chemistry
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            empty_chromosome_rejected =
                true;
        }


        expect(
            empty_chromosome_rejected,
            "Empty chromosome rejected"
        );


        /*
         * Serialization delimiter safety.
         */
        bool bad_chromosome_rejected =
            false;


        try {

            static_cast<void>(
                build_primer_pair_report_record(
                    "chr22\tbad",
                    hit,
                    features,
                    risk,
                    chemistry
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            bad_chromosome_rejected =
                true;
        }


        expect(
            bad_chromosome_rejected,
            "Chromosome delimiter contamination rejected"
        );


        /*
         * Hit geometry invariant.
         */
        bool hit_geometry_rejected =
            false;


        try {

            auto invalid_hit =
                hit;


            invalid_hit.amplicon_length =
                219;


            static_cast<void>(
                build_primer_pair_report_record(
                    "chr22",
                    invalid_hit,
                    features,
                    risk,
                    chemistry
                )
            );

        } catch (
            const std::logic_error&
        ) {

            hit_geometry_rejected =
                true;
        }


        expect(
            hit_geometry_rejected,
            "Invalid hit geometry rejected"
        );


        /*
         * Feature geometry invariant.
         */
        bool feature_geometry_rejected =
            false;


        try {

            auto invalid_features =
                features;


            invalid_features
                .amplicon_length =
                    219;


            static_cast<void>(
                build_primer_pair_report_record(
                    "chr22",
                    hit,
                    invalid_features,
                    risk,
                    chemistry
                )
            );

        } catch (
            const std::logic_error&
        ) {

            feature_geometry_rejected =
                true;
        }


        expect(
            feature_geometry_rejected,
            "Feature/hit geometry disagreement rejected"
        );


        /*
         * Mismatch-count invariant.
         */
        bool mismatch_count_rejected =
            false;


        try {

            auto invalid_features =
                features;


            invalid_features
                .right
                .mismatch_count =
                    2;


            static_cast<void>(
                build_primer_pair_report_record(
                    "chr22",
                    hit,
                    invalid_features,
                    risk,
                    chemistry
                )
            );

        } catch (
            const std::logic_error&
        ) {

            mismatch_count_rejected =
                true;
        }


        expect(
            mismatch_count_rejected,
            "Feature/hit mismatch-count disagreement rejected"
        );


        /*
         * Strand invariant.
         */
        bool strand_rejected =
            false;


        try {

            auto invalid_features =
                features;


            invalid_features
                .right
                .reverse_strand =
                    false;


            static_cast<void>(
                build_primer_pair_report_record(
                    "chr22",
                    hit,
                    invalid_features,
                    risk,
                    chemistry
                )
            );

        } catch (
            const std::logic_error&
        ) {

            strand_rejected =
                true;
        }


        expect(
            strand_rejected,
            "Invalid pair strand metadata rejected"
        );


        /*
         * Unsupported score range.
         */
        bool risk_range_rejected =
            false;


        try {

            auto invalid_risk =
                risk;


            invalid_risk
                .uncalibrated_ranking_score =
                    1.1;


            static_cast<void>(
                build_primer_pair_report_record(
                    "chr22",
                    hit,
                    features,
                    invalid_risk,
                    chemistry
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            risk_range_rejected =
                true;
        }


        expect(
            risk_range_rejected,
            "Out-of-range ranking score rejected"
        );


        /*
         * Test #75 must reject unsupported
         * calibration claims.
         */
        bool calibrated_risk_claim_rejected =
            false;


        try {

            auto invalid_chemistry =
                chemistry;


            invalid_chemistry
                .biological_risk_weights_calibrated =
                    true;


            static_cast<void>(
                build_primer_pair_report_record(
                    "chr22",
                    hit,
                    features,
                    risk,
                    invalid_chemistry
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            calibrated_risk_claim_rejected =
                true;
        }


        expect(
            calibrated_risk_claim_rejected,
            "Unsupported calibrated-risk claim rejected"
        );


        bool pcr_probability_claim_rejected =
            false;


        try {

            auto invalid_chemistry =
                chemistry;


            invalid_chemistry
                .pcr_probability_calibrated =
                    true;


            static_cast<void>(
                build_primer_pair_report_record(
                    "chr22",
                    hit,
                    features,
                    risk,
                    invalid_chemistry
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            pcr_probability_claim_rejected =
                true;
        }


        expect(
            pcr_probability_claim_rejected,
            "Unsupported PCR-probability claim rejected"
        );


        std::cout
            << "report_schema_version\t"
            << record.schema_version
            << '\n';


        std::cout
            << "report_chromosome\t"
            << record.chromosome
            << '\n';


        std::cout
            << "report_amplicon_start\t"
            << record.hit.amplicon_start
            << '\n';


        std::cout
            << "report_amplicon_end_exclusive\t"
            << record.hit
                   .amplicon_end_exclusive
            << '\n';


        std::cout
            << "report_amplicon_length\t"
            << record.hit.amplicon_length
            << '\n';


        std::cout
            << "report_right_mismatch_count\t"
            << record
                   .biological_features
                   .right
                   .mismatch_count
            << '\n';


        std::cout
            << "report_right_delta_tm_c\t"
            << record
                   .biological_features
                   .right
                   .delta_tm_celsius
            << '\n';


        std::cout
            << "report_uncalibrated_score\t"
            << record
                   .biological_risk
                   .uncalibrated_ranking_score
            << '\n';


        std::cout
            << "report_chemistry\t"
            << record.chemistry.name
            << '\n';


        std::cout
            << "REPORT_MODEL_SEARCH_HIT_PRESERVED\tYES\n";

        std::cout
            << "REPORT_MODEL_BIOLOGICAL_FEATURES_PRESERVED\tYES\n";

        std::cout
            << "REPORT_MODEL_RISK_SCORE_PRESERVED\tYES\n";

        std::cout
            << "REPORT_MODEL_CHEMISTRY_METADATA_PRESERVED\tYES\n";

        std::cout
            << "REPORT_MODEL_OWNING_STRINGS\tYES\n";

        std::cout
            << "REPORT_MODEL_SCHEMA_VERSIONED\tYES\n";

        std::cout
            << "REPORT_MODEL_SCORE_UNCALIBRATED\tYES\n";

        std::cout
            << "REPORT_MODEL_PCR_PROBABILITY\tNO\n";

        std::cout
            << "REPORT_MODEL_SEARCH_FILTERING\tNO\n";

        std::cout
            << "STRUCTURED_REPORT_MODEL_V1_COMPLETE\tYES\n";

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
