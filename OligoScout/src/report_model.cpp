/*
 * STRUCTURED_REPORT_MODEL_V1
 */

#include "primerpair/report_model.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace primerpair {

namespace {


void validate_chromosome(
    const std::string_view chromosome
)
{
    if (chromosome.empty()) {
        throw std::invalid_argument(
            "Report chromosome must not be empty."
        );
    }


    for (const char character : chromosome) {

        if (
            character == '\t' ||
            character == '\n' ||
            character == '\r'
        ) {
            throw std::invalid_argument(
                "Report chromosome contains "
                "a forbidden control delimiter."
            );
        }
    }
}


void require_unit_interval(
    const double value,
    const char* name
)
{
    if (
        !std::isfinite(value) ||
        value < 0.0 ||
        value > 1.0
    ) {
        throw std::invalid_argument(
            std::string(name) +
            " must be finite and in [0,1]."
        );
    }
}


void validate_hit_geometry(
    const PrimerPairHit& hit
)
{
    if (
        hit.amplicon_end_exclusive <=
        hit.amplicon_start
    ) {
        throw std::invalid_argument(
            "PrimerPairHit has invalid "
            "amplicon coordinates."
        );
    }


    const std::uint64_t expected_length =
        hit.amplicon_end_exclusive -
        hit.amplicon_start;


    if (
        hit.amplicon_length !=
        expected_length
    ) {
        throw std::logic_error(
            "PrimerPairHit amplicon length "
            "disagrees with its coordinates."
        );
    }
}


void validate_feature_consistency(
    const PrimerPairHit& hit,
    const PrimerPairCombinedBiologicalFeatures&
        features
)
{
    if (
        features.amplicon_start !=
            hit.amplicon_start ||
        features.amplicon_end_exclusive !=
            hit.amplicon_end_exclusive ||
        features.amplicon_length !=
            hit.amplicon_length
    ) {
        throw std::logic_error(
            "Combined biological feature geometry "
            "disagrees with PrimerPairHit."
        );
    }


    if (
        features.left.primer !=
            hit.left_primer ||
        features.right.primer !=
            hit.right_primer
    ) {
        throw std::logic_error(
            "Combined biological primer identity "
            "disagrees with PrimerPairHit."
        );
    }


    /*
     * Pair search invariant:
     *
     * left side of PCR product = forward
     * right side of PCR product = reverse
     */
    if (
        features.left.reverse_strand ||
        !features.right.reverse_strand
    ) {
        throw std::logic_error(
            "Combined biological strand metadata "
            "violates PCR pair orientation."
        );
    }


    if (
        features.left.genomic_start !=
            hit.left_position ||
        features.right.genomic_start !=
            hit.right_position
    ) {
        throw std::logic_error(
            "Combined biological primer position "
            "disagrees with PrimerPairHit."
        );
    }


    if (
        features.left.mismatch_count !=
            hit.left_mismatches ||
        features.right.mismatch_count !=
            hit.right_mismatches
    ) {
        throw std::logic_error(
            "Combined biological mismatch count "
            "disagrees with PrimerPairHit."
        );
    }
}


void validate_risk_score(
    const PrimerPairBiologicalRiskScore&
        score
)
{
    require_unit_interval(
        score.left
            .uncalibrated_ranking_score,
        "Left primer ranking score"
    );


    require_unit_interval(
        score.right
            .uncalibrated_ranking_score,
        "Right primer ranking score"
    );


    require_unit_interval(
        score.primer_mean_score,
        "Primer mean score"
    );


    require_unit_interval(
        score.primer_min_score,
        "Primer minimum score"
    );


    require_unit_interval(
        score.uncalibrated_ranking_score,
        "Primer-pair ranking score"
    );


    constexpr double tolerance =
        1.0e-12;


    if (
        score.primer_min_score >
            score.primer_mean_score +
            tolerance
    ) {
        throw std::logic_error(
            "Primer minimum score exceeds "
            "primer mean score."
        );
    }


    if (
        score.uncalibrated_ranking_score <
            score.primer_min_score -
                tolerance ||
        score.uncalibrated_ranking_score >
            score.primer_mean_score +
                tolerance
    ) {
        throw std::logic_error(
            "Primer-pair ranking score is "
            "outside the validated min/mean "
            "aggregation bounds."
        );
    }
}


void validate_chemistry_claims(
    const PcrChemistryProfile& chemistry
)
{
    if (
        !std::isfinite(
            chemistry.annealing_offset_celsius
        ) ||
        chemistry.annealing_offset_celsius <
            0.0
    ) {
        throw std::invalid_argument(
            "Chemistry annealing offset must "
            "be finite and non-negative."
        );
    }


    /*
     * Test #75 does not support calibrated-risk
     * or PCR-probability claims.
     *
     * Reject them rather than silently reporting
     * unsupported semantics.
     */
    if (
        chemistry
            .biological_risk_weights_calibrated
    ) {
        throw std::invalid_argument(
            "Test #75 report schema does not "
            "accept calibrated biological-risk "
            "claims."
        );
    }


    if (
        chemistry
            .pcr_probability_calibrated
    ) {
        throw std::invalid_argument(
            "Test #75 report schema does not "
            "accept calibrated PCR-probability "
            "claims."
        );
    }
}


}  // namespace


ReportChemistryMetadata
make_report_chemistry_metadata(
    const PcrChemistryProfile& chemistry
)
{
    validate_chemistry_claims(
        chemistry
    );


    ReportChemistryMetadata output;


    output.kind =
        chemistry.kind;


    /*
     * Intentional owning copy.
     */
    output.name =
        std::string(
            chemistry.name
        );


    output.annealing_rule =
        chemistry.annealing_rule;


    output.evidence_source =
        chemistry.evidence_source;


    output.annealing_offset_celsius =
        chemistry.annealing_offset_celsius;


    output.requires_chemistry_specific_tm =
        chemistry
            .requires_chemistry_specific_tm;


    output.source_backed =
        chemistry.source_backed;


    return output;
}


PrimerPairReportRecord
build_primer_pair_report_record(
    const std::string_view chromosome,
    const PrimerPairHit& hit,
    const PrimerPairCombinedBiologicalFeatures&
        biological_features,
    const PrimerPairBiologicalRiskScore&
        biological_risk,
    const PcrChemistryProfile& chemistry
)
{
    validate_chromosome(
        chromosome
    );


    validate_hit_geometry(
        hit
    );


    validate_feature_consistency(
        hit,
        biological_features
    );


    validate_risk_score(
        biological_risk
    );


    const ReportChemistryMetadata
        chemistry_metadata =
            make_report_chemistry_metadata(
                chemistry
            );


    PrimerPairReportRecord output;


    output.schema_version =
        kPrimerPairReportSchemaVersion;


    /*
     * Intentional owning copy.
     */
    output.chromosome =
        std::string(
            chromosome
        );


    output.hit =
        hit;


    output.biological_features =
        biological_features;


    output.biological_risk =
        biological_risk;


    output.score_semantics =
        ReportScoreSemantics::
            UncalibratedRanking;


    output.chemistry =
        chemistry_metadata;


    /*
     * Explicitly false in V1.
     */
    output.calibration
        .ranking_score_calibrated =
            false;


    output.calibration
        .empirical_calibration_applied =
            false;


    output.calibration
        .pcr_probability_available =
            false;


    output.source_search_hit_preserved =
        true;


    return output;
}


}  // namespace primerpair
