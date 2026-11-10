#pragma once

/*
 * STRUCTURED_REPORT_MODEL_V1
 *
 * OligoScout
 *
 * Stable per-amplicon reporting representation.
 *
 * This layer performs NO:
 *
 * - primer search
 * - hit filtering
 * - mismatch calculation
 * - thermodynamic calculation
 * - risk recalculation
 * - PCR probability estimation
 *
 * It collects already validated results into an
 * owning record suitable for later serialization.
 */

#include <cstdint>
#include <string>
#include <string_view>

#include "primerpair/biological_risk_scoring.hpp"
#include "primerpair/combined_biological_features.hpp"
#include "primerpair/pcr_chemistry_profile.hpp"
#include "primerpair/primer_pair_search.hpp"

namespace primerpair {


inline constexpr std::uint32_t
    kPrimerPairReportSchemaVersion = 1;


enum class ReportScoreSemantics : std::uint8_t {

    /*
     * Test #73 scoring semantics.
     *
     * Explicitly NOT probability.
     */
    UncalibratedRanking = 0
};


struct ReportCalibrationStatus {

    /*
     * Must remain false until an external
     * calibration milestone validates otherwise.
     */
    bool ranking_score_calibrated{
        false
    };


    bool empirical_calibration_applied{
        false
    };


    bool pcr_probability_available{
        false
    };


    bool operator==(
        const ReportCalibrationStatus&
    ) const = default;
};


struct ReportChemistryMetadata {

    PcrChemistryKind kind{
        PcrChemistryKind::Custom
    };


    /*
     * Owning string.
     *
     * PcrChemistryProfile uses string_view.
     * A report record must remain valid after
     * the source profile object goes out of scope.
     */
    std::string name{};


    AnnealingTemperatureRule
        annealing_rule{
            AnnealingTemperatureRule::
                UserSpecified
        };


    ChemistryEvidenceSource
        evidence_source{
            ChemistryEvidenceSource::
                UserDefined
        };


    double annealing_offset_celsius{
        0.0
    };


    bool requires_chemistry_specific_tm{
        false
    };


    bool source_backed{
        false
    };


    bool operator==(
        const ReportChemistryMetadata&
    ) const = default;
};


struct PrimerPairReportRecord {

    std::uint32_t schema_version{
        kPrimerPairReportSchemaVersion
    };


    /*
     * Owning chromosome/shard name.
     *
     * Examples:
     * chr1
     * chr22
     * chrX
     */
    std::string chromosome{};


    /*
     * Exact search/pairing result.
     *
     * No search-field transformation.
     */
    PrimerPairHit hit{};


    /*
     * Test #72 validated feature representation.
     */
    PrimerPairCombinedBiologicalFeatures
        biological_features{};


    /*
     * Test #73 validated uncalibrated ranking.
     */
    PrimerPairBiologicalRiskScore
        biological_risk{};


    ReportScoreSemantics
        score_semantics{
            ReportScoreSemantics::
                UncalibratedRanking
        };


    /*
     * Test #74 chemistry metadata.
     *
     * Stored independently of search strategy.
     */
    ReportChemistryMetadata chemistry{};


    ReportCalibrationStatus calibration{};


    /*
     * Explicit provenance statement for this
     * schema generation.
     */
    bool source_search_hit_preserved{
        true
    };
};


[[nodiscard]]
ReportChemistryMetadata
make_report_chemistry_metadata(
    const PcrChemistryProfile& chemistry
);


[[nodiscard]]
PrimerPairReportRecord
build_primer_pair_report_record(
    std::string_view chromosome,
    const PrimerPairHit& hit,
    const PrimerPairCombinedBiologicalFeatures&
        biological_features,
    const PrimerPairBiologicalRiskScore&
        biological_risk,
    const PcrChemistryProfile& chemistry
);


}  // namespace primerpair
