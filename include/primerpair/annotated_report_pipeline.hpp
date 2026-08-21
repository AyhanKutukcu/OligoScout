#pragma once

/*
 * ANNOTATED_REPORT_PIPELINE_V1
 *
 * OligoScout
 *
 * Integration layer only.
 *
 * Existing validated components:
 *
 * PrimerPairHit
 *      ↓
 * Combined Biological Features
 *      ↓
 * Biological Risk Scoring
 *      ↓
 * PrimerPairReportRecord
 *
 * This layer introduces NO new:
 *
 * - search algorithm
 * - mismatch formula
 * - thermodynamic formula
 * - biological scoring formula
 * - PCR probability
 * - empirical calibration
 */

#include <string_view>

#include "primerpair/biological_risk_scoring.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/pcr_chemistry_profile.hpp"
#include "primerpair/report_model.hpp"
#include "primerpair/thermodynamic_backend.hpp"

namespace primerpair {


[[nodiscard]]
PrimerPairReportRecord
build_annotated_primer_pair_report_record(
    const PackedReference& reference,
    std::string_view chromosome,
    const PrimerPairHit& hit,
    std::string_view primer1,
    std::string_view primer2,
    const Primer3ThermodynamicBackend& backend,
    const BiologicalRiskScoringConfig& risk_config,
    const PcrChemistryProfile& chemistry,
    const ThermodynamicConditions& conditions = {}
);


}  // namespace primerpair
