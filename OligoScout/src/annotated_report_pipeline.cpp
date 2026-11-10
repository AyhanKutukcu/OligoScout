/*
 * ANNOTATED_REPORT_PIPELINE_V1
 */

#include "primerpair/annotated_report_pipeline.hpp"

#include <stdexcept>

#include "primerpair/combined_biological_features.hpp"
#include "primerpair/primer_binding_site.hpp"

namespace primerpair {


PrimerPairReportRecord
build_annotated_primer_pair_report_record(
    const PackedReference& reference,
    const std::string_view chromosome,
    const PrimerPairHit& hit,
    const std::string_view primer1,
    const std::string_view primer2,
    const Primer3ThermodynamicBackend& backend,
    const BiologicalRiskScoringConfig& risk_config,
    const PcrChemistryProfile& chemistry,
    const ThermodynamicConditions& conditions
)
{
    if (chromosome.empty()) {

        throw std::invalid_argument(
            "Annotated report chromosome must not be empty."
        );
    }


    const auto binding_sites =
        extract_pair_binding_sites(
            reference,
            hit,
            primer1,
            primer2
        );


    const auto biological_features =
        build_combined_biological_features(
            hit,
            binding_sites,
            backend,
            conditions
        );


    const auto biological_risk =
        score_primer_pair_biological_risk(
            biological_features,
            risk_config
        );


    return build_primer_pair_report_record(
        chromosome,
        hit,
        biological_features,
        biological_risk,
        chemistry
    );
}


}  // namespace primerpair
