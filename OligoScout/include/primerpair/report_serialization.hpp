#pragma once

/*
 * REPORT_SERIALIZATION_V1
 *
 * OligoScout
 *
 * Deterministic serialization of the validated
 * PrimerPairReportRecord.
 *
 * This layer performs NO:
 *
 * - search
 * - filtering
 * - biological scoring
 * - thermodynamic calculation
 * - PCR-probability inference
 *
 * 64-bit genomic coordinates are represented as
 * decimal strings in JSON so consumers based on
 * IEEE-754 doubles cannot silently lose precision.
 *
 * 64-bit mismatch masks are represented as fixed
 * 16-digit hexadecimal strings.
 */

#include <iosfwd>
#include <span>
#include <string>

#include "primerpair/report_model.hpp"

namespace primerpair {


[[nodiscard]]
std::string report_tsv_header();


[[nodiscard]]
std::string serialize_report_tsv_row(
    const PrimerPairReportRecord& record
);


[[nodiscard]]
std::string serialize_report_json_object(
    const PrimerPairReportRecord& record
);


void write_report_tsv(
    std::ostream& output,
    std::span<
        const PrimerPairReportRecord
    > records
);


void write_report_json(
    std::ostream& output,
    std::span<
        const PrimerPairReportRecord
    > records
);


}  // namespace primerpair
