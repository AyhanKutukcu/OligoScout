#pragma once

/*
 * HTML_REPORT_V1
 *
 * OligoScout
 *
 * Self-contained human-readable HTML rendering
 * of validated PrimerPairReportRecord objects.
 *
 * This layer performs NO:
 *
 * - primer search
 * - hit filtering
 * - hit repruning
 * - biological feature calculation
 * - thermodynamic calculation
 * - ranking-score calculation
 * - PCR probability calculation
 * - empirical calibration
 *
 * Machine-readable TSV / JSON remain the
 * lossless interchange representations.
 *
 * HTML floating-point formatting is display-only.
 */

#include <cstddef>
#include <iosfwd>
#include <span>
#include <string>

#include "primerpair/report_model.hpp"

namespace primerpair {


struct HtmlReportOptions {

    std::string title{
        "OligoScout Off-Target Report"
    };

    /*
     * Human-readable display precision only.
     *
     * Does not alter values stored in the
     * underlying report model.
     */
    std::size_t display_decimal_places{
        4
    };

    bool show_machine_readable_note{
        true
    };
};


[[nodiscard]]
std::string render_html_report(
    std::span<
        const PrimerPairReportRecord
    > records,
    const HtmlReportOptions& options = {}
);


void write_html_report(
    std::ostream& output,
    std::span<
        const PrimerPairReportRecord
    > records,
    const HtmlReportOptions& options = {}
);


}  // namespace primerpair
