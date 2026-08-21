/*
 * HTML_REPORT_V1
 */

#include "primerpair/html_report.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace primerpair {

namespace {


std::string html_escape(
    const std::string_view value
)
{
    std::string output;

    output.reserve(
        value.size() +
        16
    );


    for (
        const char character :
        value
    ) {

        switch (character) {

            case '&':
                output += "&amp;";
                break;

            case '<':
                output += "&lt;";
                break;

            case '>':
                output += "&gt;";
                break;

            case '"':
                output += "&quot;";
                break;

            case '\'':
                output += "&#39;";
                break;

            default:
                output.push_back(
                    character
                );
                break;
        }
    }


    return output;
}


std::string decimal_u64(
    const std::uint64_t value
)
{
    std::array<char, 32> buffer{};


    const auto result =
        std::to_chars(
            buffer.data(),
            buffer.data() +
                buffer.size(),
            value
        );


    if (
        result.ec !=
        std::errc{}
    ) {

        throw std::runtime_error(
            "Could not format uint64 HTML value."
        );
    }


    return std::string(
        buffer.data(),
        result.ptr
    );
}


std::string decimal_size(
    const std::size_t value
)
{
    std::array<char, 32> buffer{};


    const auto result =
        std::to_chars(
            buffer.data(),
            buffer.data() +
                buffer.size(),
            value
        );


    if (
        result.ec !=
        std::errc{}
    ) {

        throw std::runtime_error(
            "Could not format size HTML value."
        );
    }


    return std::string(
        buffer.data(),
        result.ptr
    );
}


std::string decimal_u32(
    const std::uint32_t value
)
{
    std::array<char, 16> buffer{};


    const auto result =
        std::to_chars(
            buffer.data(),
            buffer.data() +
                buffer.size(),
            value
        );


    if (
        result.ec !=
        std::errc{}
    ) {

        throw std::runtime_error(
            "Could not format uint32 HTML value."
        );
    }


    return std::string(
        buffer.data(),
        result.ptr
    );
}


std::string mask_hex(
    const std::uint64_t value
)
{
    std::array<char, 32> buffer{};


    const auto result =
        std::to_chars(
            buffer.data(),
            buffer.data() +
                buffer.size(),
            value,
            16
        );


    if (
        result.ec !=
        std::errc{}
    ) {

        throw std::runtime_error(
            "Could not format mismatch mask."
        );
    }


    const std::size_t digits =
        static_cast<std::size_t>(
            result.ptr -
            buffer.data()
        );


    std::string output =
        "0x";


    if (digits < 16) {

        output.append(
            16 - digits,
            '0'
        );
    }


    output.append(
        buffer.data(),
        result.ptr
    );


    return output;
}


std::string display_double(
    const double value,
    const std::size_t decimal_places
)
{
    if (!std::isfinite(value)) {

        throw std::invalid_argument(
            "HTML report cannot display "
            "non-finite numeric values."
        );
    }


    if (decimal_places > 9) {

        throw std::invalid_argument(
            "HTML display precision must "
            "not exceed 9 decimal places."
        );
    }


    std::array<char, 128> buffer{};


    const auto result =
        std::to_chars(
            buffer.data(),
            buffer.data() +
                buffer.size(),
            value,
            std::chars_format::fixed,
            static_cast<int>(
                decimal_places
            )
        );


    if (
        result.ec !=
        std::errc{}
    ) {

        throw std::runtime_error(
            "Could not format HTML floating "
            "point value."
        );
    }


    return std::string(
        buffer.data(),
        result.ptr
    );
}


const char* primer_label(
    const PrimerIdentity identity
)
{
    switch (identity) {

        case PrimerIdentity::Primer1:
            return "Primer 1";

        case PrimerIdentity::Primer2:
            return "Primer 2";
    }


    throw std::invalid_argument(
        "Unknown primer identity."
    );
}


const char* yes_no(
    const bool value
) noexcept
{
    return value
        ? "Yes"
        : "No";
}


void validate_score(
    const double value,
    const std::string_view name
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


void validate_record(
    const PrimerPairReportRecord& record
)
{
    if (
        record.schema_version !=
        kPrimerPairReportSchemaVersion
    ) {

        throw std::invalid_argument(
            "Unsupported HTML report "
            "schema version."
        );
    }


    if (record.chromosome.empty()) {

        throw std::invalid_argument(
            "HTML report chromosome "
            "must not be empty."
        );
    }


    if (
        record.hit
            .amplicon_end_exclusive <
        record.hit
            .amplicon_start
    ) {

        throw std::invalid_argument(
            "Invalid HTML report "
            "amplicon interval."
        );
    }


    if (
        record.hit
            .amplicon_end_exclusive -
        record.hit
            .amplicon_start !=
        record.hit
            .amplicon_length
    ) {

        throw std::invalid_argument(
            "HTML report amplicon "
            "geometry disagreement."
        );
    }


    if (
        record
            .biological_features
            .amplicon_start !=
            record.hit
                .amplicon_start ||
        record
            .biological_features
            .amplicon_end_exclusive !=
            record.hit
                .amplicon_end_exclusive ||
        record
            .biological_features
            .amplicon_length !=
            record.hit
                .amplicon_length
    ) {

        throw std::invalid_argument(
            "HTML report feature/hit "
            "geometry disagreement."
        );
    }


    if (
        record
            .biological_features
            .left
            .primer !=
            record.hit
                .left_primer ||
        record
            .biological_features
            .right
            .primer !=
            record.hit
                .right_primer
    ) {

        throw std::invalid_argument(
            "HTML report primer identity "
            "disagreement."
        );
    }


    if (
        record
            .biological_features
            .left
            .reverse_strand ||
        !record
            .biological_features
            .right
            .reverse_strand
    ) {

        throw std::invalid_argument(
            "HTML report requires "
            "forward-left/reverse-right geometry."
        );
    }


    if (
        record
            .biological_features
            .left
            .genomic_start !=
            record.hit
                .left_position ||
        record
            .biological_features
            .right
            .genomic_start !=
            record.hit
                .right_position
    ) {

        throw std::invalid_argument(
            "HTML report genomic position "
            "disagreement."
        );
    }


    if (
        record
            .biological_features
            .left
            .mismatch_count !=
            record.hit
                .left_mismatches ||
        record
            .biological_features
            .right
            .mismatch_count !=
            record.hit
                .right_mismatches
    ) {

        throw std::invalid_argument(
            "HTML report mismatch-count "
            "disagreement."
        );
    }


    if (
        record.score_semantics !=
        ReportScoreSemantics::
            UncalibratedRanking
    ) {

        throw std::invalid_argument(
            "HTML V1 supports only "
            "UNCALIBRATED_RANKING semantics."
        );
    }


    if (
        record
            .calibration
            .ranking_score_calibrated ||
        record
            .calibration
            .empirical_calibration_applied ||
        record
            .calibration
            .pcr_probability_available
    ) {

        throw std::invalid_argument(
            "HTML V1 refuses unsupported "
            "calibration/probability claims."
        );
    }


    if (
        !record
            .source_search_hit_preserved
    ) {

        throw std::invalid_argument(
            "HTML V1 requires preserved "
            "source search-hit semantics."
        );
    }


    validate_score(
        record
            .biological_risk
            .left
            .uncalibrated_ranking_score,
        "Left primer ranking score"
    );


    validate_score(
        record
            .biological_risk
            .right
            .uncalibrated_ranking_score,
        "Right primer ranking score"
    );


    validate_score(
        record
            .biological_risk
            .primer_min_score,
        "Primer minimum ranking score"
    );


    validate_score(
        record
            .biological_risk
            .primer_mean_score,
        "Primer mean ranking score"
    );


    validate_score(
        record
            .biological_risk
            .uncalibrated_ranking_score,
        "Pair ranking score"
    );
}


std::string table_cell(
    const std::string_view value,
    const std::string_view css_class = {}
)
{
    std::string output =
        "<td";


    if (!css_class.empty()) {

        output +=
            " class=\"";

        output +=
            html_escape(
                css_class
            );

        output +=
            "\"";
    }


    output +=
        ">";

    output +=
        value;

    output +=
        "</td>";


    return output;
}


}  // namespace


std::string render_html_report(
    const std::span<
        const PrimerPairReportRecord
    > records,
    const HtmlReportOptions& options
)
{
    if (options.title.empty()) {

        throw std::invalid_argument(
            "HTML report title must not be empty."
        );
    }


    if (
        options.display_decimal_places >
        9
    ) {

        throw std::invalid_argument(
            "HTML report display precision "
            "must not exceed 9."
        );
    }


    for (
        const auto& record :
        records
    ) {

        validate_record(
            record
        );
    }


    std::string output;

    output.reserve(
        8192 +
        records.size() *
            2048
    );


    output +=
        "<!doctype html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" "
        "content=\"width=device-width,"
        "initial-scale=1\">\n"
        "<title>";


    output +=
        html_escape(
            options.title
        );


    output +=
        "</title>\n"
        "<style>\n"
        ":root{"
        "font-family:Inter,system-ui,"
        "-apple-system,BlinkMacSystemFont,"
        "\"Segoe UI\",sans-serif;"
        "line-height:1.45;"
        "}\n"
        "body{margin:0;padding:24px;"
        "background:#f6f7f9;"
        "color:#17191c;}\n"
        "main{max-width:1600px;"
        "margin:0 auto;}\n"
        "h1{margin:0 0 6px 0;"
        "font-size:28px;}\n"
        ".subtitle{margin:0 0 24px 0;"
        "opacity:.72;}\n"
        ".notice{padding:14px 16px;"
        "border:1px solid #bbb;"
        "border-radius:8px;"
        "background:#fff;"
        "margin:0 0 20px 0;}\n"
        ".grid{display:grid;"
        "grid-template-columns:"
        "repeat(auto-fit,minmax(190px,1fr));"
        "gap:10px;margin-bottom:20px;}\n"
        ".card{background:#fff;"
        "border:1px solid #ddd;"
        "border-radius:8px;"
        "padding:12px;}\n"
        ".label{font-size:12px;"
        "text-transform:uppercase;"
        "letter-spacing:.04em;"
        "opacity:.65;}\n"
        ".value{font-size:17px;"
        "font-weight:650;"
        "margin-top:3px;}\n"
        ".table-wrap{overflow:auto;"
        "background:#fff;"
        "border:1px solid #ddd;"
        "border-radius:8px;}\n"
        "table{border-collapse:collapse;"
        "width:100%;"
        "font-size:13px;}\n"
        "th,td{padding:9px 10px;"
        "border-bottom:1px solid #e5e5e5;"
        "white-space:nowrap;"
        "text-align:left;}\n"
        "th{position:sticky;top:0;"
        "background:#f0f1f3;"
        "font-weight:700;}\n"
        "tbody tr:hover{background:#f7f7f7;}\n"
        ".mono{font-family:"
        "ui-monospace,SFMono-Regular,"
        "Menlo,Consolas,monospace;}\n"
        ".score{font-weight:750;}\n"
        ".no-data{padding:24px;"
        "text-align:center;"
        "opacity:.7;}\n"
        "footer{margin-top:18px;"
        "font-size:12px;"
        "opacity:.7;}\n"
        "</style>\n"
        "</head>\n"
        "<body>\n"
        "<main>\n"
        "<h1>";


    output +=
        html_escape(
            options.title
        );


    output +=
        "</h1>\n"
        "<p class=\"subtitle\">"
        "OligoScout &mdash; "
        "human-readable off-target report"
        "</p>\n";


    output +=
        "<div class=\"notice\">"
        "<strong>Interpretation notice:</strong> "
        "The displayed ranking value is an "
        "<strong>UNCALIBRATED RANKING SCORE</strong> "
        "for relative ordering only. "
        "No empirical calibration has been applied "
        "and no probability estimate is available."
        "</div>\n";


    output +=
        "<section class=\"grid\">\n";


    output +=
        "<div class=\"card\">"
        "<div class=\"label\">Schema version</div>"
        "<div class=\"value\">";


    output +=
        decimal_u32(
            kPrimerPairReportSchemaVersion
        );


    output +=
        "</div></div>\n";


    output +=
        "<div class=\"card\">"
        "<div class=\"label\">Amplicon records</div>"
        "<div class=\"value\">";


    output +=
        decimal_size(
            records.size()
        );


    output +=
        "</div></div>\n";


    output +=
        "<div class=\"card\">"
        "<div class=\"label\">Score semantics</div>"
        "<div class=\"value\">"
        "UNCALIBRATED RANKING"
        "</div></div>\n";


    output +=
        "<div class=\"card\">"
        "<div class=\"label\">Empirical calibration</div>"
        "<div class=\"value\">No</div>"
        "</div>\n";


    output +=
        "<div class=\"card\">"
        "<div class=\"label\">Probability available</div>"
        "<div class=\"value\">No</div>"
        "</div>\n";


    output +=
        "</section>\n";


    output +=
        "<div class=\"table-wrap\">\n"
        "<table id=\"primerpair-results\">\n"
        "<thead>\n"
        "<tr>"
        "<th>Chromosome</th>"
        "<th>Amplicon start</th>"
        "<th>Amplicon end</th>"
        "<th>Length</th>"
        "<th>Left primer</th>"
        "<th>Right primer</th>"
        "<th>Left mismatches</th>"
        "<th>Right mismatches</th>"
        "<th>Left mask</th>"
        "<th>Right mask</th>"
        "<th>Left exact 3&#8242; run</th>"
        "<th>Right exact 3&#8242; run</th>"
        "<th>Left 3&#8242; burden</th>"
        "<th>Right 3&#8242; burden</th>"
        "<th>Left perfect Tm &deg;C</th>"
        "<th>Left observed Tm &deg;C</th>"
        "<th>Left &Delta;Tm &deg;C</th>"
        "<th>Right perfect Tm &deg;C</th>"
        "<th>Right observed Tm &deg;C</th>"
        "<th>Right &Delta;Tm &deg;C</th>"
        "<th>UNCALIBRATED RANKING SCORE</th>"
        "<th>Chemistry</th>"
        "<th>Calibrated</th>"
        "<th>Probability available</th>"
        "</tr>\n"
        "</thead>\n"
        "<tbody>\n";


    if (records.empty()) {

        output +=
            "<tr>"
            "<td colspan=\"24\" "
            "class=\"no-data\">"
            "No amplicon records."
            "</td>"
            "</tr>\n";
    }


    for (
        const auto& record :
        records
    ) {

        const auto& left =
            record
                .biological_features
                .left;


        const auto& right =
            record
                .biological_features
                .right;


        output +=
            "<tr>";


        output +=
            table_cell(
                html_escape(
                    record.chromosome
                )
            );


        output +=
            table_cell(
                decimal_u64(
                    record.hit
                        .amplicon_start
                ),
                "mono"
            );


        output +=
            table_cell(
                decimal_u64(
                    record.hit
                        .amplicon_end_exclusive
                ),
                "mono"
            );


        output +=
            table_cell(
                decimal_u64(
                    record.hit
                        .amplicon_length
                ),
                "mono"
            );


        output +=
            table_cell(
                primer_label(
                    record.hit
                        .left_primer
                )
            );


        output +=
            table_cell(
                primer_label(
                    record.hit
                        .right_primer
                )
            );


        output +=
            table_cell(
                decimal_size(
                    record.hit
                        .left_mismatches
                )
            );


        output +=
            table_cell(
                decimal_size(
                    record.hit
                        .right_mismatches
                )
            );


        output +=
            table_cell(
                mask_hex(
                    record.hit
                        .left_mismatch_mask
                ),
                "mono"
            );


        output +=
            table_cell(
                mask_hex(
                    record.hit
                        .right_mismatch_mask
                ),
                "mono"
            );


        output +=
            table_cell(
                decimal_size(
                    left
                        .exact_3prime_run_length
                )
            );


        output +=
            table_cell(
                decimal_size(
                    right
                        .exact_3prime_run_length
                )
            );


        output +=
            table_cell(
                display_double(
                    left
                        .normalized_3prime_positional_burden,
                    options
                        .display_decimal_places
                )
            );


        output +=
            table_cell(
                display_double(
                    right
                        .normalized_3prime_positional_burden,
                    options
                        .display_decimal_places
                )
            );


        output +=
            table_cell(
                display_double(
                    left
                        .perfect_match_tm_celsius,
                    options
                        .display_decimal_places
                )
            );


        output +=
            table_cell(
                display_double(
                    left
                        .observed_binding_tm_celsius,
                    options
                        .display_decimal_places
                )
            );


        output +=
            table_cell(
                display_double(
                    left
                        .delta_tm_celsius,
                    options
                        .display_decimal_places
                )
            );


        output +=
            table_cell(
                display_double(
                    right
                        .perfect_match_tm_celsius,
                    options
                        .display_decimal_places
                )
            );


        output +=
            table_cell(
                display_double(
                    right
                        .observed_binding_tm_celsius,
                    options
                        .display_decimal_places
                )
            );


        output +=
            table_cell(
                display_double(
                    right
                        .delta_tm_celsius,
                    options
                        .display_decimal_places
                )
            );


        output +=
            table_cell(
                display_double(
                    record
                        .biological_risk
                        .uncalibrated_ranking_score,
                    options
                        .display_decimal_places
                ),
                "score"
            );


        output +=
            table_cell(
                html_escape(
                    record
                        .chemistry
                        .name
                )
            );


        output +=
            table_cell(
                yes_no(
                    record
                        .calibration
                        .ranking_score_calibrated
                )
            );


        output +=
            table_cell(
                yes_no(
                    record
                        .calibration
                        .pcr_probability_available
                )
            );


        output +=
            "</tr>\n";
    }


    output +=
        "</tbody>\n"
        "</table>\n"
        "</div>\n";


    if (
        options
            .show_machine_readable_note
    ) {

        output +=
            "<footer>"
            "Machine-readable TSV and JSON outputs "
            "remain the canonical lossless interchange "
            "formats. HTML numeric rounding is display-only."
            "</footer>\n";
    }


    output +=
        "</main>\n"
        "</body>\n"
        "</html>\n";


    return output;
}


void write_html_report(
    std::ostream& output,
    const std::span<
        const PrimerPairReportRecord
    > records,
    const HtmlReportOptions& options
)
{
    output
        << render_html_report(
               records,
               options
           );


    if (!output) {

        throw std::runtime_error(
            "Failed to write HTML report."
        );
    }
}


}  // namespace primerpair
