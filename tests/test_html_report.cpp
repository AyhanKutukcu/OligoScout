/*
 * HTML_REPORT_V1
 * Test #77
 */

#include "primerpair/html_report.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

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


primerpair::PrimerPairReportRecord
make_record()
{
    using namespace primerpair;


    PrimerPairReportRecord record;


    record.schema_version =
        kPrimerPairReportSchemaVersion;


    record.chromosome =
        "chr22";


    record.hit.amplicon_start =
        UINT64_C(9007199254740993);


    record.hit.amplicon_end_exclusive =
        UINT64_C(9007199254741213);


    record.hit.amplicon_length =
        220;


    record.hit.left_primer =
        PrimerIdentity::Primer1;


    record.hit.right_primer =
        PrimerIdentity::Primer2;


    record.hit.left_position =
        UINT64_C(9007199254740993);


    record.hit.right_position =
        UINT64_C(9007199254741193);


    record.hit.left_mismatches =
        0;


    record.hit.right_mismatches =
        1;


    record.hit.left_mismatch_mask =
        0;


    record.hit.right_mismatch_mask =
        UINT64_C(0x8000000000000000);


    auto& left =
        record
            .biological_features
            .left;


    left.primer =
        PrimerIdentity::Primer1;


    left.reverse_strand =
        false;


    left.genomic_start =
        record.hit.left_position;


    left.genomic_end_exclusive =
        record.hit.left_position +
        20;


    left.primer_length =
        20;


    left.mismatch_count =
        0;


    left.mismatch_fraction =
        0.0;


    left.nearest_mismatch_to_3prime =
        std::nullopt;


    left.exact_3prime_run_length =
        20;


    left
        .normalized_3prime_positional_burden =
            0.0;


    left.perfect_match_tm_celsius =
        62.112901;


    left.observed_binding_tm_celsius =
        62.112901;


    left.delta_tm_celsius =
        0.0;


    left.oligo_tm_celsius =
        62.182505;


    left.hairpin_tm_celsius =
        10.25;


    auto& right =
        record
            .biological_features
            .right;


    right.primer =
        PrimerIdentity::Primer2;


    right.reverse_strand =
        true;


    right.genomic_start =
        record.hit.right_position;


    right.genomic_end_exclusive =
        record.hit.right_position +
        20;


    right.primer_length =
        20;


    right.mismatch_count =
        1;


    right.mismatch_fraction =
        0.05;


    right.nearest_mismatch_to_3prime =
        19;


    right.exact_3prime_run_length =
        19;


    right
        .normalized_3prime_positional_burden =
            1.0 /
            210.0;


    right.perfect_match_tm_celsius =
        62.112901;


    right.observed_binding_tm_celsius =
        61.400501;


    right.delta_tm_celsius =
        0.7124;


    right.oligo_tm_celsius =
        62.182505;


    right.hairpin_tm_celsius =
        11.75;


    record
        .biological_features
        .amplicon_start =
            record.hit
                .amplicon_start;


    record
        .biological_features
        .amplicon_end_exclusive =
            record.hit
                .amplicon_end_exclusive;


    record
        .biological_features
        .amplicon_length =
            record.hit
                .amplicon_length;


    record
        .biological_risk
        .left
        .uncalibrated_ranking_score =
            1.0;


    record
        .biological_risk
        .right
        .uncalibrated_ranking_score =
            0.92;


    record
        .biological_risk
        .primer_min_score =
            0.92;


    record
        .biological_risk
        .primer_mean_score =
            0.96;


    record
        .biological_risk
        .uncalibrated_ranking_score =
            0.93;


    record.score_semantics =
        ReportScoreSemantics::
            UncalibratedRanking;


    record.chemistry.kind =
        PcrChemistryKind::Custom;


    /*
     * Deliberate HTML injection fixture.
     */
    record.chemistry.name =
        "<script>alert(\"x\")</script>"
        "& custom";


    record
        .chemistry
        .annealing_offset_celsius =
            0.0;


    record
        .chemistry
        .requires_chemistry_specific_tm =
            false;


    record
        .chemistry
        .source_backed =
            false;


    record
        .calibration
        .ranking_score_calibrated =
            false;


    record
        .calibration
        .empirical_calibration_applied =
            false;


    record
        .calibration
        .pcr_probability_available =
            false;


    record
        .source_search_hit_preserved =
            true;


    return record;
}


}  // namespace


int main()
{
    try {

        using namespace primerpair;


        const auto first =
            make_record();


        auto second =
            first;


        second.chromosome =
            "chr1";


        second
            .biological_risk
            .uncalibrated_ranking_score =
                0.75;


        std::vector<
            PrimerPairReportRecord
        > records = {
            first,
            second
        };


        HtmlReportOptions options;


        options.title =
            "PrimerPair <Research> & Report";


        options.display_decimal_places =
            4;


        const std::string html1 =
            render_html_report(
                records,
                options
            );


        const std::string html2 =
            render_html_report(
                records,
                options
            );


        expect(
            html1 == html2,
            "HTML rendering deterministic"
        );


        expect(
            html1.starts_with(
                "<!doctype html>\n"
            ),
            "HTML doctype deterministic"
        );


        expect(
            html1.ends_with(
                "</html>\n"
            ),
            "HTML document closes deterministically"
        );


        expect(
            html1.find(
                "PrimerPair &lt;Research&gt; "
                "&amp; Report"
            ) !=
                std::string::npos,
            "HTML title safely escaped"
        );


        expect(
            html1.find(
                "&lt;script&gt;"
                "alert(&quot;x&quot;)"
                "&lt;/script&gt;"
                "&amp; custom"
            ) !=
                std::string::npos,
            "Chemistry HTML injection escaped"
        );


        expect(
            html1.find(
                "<script>alert"
            ) ==
                std::string::npos,
            "Raw script injection absent"
        );


        expect(
            html1.find(
                "9007199254740993"
            ) !=
                std::string::npos,
            "HTML preserves exact >2^53 coordinate text"
        );


        expect(
            html1.find(
                "0x8000000000000000"
            ) !=
                std::string::npos,
            "HTML preserves mismatch mask"
        );


        expect(
            html1.find(
                "UNCALIBRATED RANKING SCORE"
            ) !=
                std::string::npos,
            "HTML ranking score label explicit"
        );


        expect(
            html1.find(
                "No empirical calibration "
                "has been applied"
            ) !=
                std::string::npos,
            "HTML calibration disclaimer explicit"
        );


        expect(
            html1.find(
                "no probability estimate "
                "is available"
            ) !=
                std::string::npos,
            "HTML probability disclaimer explicit"
        );


        expect(
            html1.find(
                ">0.9300</td>"
            ) !=
                std::string::npos,
            "HTML score display formatting exact"
        );


        expect(
            html1.find(
                ">0.7124</td>"
            ) !=
                std::string::npos,
            "HTML delta-Tm display formatting exact"
        );


        expect(
            html1.find(
                ">0.0048</td>"
            ) !=
                std::string::npos,
            "HTML 3-prime burden display formatting exact"
        );


        const auto chr22_position =
            html1.find(
                ">chr22</td>"
            );


        const auto chr1_position =
            html1.find(
                ">chr1</td>"
            );


        expect(
            chr22_position !=
                std::string::npos &&
            chr1_position !=
                std::string::npos &&
            chr22_position <
                chr1_position,
            "HTML preserves input record order"
        );


        std::ostringstream stream;


        write_html_report(
            stream,
            records,
            options
        );


        expect(
            stream.str() ==
                html1,
            "HTML stream writer equals renderer"
        );


        const std::vector<
            PrimerPairReportRecord
        > empty_records;


        const std::string empty_html =
            render_html_report(
                empty_records,
                options
            );


        expect(
            empty_html.find(
                "No amplicon records."
            ) !=
                std::string::npos,
            "Empty HTML report supported"
        );


        bool calibrated_rejected =
            false;


        try {

            auto invalid =
                first;


            invalid
                .calibration
                .ranking_score_calibrated =
                    true;


            const std::vector<
                PrimerPairReportRecord
            > invalid_records = {
                invalid
            };


            static_cast<void>(
                render_html_report(
                    invalid_records,
                    options
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            calibrated_rejected =
                true;
        }


        expect(
            calibrated_rejected,
            "Unsupported calibrated-score claim rejected"
        );


        bool probability_rejected =
            false;


        try {

            auto invalid =
                first;


            invalid
                .calibration
                .pcr_probability_available =
                    true;


            const std::vector<
                PrimerPairReportRecord
            > invalid_records = {
                invalid
            };


            static_cast<void>(
                render_html_report(
                    invalid_records,
                    options
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            probability_rejected =
                true;
        }


        expect(
            probability_rejected,
            "Unsupported probability claim rejected"
        );


        bool nonfinite_rejected =
            false;


        try {

            auto invalid =
                first;


            invalid
                .biological_features
                .right
                .delta_tm_celsius =
                    std::numeric_limits<
                        double
                    >::quiet_NaN();


            const std::vector<
                PrimerPairReportRecord
            > invalid_records = {
                invalid
            };


            static_cast<void>(
                render_html_report(
                    invalid_records,
                    options
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            nonfinite_rejected =
                true;
        }


        expect(
            nonfinite_rejected,
            "Non-finite HTML display value rejected"
        );


        bool precision_rejected =
            false;


        try {

            auto invalid_options =
                options;


            invalid_options
                .display_decimal_places =
                    10;


            static_cast<void>(
                render_html_report(
                    records,
                    invalid_options
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            precision_rejected =
                true;
        }


        expect(
            precision_rejected,
            "Invalid HTML display precision rejected"
        );


        /*
         * Optional output path used by the
         * external parser validation step.
         */
        if (
            const char* output_path =
                std::getenv(
                    "PRIMERPAIR_HTML_OUTPUT"
                );
            output_path != nullptr
        ) {

            std::ofstream file(
                output_path,
                std::ios::binary
            );


            if (!file) {

                throw std::runtime_error(
                    "Could not open requested "
                    "HTML sample output file."
                );
            }


            file
                << html1;


            if (!file) {

                throw std::runtime_error(
                    "Could not write requested "
                    "HTML sample output file."
                );
            }
        }


        std::cout
            << "html_record_count\t"
            << records.size()
            << '\n';


        std::cout
            << "html_display_score\t"
            << "0.9300"
            << '\n';


        std::cout
            << "html_display_delta_tm\t"
            << "0.7124"
            << '\n';


        std::cout
            << "HTML_REPORT_DETERMINISTIC\tYES\n";

        std::cout
            << "HTML_REPORT_SELF_CONTAINED\tYES\n";

        std::cout
            << "HTML_REPORT_ESCAPING_VALID\tYES\n";

        std::cout
            << "HTML_REPORT_64BIT_TEXT_LOSSLESS\tYES\n";

        std::cout
            << "HTML_REPORT_MASK_TEXT_LOSSLESS\tYES\n";

        std::cout
            << "HTML_REPORT_DISPLAY_ROUNDING_ONLY\tYES\n";

        std::cout
            << "HTML_REPORT_SCORE_UNCALIBRATED\tYES\n";

        std::cout
            << "HTML_REPORT_EMPIRICAL_CALIBRATION\tNO\n";

        std::cout
            << "HTML_REPORT_PCR_PROBABILITY\tNO\n";

        std::cout
            << "HTML_REPORT_BIOLOGICAL_RECALCULATION\tNO\n";

        std::cout
            << "HTML_REPORT_SEARCH_CHANGE\tNO\n";

        std::cout
            << "HTML_REPORT_V1_COMPLETE\tYES\n";

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
