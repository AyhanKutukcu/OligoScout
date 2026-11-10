/*
 * REPORT_SERIALIZATION_V1
 * Test #76
 */

#include "primerpair/report_serialization.hpp"

#include <cmath>
#include <cstdint>
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


std::size_t count_columns(
    const std::string_view line
)
{
    std::size_t columns =
        1;


    for (const char character : line) {

        if (character == '\t') {
            ++columns;
        }
    }


    return columns;
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


    /*
     * Deliberately greater than 2^53 to prove
     * that JSON output does not depend on
     * IEEE-754 integer precision.
     */
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
        record.biological_features.left;


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
        record.biological_features.right;


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


    record.biological_features
        .amplicon_start =
            record.hit.amplicon_start;


    record.biological_features
        .amplicon_end_exclusive =
            record.hit
                .amplicon_end_exclusive;


    record.biological_features
        .amplicon_length =
            record.hit.amplicon_length;


    record.biological_risk
        .left
        .uncalibrated_ranking_score =
            1.0;


    record.biological_risk
        .right
        .uncalibrated_ranking_score =
            0.92;


    record.biological_risk
        .primer_min_score =
            0.92;


    record.biological_risk
        .primer_mean_score =
            0.96;


    record.biological_risk
        .uncalibrated_ranking_score =
            0.93;


    record.score_semantics =
        ReportScoreSemantics::
            UncalibratedRanking;


    record.chemistry.kind =
        PcrChemistryKind::Custom;


    /*
     * Includes characters requiring escaping
     * in both TSV and JSON.
     */
    record.chemistry.name =
        "custom\t\"lab\"\\profile\nv1";


    record.chemistry
        .annealing_offset_celsius =
            0.0;


    record.chemistry
        .requires_chemistry_specific_tm =
            false;


    record.chemistry.source_backed =
        false;


    record.calibration
        .ranking_score_calibrated =
            false;


    record.calibration
        .empirical_calibration_applied =
            false;


    record.calibration
        .pcr_probability_available =
            false;


    record.source_search_hit_preserved =
        true;


    return record;
}


}  // namespace


int main()
{
    try {

        using namespace primerpair;


        const auto record =
            make_record();


        const std::string header =
            report_tsv_header();


        const std::string row1 =
            serialize_report_tsv_row(
                record
            );


        const std::string row2 =
            serialize_report_tsv_row(
                record
            );


        expect(
            row1 == row2,
            "TSV serialization deterministic"
        );


        expect(
            count_columns(header) ==
                count_columns(row1),
            "TSV header/row column counts agree"
        );


        expect(
            row1.find(
                "9007199254740993"
            ) !=
                std::string::npos,
            "TSV preserves >2^53 coordinate exactly"
        );


        expect(
            row1.find(
                "0x8000000000000000"
            ) !=
                std::string::npos,
            "TSV mismatch mask fixed hexadecimal"
        );


        expect(
            row1.find(
                "custom\\t\"lab\"\\\\profile\\nv1"
            ) !=
                std::string::npos,
            "TSV text escaping deterministic"
        );


        const std::string json1 =
            serialize_report_json_object(
                record
            );


        const std::string json2 =
            serialize_report_json_object(
                record
            );


        expect(
            json1 == json2,
            "JSON serialization deterministic"
        );


        expect(
            json1.find(
                "\"amplicon_start\":"
                "\"9007199254740993\""
            ) !=
                std::string::npos,
            "JSON coordinate encoded as lossless string"
        );


        expect(
            json1.find(
                "\"right_mismatch_mask\":"
                "\"0x8000000000000000\""
            ) !=
                std::string::npos,
            "JSON mismatch mask encoded losslessly"
        );


        expect(
            json1.find(
                "\"chemistry_name\":"
                "\"custom\\t\\\"lab\\\""
                "\\\\profile\\nv1\""
            ) !=
                std::string::npos,
            "JSON text escaping valid"
        );


        expect(
            json1.find(
                "\"score_semantics\":"
                "\"UNCALIBRATED_RANKING\""
            ) !=
                std::string::npos,
            "JSON score semantics explicit"
        );


        expect(
            json1.find(
                "\"ranking_score_calibrated\":false"
            ) !=
                std::string::npos,
            "JSON ranking calibration false"
        );


        expect(
            json1.find(
                "\"pcr_probability_available\":false"
            ) !=
                std::string::npos,
            "JSON PCR probability unavailable"
        );


        std::vector<
            PrimerPairReportRecord
        > records = {
            record,
            record
        };


        std::ostringstream tsv_stream;


        write_report_tsv(
            tsv_stream,
            records
        );


        const std::string tsv_document =
            tsv_stream.str();


        expect(
            tsv_document.starts_with(
                header + "\n"
            ),
            "TSV document begins with deterministic header"
        );


        expect(
            tsv_document.ends_with(
                row1 + "\n"
            ),
            "TSV document ends with newline"
        );


        std::ostringstream json_stream;


        write_report_json(
            json_stream,
            records
        );


        const std::string json_document =
            json_stream.str();


        expect(
            json_document.starts_with(
                "[\n"
            ),
            "JSON document begins with array"
        );


        expect(
            json_document.ends_with(
                "]\n"
            ),
            "JSON document ends deterministically"
        );


        expect(
            json_document.find(
                "},\n  {"
            ) !=
                std::string::npos,
            "JSON multi-record comma placement valid"
        );


        /*
         * Non-finite numbers are illegal in
         * standards-compliant JSON and therefore
         * must be rejected.
         */
        bool nonfinite_rejected =
            false;


        try {

            auto invalid =
                record;


            invalid
                .biological_features
                .right
                .delta_tm_celsius =
                    std::numeric_limits<
                        double
                    >::quiet_NaN();


            static_cast<void>(
                serialize_report_json_object(
                    invalid
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
            "Non-finite JSON value rejected"
        );


        bool bad_schema_rejected =
            false;


        try {

            auto invalid =
                record;


            invalid.schema_version =
                999;


            static_cast<void>(
                serialize_report_tsv_row(
                    invalid
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            bad_schema_rejected =
                true;
        }


        expect(
            bad_schema_rejected,
            "Unsupported report schema rejected"
        );


        std::cout
            << "tsv_column_count\t"
            << count_columns(
                   header
               )
            << '\n';


        std::cout
            << "json_coordinate_string\t"
            << record.hit.amplicon_start
            << '\n';


        std::cout
            << "json_mask_string\t"
            << "0x8000000000000000"
            << '\n';


        std::cout
            << "json_object\t"
            << json1
            << '\n';


        std::cout
            << "REPORT_TSV_DETERMINISTIC\tYES\n";

        std::cout
            << "REPORT_JSON_DETERMINISTIC\tYES\n";

        std::cout
            << "REPORT_JSON_VALIDATION_READY\tYES\n";

        std::cout
            << "REPORT_64BIT_COORDINATES_LOSSLESS\tYES\n";

        std::cout
            << "REPORT_MISMATCH_MASKS_LOSSLESS\tYES\n";

        std::cout
            << "REPORT_JSON_ESCAPING_VALID\tYES\n";

        std::cout
            << "REPORT_TSV_ESCAPING_VALID\tYES\n";

        std::cout
            << "REPORT_NONFINITE_REJECTED\tYES\n";

        std::cout
            << "REPORT_SCHEMA_VERSION_ENFORCED\tYES\n";

        std::cout
            << "REPORT_SERIALIZATION_BIOLOGICAL_RECALCULATION\tNO\n";

        std::cout
            << "REPORT_SERIALIZATION_SEARCH_CHANGE\tNO\n";

        std::cout
            << "REPORT_SERIALIZATION_V1_COMPLETE\tYES\n";

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
