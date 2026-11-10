/*
 * REPORT_SERIALIZATION_V1
 */

#include "primerpair/report_serialization.hpp"

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
            "Could not serialize uint64 value."
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
            "Could not serialize size value."
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
            "Could not serialize uint32 value."
        );
    }

    return std::string(
        buffer.data(),
        result.ptr
    );
}


std::string decimal_enum(
    const auto value
)
{
    const auto numeric =
        static_cast<
            std::uint64_t
        >(value);

    return decimal_u64(
        numeric
    );
}


std::string format_double(
    const double value
)
{
    if (!std::isfinite(value)) {

        throw std::invalid_argument(
            "Report serialization does not "
            "accept NaN or infinity."
        );
    }


    std::array<char, 128> buffer{};


    const auto result =
        std::to_chars(
            buffer.data(),
            buffer.data() +
                buffer.size(),
            value,
            std::chars_format::general,
            std::numeric_limits<
                double
            >::max_digits10
        );


    if (
        result.ec !=
        std::errc{}
    ) {

        throw std::runtime_error(
            "Could not serialize floating "
            "point value."
        );
    }


    return std::string(
        buffer.data(),
        result.ptr
    );
}


std::string bool_text(
    const bool value
)
{
    return value
        ? "true"
        : "false";
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
            "Could not serialize mismatch mask."
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


const char* primer_name(
    const PrimerIdentity identity
)
{
    switch (identity) {

        case PrimerIdentity::Primer1:
            return "PRIMER1";

        case PrimerIdentity::Primer2:
            return "PRIMER2";
    }

    throw std::invalid_argument(
        "Unknown PrimerIdentity."
    );
}


std::string optional_size_tsv(
    const std::optional<
        std::size_t
    >& value
)
{
    if (!value.has_value()) {
        return "";
    }

    return decimal_size(
        *value
    );
}


std::string optional_size_json(
    const std::optional<
        std::size_t
    >& value
)
{
    if (!value.has_value()) {
        return "null";
    }

    return decimal_size(
        *value
    );
}


std::string tsv_escape(
    const std::string_view value
)
{
    std::string output;

    output.reserve(
        value.size()
    );


    for (const char character : value) {

        switch (character) {

            case '\\':
                output += "\\\\";
                break;

            case '\t':
                output += "\\t";
                break;

            case '\n':
                output += "\\n";
                break;

            case '\r':
                output += "\\r";
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


std::string json_escape(
    const std::string_view value
)
{
    static constexpr char hex[] =
        "0123456789abcdef";


    std::string output;

    output.reserve(
        value.size() +
        2
    );


    output.push_back('"');


    for (
        const unsigned char character :
        value
    ) {

        switch (character) {

            case '"':
                output += "\\\"";
                break;

            case '\\':
                output += "\\\\";
                break;

            case '\b':
                output += "\\b";
                break;

            case '\f':
                output += "\\f";
                break;

            case '\n':
                output += "\\n";
                break;

            case '\r':
                output += "\\r";
                break;

            case '\t':
                output += "\\t";
                break;

            default:

                if (character < 0x20) {

                    output += "\\u00";

                    output.push_back(
                        hex[
                            (
                                character >>
                                4
                            ) &
                            0x0f
                        ]
                    );

                    output.push_back(
                        hex[
                            character &
                            0x0f
                        ]
                    );

                } else {

                    output.push_back(
                        static_cast<char>(
                            character
                        )
                    );
                }

                break;
        }
    }


    output.push_back('"');

    return output;
}


void append_tsv_field(
    std::string& output,
    const std::string_view value,
    const bool first
)
{
    if (!first) {
        output.push_back('\t');
    }

    output += value;
}


void append_json_key(
    std::string& output,
    const std::string_view key,
    const bool first
)
{
    if (!first) {
        output.push_back(',');
    }

    output += json_escape(
        key
    );

    output.push_back(':');
}


void append_json_string_field(
    std::string& output,
    const std::string_view key,
    const std::string_view value,
    const bool first
)
{
    append_json_key(
        output,
        key,
        first
    );

    output += json_escape(
        value
    );
}


void append_json_raw_field(
    std::string& output,
    const std::string_view key,
    const std::string_view value,
    const bool first
)
{
    append_json_key(
        output,
        key,
        first
    );

    output += value;
}


void append_json_u64_string_field(
    std::string& output,
    const std::string_view key,
    const std::uint64_t value,
    const bool first
)
{
    append_json_string_field(
        output,
        key,
        decimal_u64(value),
        first
    );
}


void append_json_mask_field(
    std::string& output,
    const std::string_view key,
    const std::uint64_t value,
    const bool first
)
{
    append_json_string_field(
        output,
        key,
        mask_hex(value),
        first
    );
}


void validate_serializable_record(
    const PrimerPairReportRecord& record
)
{
    if (
        record.schema_version !=
        kPrimerPairReportSchemaVersion
    ) {

        throw std::invalid_argument(
            "Unsupported report schema version."
        );
    }


    if (record.chromosome.empty()) {

        throw std::invalid_argument(
            "Report chromosome must not be empty."
        );
    }


    /*
     * Float fields emitted by V1.
     *
     * format_double() performs the final
     * finite-value validation as well.
     */
    static_cast<void>(
        format_double(
            record
                .biological_features
                .left
                .mismatch_fraction
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .right
                .mismatch_fraction
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .left
                .normalized_3prime_positional_burden
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .right
                .normalized_3prime_positional_burden
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .left
                .perfect_match_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .right
                .perfect_match_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .left
                .observed_binding_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .right
                .observed_binding_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .left
                .delta_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .right
                .delta_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .left
                .oligo_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .right
                .oligo_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .left
                .hairpin_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_features
                .right
                .hairpin_tm_celsius
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_risk
                .left
                .uncalibrated_ranking_score
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_risk
                .right
                .uncalibrated_ranking_score
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_risk
                .primer_min_score
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_risk
                .primer_mean_score
        )
    );

    static_cast<void>(
        format_double(
            record
                .biological_risk
                .uncalibrated_ranking_score
        )
    );

    static_cast<void>(
        format_double(
            record
                .chemistry
                .annealing_offset_celsius
        )
    );
}


}  // namespace


std::string report_tsv_header()
{
    return
        "schema_version"
        "\tchromosome"
        "\tamplicon_start"
        "\tamplicon_end_exclusive"
        "\tamplicon_length"
        "\tleft_primer"
        "\tright_primer"
        "\tleft_position"
        "\tright_position"
        "\tleft_mismatch_count"
        "\tright_mismatch_count"
        "\tleft_mismatch_mask"
        "\tright_mismatch_mask"
        "\tleft_reverse"
        "\tright_reverse"
        "\tleft_mismatch_fraction"
        "\tright_mismatch_fraction"
        "\tleft_3prime_burden"
        "\tright_3prime_burden"
        "\tleft_exact_3prime_run"
        "\tright_exact_3prime_run"
        "\tleft_nearest_mismatch_3prime"
        "\tright_nearest_mismatch_3prime"
        "\tleft_perfect_tm_c"
        "\tright_perfect_tm_c"
        "\tleft_observed_tm_c"
        "\tright_observed_tm_c"
        "\tleft_delta_tm_c"
        "\tright_delta_tm_c"
        "\tleft_oligo_tm_c"
        "\tright_oligo_tm_c"
        "\tleft_hairpin_tm_c"
        "\tright_hairpin_tm_c"
        "\tleft_uncalibrated_score"
        "\tright_uncalibrated_score"
        "\tpair_min_score"
        "\tpair_mean_score"
        "\tpair_uncalibrated_score"
        "\tscore_semantics"
        "\tchemistry_kind_id"
        "\tchemistry_name"
        "\tannealing_rule_id"
        "\tevidence_source_id"
        "\tannealing_offset_c"
        "\trequires_chemistry_specific_tm"
        "\tsource_backed"
        "\tranking_score_calibrated"
        "\tempirical_calibration_applied"
        "\tpcr_probability_available"
        "\tsource_search_hit_preserved";
}


std::string serialize_report_tsv_row(
    const PrimerPairReportRecord& record
)
{
    validate_serializable_record(
        record
    );


    std::string output;

    output.reserve(
        1024
    );


    bool first = true;


    const auto append =
        [&output, &first](
            const std::string_view value
        ) {

            append_tsv_field(
                output,
                value,
                first
            );

            first = false;
        };


    append(
        decimal_u32(
            record.schema_version
        )
    );

    append(
        tsv_escape(
            record.chromosome
        )
    );

    append(
        decimal_u64(
            record.hit.amplicon_start
        )
    );

    append(
        decimal_u64(
            record.hit
                .amplicon_end_exclusive
        )
    );

    append(
        decimal_u64(
            record.hit.amplicon_length
        )
    );

    append(
        primer_name(
            record.hit.left_primer
        )
    );

    append(
        primer_name(
            record.hit.right_primer
        )
    );

    append(
        decimal_u64(
            record.hit.left_position
        )
    );

    append(
        decimal_u64(
            record.hit.right_position
        )
    );

    append(
        decimal_size(
            record.hit.left_mismatches
        )
    );

    append(
        decimal_size(
            record.hit.right_mismatches
        )
    );

    append(
        mask_hex(
            record.hit.left_mismatch_mask
        )
    );

    append(
        mask_hex(
            record.hit.right_mismatch_mask
        )
    );

    append(
        bool_text(
            record
                .biological_features
                .left
                .reverse_strand
        )
    );

    append(
        bool_text(
            record
                .biological_features
                .right
                .reverse_strand
        )
    );

    append(
        format_double(
            record
                .biological_features
                .left
                .mismatch_fraction
        )
    );

    append(
        format_double(
            record
                .biological_features
                .right
                .mismatch_fraction
        )
    );

    append(
        format_double(
            record
                .biological_features
                .left
                .normalized_3prime_positional_burden
        )
    );

    append(
        format_double(
            record
                .biological_features
                .right
                .normalized_3prime_positional_burden
        )
    );

    append(
        decimal_size(
            record
                .biological_features
                .left
                .exact_3prime_run_length
        )
    );

    append(
        decimal_size(
            record
                .biological_features
                .right
                .exact_3prime_run_length
        )
    );

    append(
        optional_size_tsv(
            record
                .biological_features
                .left
                .nearest_mismatch_to_3prime
        )
    );

    append(
        optional_size_tsv(
            record
                .biological_features
                .right
                .nearest_mismatch_to_3prime
        )
    );

    append(
        format_double(
            record
                .biological_features
                .left
                .perfect_match_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .right
                .perfect_match_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .left
                .observed_binding_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .right
                .observed_binding_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .left
                .delta_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .right
                .delta_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .left
                .oligo_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .right
                .oligo_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .left
                .hairpin_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_features
                .right
                .hairpin_tm_celsius
        )
    );

    append(
        format_double(
            record
                .biological_risk
                .left
                .uncalibrated_ranking_score
        )
    );

    append(
        format_double(
            record
                .biological_risk
                .right
                .uncalibrated_ranking_score
        )
    );

    append(
        format_double(
            record
                .biological_risk
                .primer_min_score
        )
    );

    append(
        format_double(
            record
                .biological_risk
                .primer_mean_score
        )
    );

    append(
        format_double(
            record
                .biological_risk
                .uncalibrated_ranking_score
        )
    );

    append(
        "UNCALIBRATED_RANKING"
    );

    append(
        decimal_enum(
            record.chemistry.kind
        )
    );

    append(
        tsv_escape(
            record.chemistry.name
        )
    );

    append(
        decimal_enum(
            record
                .chemistry
                .annealing_rule
        )
    );

    append(
        decimal_enum(
            record
                .chemistry
                .evidence_source
        )
    );

    append(
        format_double(
            record
                .chemistry
                .annealing_offset_celsius
        )
    );

    append(
        bool_text(
            record
                .chemistry
                .requires_chemistry_specific_tm
        )
    );

    append(
        bool_text(
            record
                .chemistry
                .source_backed
        )
    );

    append(
        bool_text(
            record
                .calibration
                .ranking_score_calibrated
        )
    );

    append(
        bool_text(
            record
                .calibration
                .empirical_calibration_applied
        )
    );

    append(
        bool_text(
            record
                .calibration
                .pcr_probability_available
        )
    );

    append(
        bool_text(
            record
                .source_search_hit_preserved
        )
    );


    return output;
}


std::string serialize_report_json_object(
    const PrimerPairReportRecord& record
)
{
    validate_serializable_record(
        record
    );


    std::string output;

    output.reserve(
        1400
    );

    output.push_back('{');


    bool first = true;


    const auto string_field =
        [&output, &first](
            const std::string_view key,
            const std::string_view value
        ) {

            append_json_string_field(
                output,
                key,
                value,
                first
            );

            first = false;
        };


    const auto raw_field =
        [&output, &first](
            const std::string_view key,
            const std::string_view value
        ) {

            append_json_raw_field(
                output,
                key,
                value,
                first
            );

            first = false;
        };


    const auto u64_field =
        [&output, &first](
            const std::string_view key,
            const std::uint64_t value
        ) {

            append_json_u64_string_field(
                output,
                key,
                value,
                first
            );

            first = false;
        };


    const auto mask_field =
        [&output, &first](
            const std::string_view key,
            const std::uint64_t value
        ) {

            append_json_mask_field(
                output,
                key,
                value,
                first
            );

            first = false;
        };


    raw_field(
        "schema_version",
        decimal_u32(
            record.schema_version
        )
    );

    string_field(
        "chromosome",
        record.chromosome
    );

    u64_field(
        "amplicon_start",
        record.hit.amplicon_start
    );

    u64_field(
        "amplicon_end_exclusive",
        record.hit
            .amplicon_end_exclusive
    );

    u64_field(
        "amplicon_length",
        record.hit.amplicon_length
    );

    string_field(
        "left_primer",
        primer_name(
            record.hit.left_primer
        )
    );

    string_field(
        "right_primer",
        primer_name(
            record.hit.right_primer
        )
    );

    u64_field(
        "left_position",
        record.hit.left_position
    );

    u64_field(
        "right_position",
        record.hit.right_position
    );

    raw_field(
        "left_mismatch_count",
        decimal_size(
            record.hit.left_mismatches
        )
    );

    raw_field(
        "right_mismatch_count",
        decimal_size(
            record.hit.right_mismatches
        )
    );

    mask_field(
        "left_mismatch_mask",
        record.hit.left_mismatch_mask
    );

    mask_field(
        "right_mismatch_mask",
        record.hit.right_mismatch_mask
    );

    raw_field(
        "left_reverse",
        bool_text(
            record
                .biological_features
                .left
                .reverse_strand
        )
    );

    raw_field(
        "right_reverse",
        bool_text(
            record
                .biological_features
                .right
                .reverse_strand
        )
    );

    raw_field(
        "left_mismatch_fraction",
        format_double(
            record
                .biological_features
                .left
                .mismatch_fraction
        )
    );

    raw_field(
        "right_mismatch_fraction",
        format_double(
            record
                .biological_features
                .right
                .mismatch_fraction
        )
    );

    raw_field(
        "left_3prime_burden",
        format_double(
            record
                .biological_features
                .left
                .normalized_3prime_positional_burden
        )
    );

    raw_field(
        "right_3prime_burden",
        format_double(
            record
                .biological_features
                .right
                .normalized_3prime_positional_burden
        )
    );

    raw_field(
        "left_exact_3prime_run",
        decimal_size(
            record
                .biological_features
                .left
                .exact_3prime_run_length
        )
    );

    raw_field(
        "right_exact_3prime_run",
        decimal_size(
            record
                .biological_features
                .right
                .exact_3prime_run_length
        )
    );

    raw_field(
        "left_nearest_mismatch_3prime",
        optional_size_json(
            record
                .biological_features
                .left
                .nearest_mismatch_to_3prime
        )
    );

    raw_field(
        "right_nearest_mismatch_3prime",
        optional_size_json(
            record
                .biological_features
                .right
                .nearest_mismatch_to_3prime
        )
    );

    raw_field(
        "left_perfect_tm_c",
        format_double(
            record
                .biological_features
                .left
                .perfect_match_tm_celsius
        )
    );

    raw_field(
        "right_perfect_tm_c",
        format_double(
            record
                .biological_features
                .right
                .perfect_match_tm_celsius
        )
    );

    raw_field(
        "left_observed_tm_c",
        format_double(
            record
                .biological_features
                .left
                .observed_binding_tm_celsius
        )
    );

    raw_field(
        "right_observed_tm_c",
        format_double(
            record
                .biological_features
                .right
                .observed_binding_tm_celsius
        )
    );

    raw_field(
        "left_delta_tm_c",
        format_double(
            record
                .biological_features
                .left
                .delta_tm_celsius
        )
    );

    raw_field(
        "right_delta_tm_c",
        format_double(
            record
                .biological_features
                .right
                .delta_tm_celsius
        )
    );

    raw_field(
        "left_oligo_tm_c",
        format_double(
            record
                .biological_features
                .left
                .oligo_tm_celsius
        )
    );

    raw_field(
        "right_oligo_tm_c",
        format_double(
            record
                .biological_features
                .right
                .oligo_tm_celsius
        )
    );

    raw_field(
        "left_hairpin_tm_c",
        format_double(
            record
                .biological_features
                .left
                .hairpin_tm_celsius
        )
    );

    raw_field(
        "right_hairpin_tm_c",
        format_double(
            record
                .biological_features
                .right
                .hairpin_tm_celsius
        )
    );

    raw_field(
        "left_uncalibrated_score",
        format_double(
            record
                .biological_risk
                .left
                .uncalibrated_ranking_score
        )
    );

    raw_field(
        "right_uncalibrated_score",
        format_double(
            record
                .biological_risk
                .right
                .uncalibrated_ranking_score
        )
    );

    raw_field(
        "pair_min_score",
        format_double(
            record
                .biological_risk
                .primer_min_score
        )
    );

    raw_field(
        "pair_mean_score",
        format_double(
            record
                .biological_risk
                .primer_mean_score
        )
    );

    raw_field(
        "pair_uncalibrated_score",
        format_double(
            record
                .biological_risk
                .uncalibrated_ranking_score
        )
    );

    string_field(
        "score_semantics",
        "UNCALIBRATED_RANKING"
    );

    raw_field(
        "chemistry_kind_id",
        decimal_enum(
            record.chemistry.kind
        )
    );

    string_field(
        "chemistry_name",
        record.chemistry.name
    );

    raw_field(
        "annealing_rule_id",
        decimal_enum(
            record
                .chemistry
                .annealing_rule
        )
    );

    raw_field(
        "evidence_source_id",
        decimal_enum(
            record
                .chemistry
                .evidence_source
        )
    );

    raw_field(
        "annealing_offset_c",
        format_double(
            record
                .chemistry
                .annealing_offset_celsius
        )
    );

    raw_field(
        "requires_chemistry_specific_tm",
        bool_text(
            record
                .chemistry
                .requires_chemistry_specific_tm
        )
    );

    raw_field(
        "source_backed",
        bool_text(
            record
                .chemistry
                .source_backed
        )
    );

    raw_field(
        "ranking_score_calibrated",
        bool_text(
            record
                .calibration
                .ranking_score_calibrated
        )
    );

    raw_field(
        "empirical_calibration_applied",
        bool_text(
            record
                .calibration
                .empirical_calibration_applied
        )
    );

    raw_field(
        "pcr_probability_available",
        bool_text(
            record
                .calibration
                .pcr_probability_available
        )
    );

    raw_field(
        "source_search_hit_preserved",
        bool_text(
            record
                .source_search_hit_preserved
        )
    );


    output.push_back('}');

    return output;
}


void write_report_tsv(
    std::ostream& output,
    const std::span<
        const PrimerPairReportRecord
    > records
)
{
    output
        << report_tsv_header()
        << '\n';


    for (
        const auto& record :
        records
    ) {

        output
            << serialize_report_tsv_row(
                   record
               )
            << '\n';
    }


    if (!output) {

        throw std::runtime_error(
            "Failed to write TSV report."
        );
    }
}


void write_report_json(
    std::ostream& output,
    const std::span<
        const PrimerPairReportRecord
    > records
)
{
    output
        << "[\n";


    for (
        std::size_t i = 0;
        i < records.size();
        ++i
    ) {

        output
            << "  "
            << serialize_report_json_object(
                   records[i]
               );


        if (
            i + 1 <
            records.size()
        ) {

            output
                << ',';
        }


        output
            << '\n';
    }


    output
        << "]\n";


    if (!output) {

        throw std::runtime_error(
            "Failed to write JSON report."
        );
    }
}


}  // namespace primerpair
