/*
 * ANNOTATED_REPORT_PIPELINE_V1
 * Test #78
 */

#include "primerpair/annotated_report_pipeline.hpp"
#include "primerpair/html_report.hpp"
#include "primerpair/report_serialization.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
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


std::filesystem::path
required_environment_path(
    const char* name
)
{
    const char* value =
        std::getenv(
            name
        );


    if (
        value == nullptr ||
        *value == '\0'
    ) {

        throw std::runtime_error(
            std::string(
                "Missing environment variable: "
            ) +
            name
        );
    }


    return std::filesystem::path(
        value
    );
}


char complement(
    const char base
)
{
    switch (base) {

        case 'A':
            return 'T';

        case 'C':
            return 'G';

        case 'G':
            return 'C';

        case 'T':
            return 'A';

        default:
            throw std::invalid_argument(
                "Non-canonical DNA base."
            );
    }
}


std::string test_reverse_complement(
    const std::string_view sequence
)
{
    std::string output;

    output.reserve(
        sequence.size()
    );


    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {

        output.push_back(
            complement(
                *it
            )
        );
    }


    return output;
}


primerpair::BiologicalRiskScoringConfig
integration_risk_config()
{
    using namespace primerpair;


    BiologicalRiskScoringConfig config;

    config.weights.sequence_similarity = 1.0;
    config.weights.positional_3prime_integrity = 1.0;
    config.weights.exact_3prime_fraction = 1.0;
    config.weights.terminal_3_integrity = 1.0;
    config.weights.terminal_5_integrity = 1.0;
    config.weights.terminal_8_integrity = 1.0;
    config.weights.terminal_12_integrity = 1.0;
    config.weights.thermodynamic_retention = 1.0;

    /*
     * Integration fixture only.
     *
     * Not a production biological calibration.
     */
    config.delta_tm_scale_celsius =
        10.0;


    config.weakest_primer_fraction =
        0.75;


    return config;
}


}  // namespace


int main()
{
    try {

        using namespace primerpair;


        const std::string primer1 =
            "GCGTACGATCGTACGCATGC";


        const std::string primer2 =
            "ATGCCGTAACGTTAGCGTCA";


        constexpr std::uint64_t
            left_position = 10;


        constexpr std::uint64_t
            right_position = 190;


        constexpr std::uint64_t
            amplicon_end = 210;


        constexpr std::uint64_t
            amplicon_length = 200;


        std::string reference_sequence(
            240,
            'A'
        );


        reference_sequence.replace(
            static_cast<std::size_t>(
                left_position
            ),
            primer1.size(),
            primer1
        );


        reference_sequence.replace(
            static_cast<std::size_t>(
                right_position
            ),
            primer2.size(),
            test_reverse_complement(
                primer2
            )
        );


        PackedReference reference(
            reference_sequence
        );


        PrimerPairHit hit;


        hit.left_primer =
            PrimerIdentity::Primer1;


        hit.right_primer =
            PrimerIdentity::Primer2;


        hit.left_position =
            left_position;


        hit.right_position =
            right_position;


        hit.left_mismatches =
            0;


        hit.right_mismatches =
            0;


        hit.left_mismatch_mask =
            0;


        hit.right_mismatch_mask =
            0;


        hit.amplicon_start =
            left_position;


        hit.amplicon_end_exclusive =
            amplicon_end;


        hit.amplicon_length =
            amplicon_length;


        Primer3ThermodynamicBackend backend(
            required_environment_path(
                "PRIMERPAIR_NTTHAL"
            ),
            required_environment_path(
                "PRIMERPAIR_OLIGOTM"
            ),
            required_environment_path(
                "PRIMERPAIR_THERMO_CONFIG"
            )
        );


        const auto risk_config =
            integration_risk_config();


        const auto chemistry =
            standard_taq_chemistry_profile();


        const auto record =
            build_annotated_primer_pair_report_record(
                reference,
                "chrSynthetic",
                hit,
                primer1,
                primer2,
                backend,
                risk_config,
                chemistry
            );


        expect(
            record.chromosome ==
                "chrSynthetic",
            "Pipeline chromosome preserved"
        );


        expect(
            record.hit ==
                hit,
            "Pipeline source hit preserved"
        );


        expect(
            record
                .biological_features
                .amplicon_start ==
                hit.amplicon_start &&
            record
                .biological_features
                .amplicon_end_exclusive ==
                hit.amplicon_end_exclusive &&
            record
                .biological_features
                .amplicon_length ==
                hit.amplicon_length,
            "Pipeline biological geometry preserved"
        );


        expect(
            record
                .biological_features
                .left
                .mismatch_count == 0 &&
            record
                .biological_features
                .right
                .mismatch_count == 0,
            "Pipeline exact binding mismatch counts valid"
        );


        expect(
            std::abs(
                record
                    .biological_features
                    .left
                    .delta_tm_celsius
            ) <
                1.0e-12 &&
            std::abs(
                record
                    .biological_features
                    .right
                    .delta_tm_celsius
            ) <
                1.0e-12,
            "Pipeline exact binding delta-Tm zero"
        );


        expect(
            record
                .biological_risk
                .uncalibrated_ranking_score >=
                0.0 &&
            record
                .biological_risk
                .uncalibrated_ranking_score <=
                1.0,
            "Pipeline uncalibrated score bounded"
        );


        expect(
            record.score_semantics ==
                ReportScoreSemantics::
                    UncalibratedRanking,
            "Pipeline score semantics uncalibrated"
        );


        expect(
            !record
                .calibration
                .ranking_score_calibrated &&
            !record
                .calibration
                .empirical_calibration_applied &&
            !record
                .calibration
                .pcr_probability_available,
            "Pipeline calibration flags remain false"
        );


        expect(
            record
                .source_search_hit_preserved,
            "Pipeline marks source search hit preserved"
        );


        const std::string tsv =
            serialize_report_tsv_row(
                record
            );


        const std::string json =
            serialize_report_json_object(
                record
            );


        const std::vector<
            PrimerPairReportRecord
        > records{
            record
        };


        const std::string html =
            render_html_report(
                records
            );


        expect(
            tsv.find(
                "chrSynthetic"
            ) !=
                std::string::npos,
            "Pipeline TSV serialization valid"
        );


        expect(
            json.find(
                "\"chromosome\":\"chrSynthetic\""
            ) !=
                std::string::npos,
            "Pipeline JSON serialization valid"
        );


        expect(
            html.find(
                "UNCALIBRATED RANKING SCORE"
            ) !=
                std::string::npos,
            "Pipeline HTML semantics valid"
        );


        bool empty_chromosome_rejected =
            false;


        try {

            static_cast<void>(
                build_annotated_primer_pair_report_record(
                    reference,
                    "",
                    hit,
                    primer1,
                    primer2,
                    backend,
                    risk_config,
                    chemistry
                )
            );

        } catch (
            const std::invalid_argument&
        ) {

            empty_chromosome_rejected =
                true;
        }


        expect(
            empty_chromosome_rejected,
            "Pipeline empty chromosome rejected"
        );


        std::cout
            << "pipeline_amplicon_length\t"
            << record.hit.amplicon_length
            << '\n';


        std::cout
            << "pipeline_score\t"
            << record
                   .biological_risk
                   .uncalibrated_ranking_score
            << '\n';


        std::cout
            << "ANNOTATED_PIPELINE_FEATURES_VALID\tYES\n";

        std::cout
            << "ANNOTATED_PIPELINE_RISK_SCORE_VALID\tYES\n";

        std::cout
            << "ANNOTATED_PIPELINE_REPORT_MODEL_VALID\tYES\n";

        std::cout
            << "ANNOTATED_PIPELINE_TSV_VALID\tYES\n";

        std::cout
            << "ANNOTATED_PIPELINE_JSON_VALID\tYES\n";

        std::cout
            << "ANNOTATED_PIPELINE_HTML_VALID\tYES\n";

        std::cout
            << "ANNOTATED_PIPELINE_SCORE_UNCALIBRATED\tYES\n";

        std::cout
            << "ANNOTATED_PIPELINE_PCR_PROBABILITY\tNO\n";

        std::cout
            << "ANNOTATED_PIPELINE_SEARCH_CHANGE\tNO\n";

        std::cout
            << "ANNOTATED_REPORT_PIPELINE_V1_COMPLETE\tYES\n";

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
