/*
 * ANNOTATED_REPORT_PIPELINE_V1
 *
 * Whole-GRCh38 production-index integration
 * validation tool.
 *
 * Searches all manifest shards.
 *
 * The real 32-pair chr22 panel is searched
 * across the whole genome.
 *
 * The 32 known intended chr22 products are then
 * annotated through the validated biological
 * feature, scoring and report stack.
 *
 * Equal scoring weights are an integration
 * fixture only and are NOT production-calibrated
 * biological coefficients.
 */

#include "primerpair/annotated_report_pipeline.hpp"
#include "primerpair/html_report.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"
#include "primerpair/primer_pair_assembler.hpp"
#include "primerpair/report_serialization.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {


struct PanelPair {

    std::size_t pair_id{0};

    std::string primer1;

    std::string primer2;

    std::uint64_t left_start{0};

    std::uint64_t right_start{0};

    std::uint64_t
        expected_amplicon_length{0};
};


std::vector<PanelPair>
load_panel(
    const std::filesystem::path& path
)
{
    std::ifstream input(
        path
    );


    if (!input) {

        throw std::runtime_error(
            "Could not open real primer panel."
        );
    }


    std::string header;

    std::getline(
        input,
        header
    );


    if (
        header.find(
            "pair_id"
        ) ==
            std::string::npos
    ) {

        throw std::runtime_error(
            "Unexpected primer-panel header."
        );
    }


    std::vector<PanelPair>
        panel;


    PanelPair row;


    while (
        input
            >> row.pair_id
            >> row.primer1
            >> row.primer2
            >> row.left_start
            >> row.right_start
            >> row.expected_amplicon_length
    ) {

        panel.push_back(
            row
        );
    }


    if (!input.eof()) {

        throw std::runtime_error(
            "Malformed primer panel."
        );
    }


    if (panel.empty()) {

        throw std::runtime_error(
            "Primer panel is empty."
        );
    }


    return panel;
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
     * Deterministic integration fixture.
     *
     * NOT empirically calibrated.
     */
    config.delta_tm_scale_celsius =
        10.0;


    config.weakest_primer_fraction =
        0.75;


    return config;
}


std::filesystem::path
with_extension(
    const std::filesystem::path& prefix,
    const std::string_view extension
)
{
    return std::filesystem::path(
        prefix.string() +
        std::string(extension)
    );
}


}  // namespace


int main(
    const int argc,
    char** argv
)
{
    try {

        using namespace primerpair;


        if (argc != 8) {

            std::cerr
                << "Usage:\n"
                << argv[0]
                << " <manifest.tsv>"
                << " <ppfm_dir>"
                << " <panel.tsv>"
                << " <ntthal>"
                << " <oligotm>"
                << " <primer3_config>"
                << " <output_prefix>\n";

            return 2;
        }


        const std::filesystem::path
            manifest_path =
                argv[1];


        const std::filesystem::path
            index_directory =
                argv[2];


        const std::filesystem::path
            panel_path =
                argv[3];


        const std::filesystem::path
            ntthal_path =
                argv[4];


        const std::filesystem::path
            oligotm_path =
                argv[5];


        const std::filesystem::path
            thermo_config_path =
                argv[6];


        const std::filesystem::path
            output_prefix =
                argv[7];


        auto manifest =
            PpfmManifest::load(
                manifest_path,
                index_directory
            );


        if (manifest.size() != 24) {

            throw std::runtime_error(
                "Expected 24 GRCh38 manifest shards."
            );
        }


        const auto panel =
            load_panel(
                panel_path
            );


        if (panel.size() != 32) {

            throw std::runtime_error(
                "Expected exactly 32 real panel pairs."
            );
        }


        Primer3ThermodynamicBackend backend(
            ntthal_path,
            oligotm_path,
            thermo_config_path
        );


        const auto risk_config =
            integration_risk_config();


        const auto chemistry =
            standard_taq_chemistry_profile();


        PpfmShardCache cache(
            manifest,
            1,
            8
        );


        std::vector<
            PrimerPairReportRecord
        > report_records;


        report_records.reserve(
            panel.size()
        );


        std::uint64_t
            whole_genome_intended_amplicons =
                0;


        std::size_t
            searched_unique_primer_total =
                0;


        std::size_t
            expected_targets_found =
                0;


        for (
            const auto& entry :
            manifest.entries()
        ) {

            const GenomeShard& shard =
                cache.get(
                    entry.chromosome
                );


            StrandAwarePrimerSearchEngine
                searcher(
                    shard.index(),
                    shard.reference()
                );


            std::unordered_map<
                std::string,
                std::vector<
                    OrientedPrimerSearchHit
                >
            > primer_hits;


            primer_hits.reserve(
                panel.size() * 2
            );


            for (
                const auto& pair :
                panel
            ) {

                for (
                    const std::string* primer :
                    {
                        &pair.primer1,
                        &pair.primer2
                    }
                ) {

                    if (
                        primer_hits.find(
                            *primer
                        ) !=
                        primer_hits.end()
                    ) {

                        continue;
                    }


                    auto result =
                        searcher.search(
                            *primer,
                            12,
                            3
                        );


                    primer_hits.emplace(
                        *primer,
                        std::move(
                            result.hits
                        )
                    );
                }
            }


            searched_unique_primer_total +=
                primer_hits.size();


            for (
                const auto& pair :
                panel
            ) {

                const auto p1 =
                    primer_hits.find(
                        pair.primer1
                    );


                const auto p2 =
                    primer_hits.find(
                        pair.primer2
                    );


                if (
                    p1 ==
                        primer_hits.end() ||
                    p2 ==
                        primer_hits.end()
                ) {

                    throw std::runtime_error(
                        "Internal primer-hit cache failure."
                    );
                }


                const auto assembled =
                    assemble_primer_pair_hits(
                        pair.primer1,
                        p1->second,
                        pair.primer2,
                        p2->second,
                        50,
                        3000
                    );


                whole_genome_intended_amplicons +=
                    static_cast<std::uint64_t>(
                        assembled
                            .amplicons
                            .size()
                    );


                if (
                    entry.chromosome !=
                        "chr22"
                ) {

                    continue;
                }


                const std::uint64_t
                    expected_end =
                        pair.left_start +
                        pair.expected_amplicon_length;


                const PrimerPairHit*
                    expected_hit =
                        nullptr;


                for (
                    const auto& hit :
                    assembled.amplicons
                ) {

                    if (
                        hit.amplicon_start ==
                            pair.left_start &&
                        hit.amplicon_end_exclusive ==
                            expected_end &&
                        hit.amplicon_length ==
                            pair.expected_amplicon_length
                    ) {

                        if (
                            expected_hit !=
                                nullptr
                        ) {

                            throw std::runtime_error(
                                "Duplicate expected chr22 target."
                            );
                        }


                        expected_hit =
                            &hit;
                    }
                }


                if (
                    expected_hit ==
                        nullptr
                ) {

                    throw std::runtime_error(
                        "Expected chr22 target not found for pair " +
                        std::to_string(
                            pair.pair_id
                        )
                    );
                }


                report_records.push_back(
                    build_annotated_primer_pair_report_record(
                        shard.reference(),
                        entry.chromosome,
                        *expected_hit,
                        pair.primer1,
                        pair.primer2,
                        backend,
                        risk_config,
                        chemistry
                    )
                );


                ++expected_targets_found;
            }
        }


        if (
            expected_targets_found !=
                panel.size()
        ) {

            throw std::runtime_error(
                "Not all expected chr22 targets were annotated."
            );
        }


        if (
            report_records.size() !=
                panel.size()
        ) {

            throw std::runtime_error(
                "Unexpected report-record count."
            );
        }


        std::sort(
            report_records.begin(),
            report_records.end(),
            [](
                const PrimerPairReportRecord& lhs,
                const PrimerPairReportRecord& rhs
            ) {

                if (
                    lhs.hit.amplicon_start !=
                    rhs.hit.amplicon_start
                ) {

                    return
                        lhs.hit.amplicon_start <
                        rhs.hit.amplicon_start;
                }


                return
                    lhs.hit.amplicon_end_exclusive <
                    rhs.hit.amplicon_end_exclusive;
            }
        );


        for (
            const auto& record :
            report_records
        ) {

            if (
                record.chromosome !=
                    "chr22"
            ) {

                throw std::runtime_error(
                    "Unexpected chromosome in annotated report."
                );
            }


            if (
                record
                    .biological_risk
                    .uncalibrated_ranking_score <
                        0.0 ||
                record
                    .biological_risk
                    .uncalibrated_ranking_score >
                        1.0
            ) {

                throw std::runtime_error(
                    "Out-of-range integration ranking score."
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

                throw std::runtime_error(
                    "Unsupported calibration claim."
                );
            }


            if (
                !record
                    .source_search_hit_preserved
            ) {

                throw std::runtime_error(
                    "Source-search-hit preservation lost."
                );
            }
        }


        const auto tsv_path =
            with_extension(
                output_prefix,
                ".tsv"
            );


        const auto json_path =
            with_extension(
                output_prefix,
                ".json"
            );


        const auto html_path =
            with_extension(
                output_prefix,
                ".html"
            );


        {
            std::ofstream output(
                tsv_path,
                std::ios::binary
            );


            if (!output) {

                throw std::runtime_error(
                    "Could not open TSV output."
                );
            }


            write_report_tsv(
                output,
                report_records
            );
        }


        {
            std::ofstream output(
                json_path,
                std::ios::binary
            );


            if (!output) {

                throw std::runtime_error(
                    "Could not open JSON output."
                );
            }


            write_report_json(
                output,
                report_records
            );
        }


        {
            std::ofstream output(
                html_path,
                std::ios::binary
            );


            if (!output) {

                throw std::runtime_error(
                    "Could not open HTML output."
                );
            }


            HtmlReportOptions options;


            options.title =
                "OligoScout "
                "GRCh38 Panel-32 "
                "Annotated Integration Report";


            options.display_decimal_places =
                4;


            write_html_report(
                output,
                report_records,
                options
            );
        }


        std::cout
            << "manifest_shards\t"
            << manifest.size()
            << '\n';


        std::cout
            << "panel_pairs\t"
            << panel.size()
            << '\n';


        std::cout
            << "searched_unique_primer_instances\t"
            << searched_unique_primer_total
            << '\n';


        std::cout
            << "whole_genome_intended_amplicons\t"
            << whole_genome_intended_amplicons
            << '\n';


        std::cout
            << "expected_chr22_targets_found\t"
            << expected_targets_found
            << '\n';


        std::cout
            << "report_records\t"
            << report_records.size()
            << '\n';


        std::cout
            << "ppfm_loads\t"
            << cache.load_count()
            << '\n';


        std::cout
            << "ppfm_evictions\t"
            << cache.eviction_count()
            << '\n';


        std::cout
            << "tsv_output\t"
            << tsv_path.string()
            << '\n';


        std::cout
            << "json_output\t"
            << json_path.string()
            << '\n';


        std::cout
            << "html_output\t"
            << html_path.string()
            << '\n';


        std::cout
            << "WHOLE_GRCH38_ALL_24_SHARDS_SEARCHED\tYES\n";

        std::cout
            << "WHOLE_GRCH38_REAL_PANEL32_SEARCHED\tYES\n";

        std::cout
            << "WHOLE_GRCH38_EXPECTED_TARGETS_32_OF_32\tYES\n";

        std::cout
            << "WHOLE_GRCH38_REPORT_FEATURES_FROM_VALIDATED_BUILDER\tYES\n";

        std::cout
            << "WHOLE_GRCH38_SCORE_UNCALIBRATED\tYES\n";

        std::cout
            << "WHOLE_GRCH38_PRODUCTION_WEIGHTS\tNO\n";

        std::cout
            << "WHOLE_GRCH38_PCR_PROBABILITY\tNO\n";

        std::cout
            << "WHOLE_GRCH38_SEARCH_CORRECTNESS_CHANGED\tNO\n";

        std::cout
            << "WHOLE_GRCH38_ANNOTATED_REPORT_INTEGRATION_V1_COMPLETE\tYES\n";

        std::cout
            << "ALL_CHECKS\tYES\n";


        return 0;

    } catch (
        const std::exception& e
    ) {

        std::cerr
            << "INTEGRATION_FAILURE\t"
            << e.what()
            << '\n';

        return 1;
    }
}
