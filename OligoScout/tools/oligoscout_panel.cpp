/*
 * OLIGOSCOUT_PANEL_V1
 *
 * Panel-level whole-genome search entry point.
 *
 * Scope boundary (deliberate):
 *   This tool performs SEARCH and PAIR ASSEMBLY only. It introduces no
 *   thermodynamic formula, no scoring and no calibration, and it has no
 *   Primer3 dependency. Every binding site and every assembled product is
 *   emitted as a record; thermodynamic interpretation is the job of the
 *   QC layer that consumes these records.
 *
 * Differences from integrate_grch38_annotated_report (a benchmark fixture):
 *   - accepts a panel of any size, not exactly 32 pairs
 *   - accepts either primer1/primer2 or forward_primer/reverse_primer headers
 *   - requires no ntthal / oligotm / primer3 config
 *   - emits all off-target products, not only the expected chr22 targets
 *
 * Usage:
 *   oligoscout_panel <manifest.tsv> <ppfm_dir> <panel.tsv> <out_prefix>
 *       [--anchor 12] [--max-mismatch 3]
 *       [--min-amplicon 50] [--max-amplicon 3000]
 *       [--cache 1]
 *
 * Outputs:
 *   <out_prefix>.bindings.tsv
 *   <out_prefix>.products.tsv
 *   <out_prefix>.contract.tsv
 */

#include "primerpair/packed_reference.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"
#include "primerpair/primer_pair_assembler.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <cctype>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

struct PanelPair {
    std::string pair_id;
    std::string primer1;
    std::string primer2;
};


std::vector<std::string>
split_tabs(
    const std::string& line
)
{
    std::vector<std::string> fields;

    std::size_t begin = 0;

    while (true) {

        const std::size_t tab =
            line.find('\t', begin);

        if (tab == std::string::npos) {
            fields.push_back(
                line.substr(begin)
            );
            break;
        }

        fields.push_back(
            line.substr(begin, tab - begin)
        );

        begin = tab + 1;
    }

    return fields;
}


std::string
lowercase(
    std::string value
)
{
    for (char& c : value) {
        c = static_cast<char>(
            std::tolower(
                static_cast<unsigned char>(c)
            )
        );
    }
    return value;
}


std::size_t
find_column(
    const std::vector<std::string>& header,
    const std::vector<std::string>& candidates
)
{
    for (const auto& want : candidates) {
        for (std::size_t i = 0; i < header.size(); ++i) {
            if (lowercase(header[i]) == want) {
                return i;
            }
        }
    }
    return header.size();
}


std::vector<PanelPair>
load_panel(
    const std::filesystem::path& path
)
{
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Could not open primer panel."
        );
    }

    std::string header_line;

    if (!std::getline(input, header_line)) {
        throw std::runtime_error(
            "Primer panel is empty."
        );
    }

    if (!header_line.empty() &&
        header_line.back() == '\r') {
        header_line.pop_back();
    }

    const auto header = split_tabs(header_line);

    const std::size_t id_col =
        find_column(header, {"pair_id", "id", "name"});

    const std::size_t left_col =
        find_column(
            header,
            {"primer1", "forward_primer", "left_primer", "forward"}
        );

    const std::size_t right_col =
        find_column(
            header,
            {"primer2", "reverse_primer", "right_primer", "reverse"}
        );

    if (left_col == header.size() ||
        right_col == header.size()) {
        throw std::runtime_error(
            "Panel header must contain primer1/primer2 or "
            "forward_primer/reverse_primer."
        );
    }

    std::vector<PanelPair> panel;

    std::string line;
    std::size_t row_index = 0;

    while (std::getline(input, line)) {

        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line.empty()) {
            continue;
        }

        const auto fields = split_tabs(line);

        if (fields.size() <= left_col ||
            fields.size() <= right_col) {
            throw std::runtime_error(
                "Malformed primer panel row."
            );
        }

        PanelPair row;

        row.pair_id =
            (id_col < fields.size())
                ? fields[id_col]
                : std::to_string(row_index);

        row.primer1 = fields[left_col];
        row.primer2 = fields[right_col];

        panel.push_back(std::move(row));

        ++row_index;
    }

    if (panel.empty()) {
        throw std::runtime_error(
            "Primer panel has no rows."
        );
    }

    return panel;
}


std::string
render_mask(
    const std::uint64_t mask,
    const std::size_t primer_length
)
{
    /*
     * Bit i (LSB = 0) is a mismatch at primer position i,
     * counted from the biological 5' end.
     */
    std::string rendered(primer_length, '.');

    for (std::size_t i = 0; i < primer_length; ++i) {
        if (mask & (std::uint64_t{1} << i)) {
            rendered[i] = 'X';
        }
    }

    return rendered;
}


std::string
reference_segment(
    const primerpair::PackedReference& reference,
    const std::uint64_t start,
    const std::size_t length
)
{
    std::string segment;
    segment.reserve(length);

    for (std::size_t i = 0; i < length; ++i) {
        const std::uint64_t position = start + i;

        if (position >= reference.size()) {
            segment.push_back('N');
            continue;
        }

        segment.push_back(
            reference.base_at(position)
        );
    }

    return segment;
}


std::uint64_t
parse_option(
    const std::string& value,
    const char* name
)
{
    try {
        return std::stoull(value);
    } catch (...) {
        throw std::runtime_error(
            std::string("Invalid value for ") + name
        );
    }
}

}  // namespace


int main(
    const int argc,
    char** argv
)
{
    try {
        using namespace primerpair;

        if (argc < 5) {
            std::cerr
                << "Usage:\n"
                << argv[0]
                << " <manifest.tsv> <ppfm_dir> <panel.tsv> <out_prefix>\n"
                << "    [--anchor 12] [--max-mismatch 3]\n"
                << "    [--min-amplicon 50] [--max-amplicon 3000]\n"
                << "    [--cache 1]\n";
            return 2;
        }

        const std::filesystem::path manifest_path = argv[1];
        const std::filesystem::path index_directory = argv[2];
        const std::filesystem::path panel_path = argv[3];
        const std::string output_prefix = argv[4];

        std::size_t anchor_length = 12;
        std::size_t max_mismatches = 3;
        std::uint64_t min_amplicon = 50;
        std::uint64_t max_amplicon = 3000;
        std::size_t cache_capacity = 1;

        for (int i = 5; i + 1 < argc; i += 2) {
            const std::string flag = argv[i];
            const std::string value = argv[i + 1];

            if (flag == "--anchor") {
                anchor_length =
                    static_cast<std::size_t>(
                        parse_option(value, "--anchor")
                    );
            } else if (flag == "--max-mismatch") {
                max_mismatches =
                    static_cast<std::size_t>(
                        parse_option(value, "--max-mismatch")
                    );
            } else if (flag == "--min-amplicon") {
                min_amplicon =
                    parse_option(value, "--min-amplicon");
            } else if (flag == "--max-amplicon") {
                max_amplicon =
                    parse_option(value, "--max-amplicon");
            } else if (flag == "--cache") {
                cache_capacity =
                    static_cast<std::size_t>(
                        parse_option(value, "--cache")
                    );
            } else {
                throw std::runtime_error(
                    "Unknown option: " + flag
                );
            }
        }

        if (cache_capacity == 0) {
            cache_capacity = 1;
        }

        auto manifest =
            PpfmManifest::load(
                manifest_path,
                index_directory
            );

        const auto panel = load_panel(panel_path);

        std::ofstream bindings(
            output_prefix + ".bindings.tsv"
        );

        std::ofstream products(
            output_prefix + ".products.tsv"
        );

        if (!bindings || !products) {
            throw std::runtime_error(
                "Could not open output files."
            );
        }

        bindings
            << "pair_id\trole\tprimer\tchrom\tstart\tend_exclusive"
               "\tstrand\tmismatches\tmask_bits\tmask\tsite_plus_strand\n";

        products
            << "pair_id\tchrom\tamplicon_start\tamplicon_end_exclusive"
               "\tamplicon_length\tleft_position\tright_position"
               "\tleft_mismatches\tright_mismatches\ttotal_mismatches"
               "\tleft_mask\tright_mask\n";

        PpfmShardCache cache(
            manifest,
            cache_capacity,
            8
        );

        const auto started =
            std::chrono::steady_clock::now();

        std::uint64_t binding_total = 0;
        std::uint64_t product_total = 0;
        std::uint64_t unique_primer_searches = 0;

        for (const auto& entry : manifest.entries()) {

            const GenomeShard& shard =
                cache.get(entry.chromosome);

            StrandAwarePrimerSearchEngine searcher(
                shard.index(),
                shard.reference()
            );

            std::unordered_map<
                std::string,
                std::vector<OrientedPrimerSearchHit>
            > primer_hits;

            primer_hits.reserve(panel.size() * 2);

            for (const auto& pair : panel) {
                for (const std::string* primer :
                     {&pair.primer1, &pair.primer2}) {

                    if (primer_hits.count(*primer) != 0) {
                        continue;
                    }

                    auto result =
                        searcher.search(
                            *primer,
                            anchor_length,
                            max_mismatches
                        );

                    primer_hits.emplace(
                        *primer,
                        std::move(result.hits)
                    );

                    ++unique_primer_searches;
                }
            }

            /*
             * Binding records. A primer shared by several pairs is
             * reported once per pair so that every pair's record set
             * is self-contained.
             */
            for (const auto& pair : panel) {

                const std::string* sequences[2] = {
                    &pair.primer1,
                    &pair.primer2
                };

                const char* roles[2] = {"F", "R"};

                for (int side = 0; side < 2; ++side) {

                    const auto found =
                        primer_hits.find(*sequences[side]);

                    if (found == primer_hits.end()) {
                        continue;
                    }

                    const std::size_t primer_length =
                        sequences[side]->size();

                    for (const auto& hit : found->second) {

                        bindings
                            << pair.pair_id << '\t'
                            << roles[side] << '\t'
                            << *sequences[side] << '\t'
                            << entry.chromosome << '\t'
                            << hit.position << '\t'
                            << (hit.position + primer_length) << '\t'
                            << (hit.orientation ==
                                PrimerOrientation::Forward ? '+' : '-')
                            << '\t'
                            << hit.mismatches << '\t'
                            << hit.mismatch_mask << '\t'
                            << render_mask(
                                   hit.mismatch_mask,
                                   primer_length
                               )
                            << '\t'
                            << reference_segment(
                                   shard.reference(),
                                   hit.position,
                                   primer_length
                               )
                            << '\n';

                        ++binding_total;
                    }
                }

                const auto left =
                    primer_hits.find(pair.primer1);

                const auto right =
                    primer_hits.find(pair.primer2);

                if (left == primer_hits.end() ||
                    right == primer_hits.end()) {
                    continue;
                }

                const auto assembled =
                    assemble_primer_pair_hits(
                        pair.primer1,
                        left->second,
                        pair.primer2,
                        right->second,
                        min_amplicon,
                        max_amplicon
                    );

                for (const auto& hit : assembled.amplicons) {

                    products
                        << pair.pair_id << '\t'
                        << entry.chromosome << '\t'
                        << hit.amplicon_start << '\t'
                        << hit.amplicon_end_exclusive << '\t'
                        << hit.amplicon_length << '\t'
                        << hit.left_position << '\t'
                        << hit.right_position << '\t'
                        << hit.left_mismatches << '\t'
                        << hit.right_mismatches << '\t'
                        << hit.total_mismatches() << '\t'
                        << render_mask(
                               hit.left_mismatch_mask,
                               pair.primer1.size()
                           )
                        << '\t'
                        << render_mask(
                               hit.right_mismatch_mask,
                               pair.primer2.size()
                           )
                        << '\n';

                    ++product_total;
                }
            }

            std::cerr
                << entry.chromosome
                << " done, bindings="
                << binding_total
                << ", products="
                << product_total
                << '\n';
        }

        const auto elapsed =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started
            ).count();

        std::ofstream contract(
            output_prefix + ".contract.tsv"
        );

        contract
            << "field\tvalue\n"
            << "manifest\t" << manifest_path.string() << '\n'
            << "index_dir\t" << index_directory.string() << '\n'
            << "shards\t" << manifest.size() << '\n'
            << "panel\t" << panel_path.string() << '\n'
            << "panel_pairs\t" << panel.size() << '\n'
            << "anchor_length\t" << anchor_length << '\n'
            << "max_mismatches\t" << max_mismatches << '\n'
            << "min_amplicon_length\t" << min_amplicon << '\n'
            << "max_amplicon_length\t" << max_amplicon << '\n'
            << "mask_convention\tbit i = mismatch at primer position i "
               "from the 5' end\n"
            << "site_strand\tplus strand of the reference\n"
            << "binding_records\t" << binding_total << '\n'
            << "product_records\t" << product_total << '\n'
            << "unique_primer_searches\t" << unique_primer_searches << '\n'
            << "shard_loads\t" << cache.load_count() << '\n'
            << "shard_evictions\t" << cache.eviction_count() << '\n'
            << "search_seconds\t" << elapsed << '\n'
            << "thermodynamics\tnot computed in this tool (QC layer)\n";

        std::cout
            << "bindings\t" << binding_total << '\n'
            << "products\t" << product_total << '\n'
            << "seconds\t" << elapsed << '\n'
            << "outputs\t" << output_prefix << ".{bindings,products,contract}.tsv\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
