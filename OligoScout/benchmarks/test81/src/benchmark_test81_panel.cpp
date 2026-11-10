#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"
#include "primerpair/primer_pair_assembler.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

struct PanelPair {
    std::uint64_t pair_id{0};
    std::string primer1_id;
    std::string primer1;
    std::string primer2_id;
    std::string primer2;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        if (end == std::string::npos) {
            fields.emplace_back(line.substr(begin));
            break;
        }
        fields.emplace_back(line.substr(begin, end - begin));
        begin = end + 1;
    }
    return fields;
}

std::size_t column_index(
    const std::vector<std::string>& header,
    const std::string_view name
) {
    for (std::size_t index = 0; index < header.size(); ++index) {
        if (header.at(index) == name) {
            return index;
        }
    }
    throw std::runtime_error("Missing panel column: " + std::string{name});
}

std::vector<PanelPair> read_panel(const std::string& path) {
    std::ifstream input{path};
    if (!input) {
        throw std::runtime_error("Cannot open panel TSV: " + path);
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("Panel TSV is empty.");
    }
    const auto header = split_tab(line);
    const auto pair_id_col = column_index(header, "pair_id");
    const auto p1_id_col = column_index(header, "primer1_id");
    const auto p1_col = column_index(header, "primer1");
    const auto p2_id_col = column_index(header, "primer2_id");
    const auto p2_col = column_index(header, "primer2");
    const std::size_t required =
        std::max({pair_id_col, p1_id_col, p1_col, p2_id_col, p2_col}) + 1;

    std::vector<PanelPair> rows;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        const auto fields = split_tab(line);
        if (fields.size() < required) {
            throw std::runtime_error("Malformed panel row.");
        }
        PanelPair row{
            std::stoull(fields.at(pair_id_col)),
            fields.at(p1_id_col),
            fields.at(p1_col),
            fields.at(p2_id_col),
            fields.at(p2_col),
        };
        if (row.primer1.empty() || row.primer2.empty()) {
            throw std::runtime_error("Empty primer sequence.");
        }
        if (row.primer1.size() > 64 || row.primer2.size() > 64) {
            throw std::runtime_error("Primer length exceeds 64 nt.");
        }
        rows.push_back(std::move(row));
    }
    if (rows.empty()) {
        throw std::runtime_error("Panel contains no pairs.");
    }
    return rows;
}

std::string mask_hex(const std::uint64_t value) {
    std::ostringstream stream;
    stream << "0x" << std::hex << std::nouppercase
           << std::setw(16) << std::setfill('0') << value;
    return stream.str();
}

std::size_t masked_region_count(
    const std::uint64_t mask,
    const std::size_t primer_length,
    const std::size_t region_length
) {
    const std::size_t begin =
        primer_length > region_length ? primer_length - region_length : 0;
    std::size_t count = 0;
    for (std::size_t index = begin; index < primer_length; ++index) {
        count += static_cast<std::size_t>((mask >> index) & 1U);
    }
    return count;
}

std::string reference_sequence(
    const primerpair::PackedReference& reference,
    const std::uint64_t start,
    const std::size_t length
) {
    std::string sequence;
    sequence.reserve(length);
    for (std::size_t index = 0; index < length; ++index) {
        sequence.push_back(reference.base_at(start + index));
    }
    return sequence;
}

const char* identity_name(const primerpair::PrimerIdentity identity) {
    return identity == primerpair::PrimerIdentity::Primer1 ? "P1" : "P2";
}

void write_bindings(
    std::ofstream& output,
    const PanelPair& pair,
    const std::string& chromosome,
    const std::string_view identity,
    const std::string& query_id,
    const std::string& primer,
    const primerpair::StrandAwarePrimerSearchResult& result,
    const primerpair::PackedReference& reference,
    const std::size_t anchor_length,
    std::uint64_t& binding_count
) {
    for (const auto& hit : result.hits) {
        const bool forward =
            hit.orientation == primerpair::PrimerOrientation::Forward;
        const std::string genomic =
            reference_sequence(reference, hit.position, primer.size());
        const std::string binding =
            forward ? genomic : primerpair::reverse_complement(genomic);
        const std::size_t anchor_mismatches =
            masked_region_count(hit.mismatch_mask, primer.size(), anchor_length);
        output
            << pair.pair_id << '\t'
            << identity << '\t'
            << query_id << '\t'
            << chromosome << '\t'
            << (forward ? "plus" : "minus") << '\t'
            << hit.position << '\t'
            << hit.position + primer.size() << '\t'
            << primer.size() << '\t'
            << hit.mismatches << '\t'
            << mask_hex(hit.mismatch_mask) << '\t'
            << anchor_length << '\t'
            << anchor_mismatches << '\t'
            << masked_region_count(hit.mismatch_mask, primer.size(), 3) << '\t'
            << masked_region_count(hit.mismatch_mask, primer.size(), 5) << '\t'
            << primer << '\t'
            << binding << '\n';
        ++binding_count;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 10) {
            std::cerr
                << "Usage: " << argv[0]
                << " <manifest.tsv> <ppfm_directory> <panel.tsv>"
                << " <bindings.tsv> <products.tsv> <anchor_length>"
                << " <max_mismatches> <min_amplicon> <max_amplicon>\n";
            return 2;
        }

        const std::string manifest_path{argv[1]};
        const std::string ppfm_directory{argv[2]};
        const std::string panel_path{argv[3]};
        const std::string binding_path{argv[4]};
        const std::string product_path{argv[5]};
        const std::size_t anchor_length = std::stoull(argv[6]);
        const std::size_t max_mismatches = std::stoull(argv[7]);
        const std::uint64_t min_amplicon = std::stoull(argv[8]);
        const std::uint64_t max_amplicon = std::stoull(argv[9]);

        if (anchor_length == 0 || max_amplicon < min_amplicon) {
            throw std::invalid_argument("Invalid benchmark contract.");
        }

        auto panel = read_panel(panel_path);
        auto manifest = primerpair::PpfmManifest::load(
            manifest_path,
            ppfm_directory
        );
        std::vector<std::string> chromosomes;
        chromosomes.reserve(manifest.size());
        for (const auto& entry : manifest.entries()) {
            chromosomes.push_back(entry.chromosome);
        }
        primerpair::PpfmShardCache cache{std::move(manifest), 1, 8};

        std::ofstream bindings{binding_path};
        std::ofstream products{product_path};
        if (!bindings || !products) {
            throw std::runtime_error("Cannot open Test #81 output TSV.");
        }

        bindings
            << "pair_id\tprimer_identity\tquery_id\tchromosome\tstrand\t"
            << "binding_start\tbinding_end_exclusive\tprimer_length\t"
            << "mismatch_count\tmismatch_mask\tanchor_length\t"
            << "anchor_mismatches\tlast3_mismatches\tlast5_mismatches\t"
            << "primer_sequence\tbinding_sequence\n";

        products
            << "pair_id\tchromosome\tleft_primer\tright_primer\t"
            << "left_position\tright_position\tamplicon_start\t"
            << "amplicon_end_exclusive\tamplicon_length\tleft_mismatches\t"
            << "right_mismatches\ttotal_mismatches\tleft_anchor_mismatches\t"
            << "right_anchor_mismatches\tleft_last3_mismatches\t"
            << "right_last3_mismatches\tleft_mismatch_mask\t"
            << "right_mismatch_mask\n";

        std::uint64_t binding_count = 0;
        std::uint64_t product_count = 0;
        std::uint64_t search_instances = 0;

        for (const auto& chromosome : chromosomes) {
            const auto& shard = cache.get(chromosome);
            primerpair::StrandAwarePrimerSearchEngine searcher{
                shard.index(), shard.reference()
            };
            for (const auto& pair : panel) {
                if (anchor_length > pair.primer1.size() ||
                    anchor_length > pair.primer2.size()) {
                    throw std::invalid_argument("Anchor exceeds primer length.");
                }
                const auto p1 = searcher.search(
                    pair.primer1, anchor_length, max_mismatches
                );
                const auto p2 = searcher.search(
                    pair.primer2, anchor_length, max_mismatches
                );
                search_instances += 2;

                write_bindings(
                    bindings, pair, chromosome, "P1", pair.primer1_id,
                    pair.primer1, p1, shard.reference(), anchor_length,
                    binding_count
                );
                write_bindings(
                    bindings, pair, chromosome, "P2", pair.primer2_id,
                    pair.primer2, p2, shard.reference(), anchor_length,
                    binding_count
                );

                const auto assembled = primerpair::assemble_primer_pair_hits(
                    pair.primer1, p1.hits,
                    pair.primer2, p2.hits,
                    min_amplicon, max_amplicon
                );
                for (const auto& hit : assembled.amplicons) {
                    const std::size_t left_length =
                        hit.left_primer == primerpair::PrimerIdentity::Primer1
                            ? pair.primer1.size() : pair.primer2.size();
                    const std::size_t right_length =
                        hit.right_primer == primerpair::PrimerIdentity::Primer1
                            ? pair.primer1.size() : pair.primer2.size();
                    products
                        << pair.pair_id << '\t'
                        << chromosome << '\t'
                        << identity_name(hit.left_primer) << '\t'
                        << identity_name(hit.right_primer) << '\t'
                        << hit.left_position << '\t'
                        << hit.right_position << '\t'
                        << hit.amplicon_start << '\t'
                        << hit.amplicon_end_exclusive << '\t'
                        << hit.amplicon_length << '\t'
                        << hit.left_mismatches << '\t'
                        << hit.right_mismatches << '\t'
                        << hit.total_mismatches() << '\t'
                        << masked_region_count(
                            hit.left_mismatch_mask, left_length, anchor_length
                        ) << '\t'
                        << masked_region_count(
                            hit.right_mismatch_mask, right_length, anchor_length
                        ) << '\t'
                        << masked_region_count(
                            hit.left_mismatch_mask, left_length, 3
                        ) << '\t'
                        << masked_region_count(
                            hit.right_mismatch_mask, right_length, 3
                        ) << '\t'
                        << mask_hex(hit.left_mismatch_mask) << '\t'
                        << mask_hex(hit.right_mismatch_mask) << '\n';
                    ++product_count;
                }
            }
            std::cerr << "TEST81_CHROMOSOME_COMPLETE\t" << chromosome << '\n';
        }

        bindings.flush();
        products.flush();
        if (!bindings || !products) {
            throw std::runtime_error("Failed while writing Test #81 output.");
        }

        std::cout
            << "PANEL_PAIRS\t" << panel.size() << '\n'
            << "PPFM_SHARDS\t" << chromosomes.size() << '\n'
            << "SEARCH_INSTANCES\t" << search_instances << '\n'
            << "BINDING_COUNT\t" << binding_count << '\n'
            << "PRODUCT_COUNT\t" << product_count << '\n'
            << "ANCHOR_LENGTH\t" << anchor_length << '\n'
            << "MAX_MISMATCHES\t" << max_mismatches << '\n'
            << "AMPLICON_RANGE\t" << min_amplicon << ".." << max_amplicon << '\n'
            << "TEST81_OLIGOSCOUT_COMPLETE\tYES\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
