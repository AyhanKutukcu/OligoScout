#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
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

struct Probe {
    std::uint64_t pair_id{0};
    std::string identity;
    std::string query_id;
    std::string primer;
    std::string oriented;
    bool plus{true};
    std::size_t candidate_offset{0};
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const auto end = line.find('\t', begin);
        if (end == std::string::npos) {
            fields.emplace_back(line.substr(begin));
            return fields;
        }
        fields.emplace_back(line.substr(begin, end - begin));
        begin = end + 1;
    }
}

std::size_t column_index(
    const std::vector<std::string>& header,
    const std::string_view name
) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) {
        throw std::runtime_error("Missing panel column: " + std::string{name});
    }
    return static_cast<std::size_t>(found - header.begin());
}

std::string uppercase_dna(std::string sequence) {
    for (char& base : sequence) {
        base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            throw std::runtime_error("Panel primer is not strict A/C/G/T DNA.");
        }
    }
    return sequence;
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
    const auto required =
        std::max({pair_id_col, p1_id_col, p1_col, p2_id_col, p2_col}) + 1;

    std::vector<PanelPair> pairs;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        const auto fields = split_tab(line);
        if (fields.size() < required) {
            throw std::runtime_error("Malformed panel row.");
        }
        pairs.push_back(PanelPair{
            std::stoull(fields.at(pair_id_col)),
            fields.at(p1_id_col),
            uppercase_dna(fields.at(p1_col)),
            fields.at(p2_id_col),
            uppercase_dna(fields.at(p2_col)),
        });
    }
    if (pairs.empty()) {
        throw std::runtime_error("Panel contains no primer pairs.");
    }
    return pairs;
}

char complement(const char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string reverse_complement(const std::string_view sequence) {
    std::string result(sequence.size(), 'N');
    for (std::size_t index = 0; index < sequence.size(); ++index) {
        result.at(sequence.size() - 1 - index) = complement(sequence.at(index));
    }
    return result;
}

bool encode_kmer(const std::string_view sequence, std::uint64_t& value) {
    value = 0;
    for (const char base : sequence) {
        std::uint64_t code = 0;
        switch (base) {
            case 'A': code = 0; break;
            case 'C': code = 1; break;
            case 'G': code = 2; break;
            case 'T': code = 3; break;
            default: return false;
        }
        value = (value << 2U) | code;
    }
    return true;
}

std::string mask_hex(const std::uint64_t value) {
    std::ostringstream stream;
    stream << "0x" << std::hex << std::nouppercase
           << std::setw(16) << std::setfill('0') << value;
    return stream.str();
}

std::size_t region_count(
    const std::uint64_t mask,
    const std::size_t primer_length,
    const std::size_t region_length
) {
    const auto begin =
        primer_length > region_length ? primer_length - region_length : 0;
    std::size_t count = 0;
    for (std::size_t index = begin; index < primer_length; ++index) {
        count += static_cast<std::size_t>((mask >> index) & 1U);
    }
    return count;
}

class FastaReader {
public:
    explicit FastaReader(const std::string& path) : input_{path} {
        if (!input_) {
            throw std::runtime_error("Cannot open reference FASTA: " + path);
        }
    }

    bool next(std::string& name, std::string& sequence) {
        name.clear();
        sequence.clear();
        std::string line;
        if (!pending_header_.empty()) {
            line = std::move(pending_header_);
            pending_header_.clear();
        } else {
            while (std::getline(input_, line) && (line.empty() || line[0] != '>')) {
            }
            if (!input_ && (line.empty() || line[0] != '>')) {
                return false;
            }
        }
        if (line.empty() || line[0] != '>') {
            return false;
        }
        const auto end = line.find_first_of(" \t", 1);
        name = line.substr(1, end == std::string::npos ? end : end - 1);
        while (std::getline(input_, line)) {
            if (!line.empty() && line[0] == '>') {
                pending_header_ = std::move(line);
                break;
            }
            for (const char raw : line) {
                if (!std::isspace(static_cast<unsigned char>(raw))) {
                    sequence.push_back(static_cast<char>(
                        std::toupper(static_cast<unsigned char>(raw))
                    ));
                }
            }
        }
        return true;
    }

private:
    std::ifstream input_;
    std::string pending_header_;
};

void add_probe(
    std::vector<Probe>& probes,
    std::unordered_map<std::uint64_t, std::vector<std::size_t>>& lookup,
    const PanelPair& pair,
    const std::string& identity,
    const std::string& query_id,
    const std::string& primer,
    const std::size_t anchor_length,
    const bool plus
) {
    const std::string oriented = plus ? primer : reverse_complement(primer);
    const std::size_t anchor_offset = plus ? primer.size() - anchor_length : 0;
    std::uint64_t key = 0;
    if (!encode_kmer(
            std::string_view{oriented}.substr(anchor_offset, anchor_length), key
        )) {
        throw std::runtime_error("Cannot encode anchor.");
    }
    const auto index = probes.size();
    probes.push_back(Probe{
        pair.pair_id, identity, query_id, primer, oriented, plus, anchor_offset
    });
    lookup[key].push_back(index);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 6) {
            std::cerr
                << "Usage: " << argv[0]
                << " <reference.fa> <panel.tsv> <bindings.tsv>"
                << " <anchor_length> <max_mismatches>\n";
            return 2;
        }
        const std::string reference_path{argv[1]};
        const std::string panel_path{argv[2]};
        const std::string output_path{argv[3]};
        const std::size_t anchor_length = std::stoull(argv[4]);
        const std::size_t max_mismatches = std::stoull(argv[5]);
        if (anchor_length == 0 || anchor_length > 31) {
            throw std::invalid_argument("Oracle anchor must be in 1..31.");
        }

        const auto panel = read_panel(panel_path);
        std::vector<Probe> probes;
        probes.reserve(panel.size() * 4);
        std::unordered_map<std::uint64_t, std::vector<std::size_t>> lookup;
        for (const auto& pair : panel) {
            if (anchor_length > pair.primer1.size() ||
                anchor_length > pair.primer2.size() ||
                pair.primer1.size() > 64 || pair.primer2.size() > 64) {
                throw std::invalid_argument("Invalid primer or anchor length.");
            }
            add_probe(probes, lookup, pair, "P1", pair.primer1_id,
                      pair.primer1, anchor_length, true);
            add_probe(probes, lookup, pair, "P1", pair.primer1_id,
                      pair.primer1, anchor_length, false);
            add_probe(probes, lookup, pair, "P2", pair.primer2_id,
                      pair.primer2, anchor_length, true);
            add_probe(probes, lookup, pair, "P2", pair.primer2_id,
                      pair.primer2, anchor_length, false);
        }

        std::ofstream output{output_path};
        if (!output) {
            throw std::runtime_error("Cannot open oracle output.");
        }
        output
            << "pair_id\tprimer_identity\tquery_id\tchromosome\tstrand\t"
            << "binding_start\tbinding_end_exclusive\tprimer_length\t"
            << "mismatch_count\tmismatch_mask\tanchor_length\t"
            << "anchor_mismatches\tlast3_mismatches\tlast5_mismatches\t"
            << "primer_sequence\tbinding_sequence\n";

        FastaReader fasta{reference_path};
        std::string chromosome;
        std::string reference;
        std::uint64_t binding_count = 0;
        std::uint64_t reference_bp = 0;
        const std::uint64_t key_mask =
            anchor_length == 32
                ? std::numeric_limits<std::uint64_t>::max()
                : (std::uint64_t{1} << (anchor_length * 2U)) - 1U;

        while (fasta.next(chromosome, reference)) {
            reference_bp += reference.size();
            std::uint64_t rolling = 0;
            std::size_t valid = 0;
            for (std::size_t position = 0; position < reference.size(); ++position) {
                std::uint64_t code = 0;
                switch (reference.at(position)) {
                    case 'A': code = 0; break;
                    case 'C': code = 1; break;
                    case 'G': code = 2; break;
                    case 'T': code = 3; break;
                    default:
                        rolling = 0;
                        valid = 0;
                        continue;
                }
                rolling = ((rolling << 2U) | code) & key_mask;
                ++valid;
                if (valid < anchor_length) {
                    continue;
                }
                const std::size_t anchor_start = position + 1 - anchor_length;
                const auto found = lookup.find(rolling);
                if (found == lookup.end()) {
                    continue;
                }
                for (const auto probe_index : found->second) {
                    const auto& probe = probes.at(probe_index);
                    if (anchor_start < probe.candidate_offset) {
                        continue;
                    }
                    const auto start = anchor_start - probe.candidate_offset;
                    if (start + probe.oriented.size() > reference.size()) {
                        continue;
                    }
                    std::size_t mismatches = 0;
                    std::uint64_t mismatch_mask = 0;
                    for (std::size_t index = 0; index < probe.oriented.size(); ++index) {
                        const char observed = reference.at(start + index);
                        if (observed != probe.oriented.at(index)) {
                            if (observed != 'A' && observed != 'C' &&
                                observed != 'G' && observed != 'T') {
                                mismatches = max_mismatches + 1;
                                break;
                            }
                            ++mismatches;
                            const std::size_t original_index =
                                probe.plus ? index : probe.primer.size() - 1 - index;
                            mismatch_mask |= std::uint64_t{1} << original_index;
                            if (mismatches > max_mismatches) {
                                break;
                            }
                        }
                    }
                    if (mismatches > max_mismatches) {
                        continue;
                    }
                    const std::string genomic = reference.substr(start, probe.primer.size());
                    const std::string binding =
                        probe.plus ? genomic : reverse_complement(genomic);
                    output
                        << probe.pair_id << '\t'
                        << probe.identity << '\t'
                        << probe.query_id << '\t'
                        << chromosome << '\t'
                        << (probe.plus ? "plus" : "minus") << '\t'
                        << start << '\t'
                        << start + probe.primer.size() << '\t'
                        << probe.primer.size() << '\t'
                        << mismatches << '\t'
                        << mask_hex(mismatch_mask) << '\t'
                        << anchor_length << '\t'
                        << region_count(mismatch_mask, probe.primer.size(), anchor_length) << '\t'
                        << region_count(mismatch_mask, probe.primer.size(), 3) << '\t'
                        << region_count(mismatch_mask, probe.primer.size(), 5) << '\t'
                        << probe.primer << '\t'
                        << binding << '\n';
                    ++binding_count;
                }
            }
            std::cerr << "ORACLE_CHROMOSOME_COMPLETE\t" << chromosome
                      << "\t" << reference.size() << '\n';
        }
        output.flush();
        if (!output) {
            throw std::runtime_error("Failed while writing oracle output.");
        }
        std::cout
            << "PANEL_PAIRS\t" << panel.size() << '\n'
            << "ORACLE_PROBES\t" << probes.size() << '\n'
            << "REFERENCE_BP\t" << reference_bp << '\n'
            << "BINDING_COUNT\t" << binding_count << '\n'
            << "ANCHOR_LENGTH\t" << anchor_length << '\n'
            << "MAX_MISMATCHES\t" << max_mismatches << '\n'
            << "TEST81_ORACLE_COMPLETE\tYES\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
