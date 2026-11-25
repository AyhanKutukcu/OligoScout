#include "primerpair/terminal_rule_search.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

namespace primerpair {
namespace {

constexpr char kBases[4] = {'A', 'C', 'G', 'T'};

/*
 * IUPAC code -> 4-bit base set (A=1, C=2, G=4, T=8). A template base is a
 * match when it belongs to the primer position's set; degenerate primers
 * (e.g. RT-PCR consensus primers) are therefore searched natively.
 */
std::uint8_t iupac_mask(const char c) {
    switch (c) {
        case 'A': return 1; case 'C': return 2; case 'G': return 4; case 'T': return 8;
        case 'R': return 1 | 4; case 'Y': return 2 | 8; case 'S': return 2 | 4; case 'W': return 1 | 8;
        case 'K': return 4 | 8; case 'M': return 1 | 2; case 'B': return 2 | 4 | 8;
        case 'D': return 1 | 4 | 8; case 'H': return 1 | 2 | 8; case 'V': return 1 | 2 | 4;
        case 'N': return 15;
        default: return 0;
    }
}

std::uint8_t complement_mask(const std::uint8_t m) {
    // A<->T (1<->8), C<->G (2<->4)
    return static_cast<std::uint8_t>(((m & 1) << 3) | ((m & 8) >> 3) | ((m & 2) << 1) | ((m & 4) >> 1));
}

std::uint8_t base_bit(const char b) {
    switch (b) {
        case 'A': return 1; case 'C': return 2; case 'G': return 4; case 'T': return 8;
        default: return 0;
    }
}

std::string normalize(std::string_view primer) {
    std::string s;
    s.reserve(primer.size());
    for (const char c : primer) {
        const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (iupac_mask(u) == 0) {
            throw std::invalid_argument("terminal-rule search accepts A/C/G/T and IUPAC codes only");
        }
        s.push_back(u);
    }
    return s;
}

struct Node {
    BidirectionalInterval state;
    std::size_t depth;          // bases matched from the primer's 3' end
    std::size_t mismatches;
    std::size_t terminal;
};

/*
 * Re-derive the site's mismatch mask from the reference. Returns false when
 * the window hits N, runs off the sequence, or violates the rule.
 */
bool verify(const PackedReference& reference, const std::string& primer, const std::uint64_t pos,
            const bool reverse, const TerminalRuleOptions& o, OrientedPrimerSearchHit& out) {
    const std::size_t L = primer.size();
    if (pos + L > reference.size()) {
        return false;
    }
    std::uint64_t mask = 0;
    std::size_t mm = 0, term = 0;
    for (std::size_t j = 0; j < L; ++j) {
        const char obs = reference.base_at(pos + j);
        if (obs == 'N') {
            return false;
        }
        // primer position (5' origin) aligned to window position j
        const std::size_t p = reverse ? L - 1 - j : j;
        const std::uint8_t allowed = reverse ? complement_mask(iupac_mask(primer[p]))
                                             : iupac_mask(primer[p]);
        if ((base_bit(obs) & allowed) == 0) {
            mask |= (std::uint64_t{1} << p);
            ++mm;
            if (p + o.terminal_window >= L) {
                ++term;
            }
        }
    }
    if (mm > o.max_mismatches || term > o.max_terminal_mismatches) {
        return false;
    }
    out.position = pos;
    out.mismatches = mm;
    out.orientation = reverse ? PrimerOrientation::Reverse : PrimerOrientation::Forward;
    out.mismatch_mask = mask;
    return true;
}

}  // namespace

std::vector<OrientedPrimerSearchHit>
search_terminal_rule(const BidirectionalFMIndex& index, const PackedReference& reference,
                     const std::string_view raw_primer, const TerminalRuleOptions& o,
                     TerminalRuleStats* stats) {
    const std::string primer = normalize(raw_primer);
    const std::size_t L = primer.size();
    if (L == 0 || L > 64) {
        throw std::invalid_argument("primer length must be 1..64");
    }
    if (o.max_mismatches > 3) {
        throw std::invalid_argument("at most 3 mismatches per primer are supported");
    }
    const std::size_t window = std::min(o.terminal_window, L);
    TerminalRuleStats local;
    std::vector<OrientedPrimerSearchHit> hits;

    for (const bool reverse : {false, true}) {
        /*
         * Plus strand:  site == primer. The primer's 3' end is the site's
         *               right end -> start at primer[L-1], extend LEFT.
         * Minus strand: site == revcomp(primer). The primer's 3' end is the
         *               site's left end -> start at comp(primer[L-1]),
         *               extend RIGHT through comp(primer[L-2]) ...
         * In both cases depth d consumes primer position L-1-d.
         */
        std::vector<Node> stack;
        for (const char b : kBases) {
            const std::uint8_t want = reverse ? complement_mask(iupac_mask(primer[L - 1]))
                                              : iupac_mask(primer[L - 1]);
            const bool mis = (base_bit(b) & want) == 0;
            const std::size_t term = mis ? 1 : 0;
            if (term > o.max_terminal_mismatches || term > o.max_mismatches) {
                continue;
            }
            auto st = index.search(std::string_view(&b, 1));
            ++local.states_visited;
            if (!st.empty()) {
                stack.push_back(Node{st, 1, term, term});
            }
        }
        while (!stack.empty()) {
            Node n = stack.back();
            stack.pop_back();
            if (n.depth == L) {
                for (const auto pos : index.locate_unsorted(n.state)) {
                    ++local.located;
                    OrientedPrimerSearchHit h;
                    if (verify(reference, primer, pos, reverse, o, h)) {
                        hits.push_back(h);
                    } else {
                        ++local.rejected_on_verify;
                    }
                }
                continue;
            }
            const std::size_t p = L - 1 - n.depth;               // primer position
            const bool in_window = n.depth < window;
            const std::uint8_t want = reverse ? complement_mask(iupac_mask(primer[p]))
                                              : iupac_mask(primer[p]);
            if (!reverse) {
                const auto kids = index.extend_left_all(n.state);
                for (int k = 0; k < 4; ++k) {
                    ++local.states_visited;
                    if (kids[k].empty()) continue;
                    const bool mis = (base_bit(kBases[k]) & want) == 0;
                    const std::size_t mm = n.mismatches + (mis ? 1 : 0);
                    const std::size_t term = n.terminal + (mis && in_window ? 1 : 0);
                    if (mm > o.max_mismatches || term > o.max_terminal_mismatches) continue;
                    stack.push_back(Node{kids[k], n.depth + 1, mm, term});
                }
            } else {
                for (const char b : kBases) {
                    ++local.states_visited;
                    const bool mis = (base_bit(b) & want) == 0;
                    const std::size_t mm = n.mismatches + (mis ? 1 : 0);
                    const std::size_t term = n.terminal + (mis && in_window ? 1 : 0);
                    if (mm > o.max_mismatches || term > o.max_terminal_mismatches) continue;
                    auto st = index.extend_right(n.state, b);
                    if (st.empty()) continue;
                    stack.push_back(Node{st, n.depth + 1, mm, term});
                }
            }
        }
    }

    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
        if (a.position != b.position) return a.position < b.position;
        if (a.orientation != b.orientation)
            return static_cast<int>(a.orientation) < static_cast<int>(b.orientation);
        return a.mismatch_mask < b.mismatch_mask;
    });
    hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
    if (stats) *stats = local;
    return hits;
}

}  // namespace primerpair
