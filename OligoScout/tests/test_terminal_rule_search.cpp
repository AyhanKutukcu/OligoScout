// Terminal-rule (sensitive) search vs brute force on a synthetic reference
// with planted sites: both strands, mismatches inside and outside the 3'
// window, and N blocks.
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/terminal_rule_search.hpp"

#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

using namespace primerpair;

namespace {

std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
    return r;
}

char other(char b, int k) {
    const char* B = "ACGT";
    int i = 0;
    while (B[i] != b) ++i;
    return B[(i + 1 + k % 3) % 4];
}

using Key = std::tuple<std::uint64_t, int, std::uint64_t, std::size_t>;

std::set<Key> brute(const std::string& ref, const std::string& p, const TerminalRuleOptions& o) {
    std::set<Key> out;
    const std::size_t L = p.size();
    const std::string rc = revcomp(p);
    for (std::size_t pos = 0; pos + L <= ref.size(); ++pos) {
        for (int strand = 0; strand < 2; ++strand) {
            std::uint64_t mask = 0; std::size_t mm = 0, term = 0; bool n = false;
            for (std::size_t j = 0; j < L; ++j) {
                const char obs = ref[pos + j];
                if (obs == 'N') { n = true; break; }
                const std::size_t pp = strand ? L - 1 - j : j;
                const char want = strand ? rc[j] : p[j];
                if (obs != want) { mask |= std::uint64_t{1} << pp; ++mm; if (pp + o.terminal_window >= L) ++term; }
            }
            if (!n && mm <= o.max_mismatches && term <= o.max_terminal_mismatches)
                out.insert({pos, strand, mask, mm});
        }
    }
    return out;
}

}  // namespace

int main() {
    try {
        std::uint64_t st = 0x9E3779B97F4A7C15ULL;
        auto rnd = [&]() { st = st * 6364136223846793005ULL + 1442695040888963407ULL; return st >> 33; };
        std::string ref;
        for (int i = 0; i < 60000; ++i) ref.push_back("ACGT"[rnd() & 3]);
        const std::string primer = "CCTTATCAAAGCGGCTCTTG";
        const std::size_t L = primer.size();
        int planted = 0;
        // plant variants: mismatch sets at different 3' distances, both strands
        const int patterns[][3] = {{-1,-1,-1},{19,-1,-1},{15,-1,-1},{10,-1,-1},{19,18,-1},{3,8,17},
                                   {0,5,12},{16,-1,-1},{14,15,-1},{2,9,19}};
        for (int s = 0; s < 2; ++s) {
            for (const auto& pat : patterns) {
                std::string site = primer;
                for (int k = 0; k < 3; ++k) if (pat[k] >= 0) site[pat[k]] = other(site[pat[k]], k + s);
                if (s) site = revcomp(site);
                const std::size_t pos = 500 + planted * 2000;
                ref.replace(pos, L, site);
                ++planted;
            }
        }
        ref.replace(45000, 30, std::string(30, 'N'));
        // a site straddling the N block must be rejected
        std::string straddle = primer; ref.replace(45000 - 10, L, straddle);
        ref.replace(45000, 30, std::string(30, 'N'));

        BidirectionalFMIndex index(ref);
        PackedReference packed(ref);

        for (std::size_t tmax : {0u, 1u, 2u}) {
            TerminalRuleOptions o; o.max_mismatches = 3; o.terminal_window = 5; o.max_terminal_mismatches = tmax;
            const auto hits = search_terminal_rule(index, packed, primer, o);
            std::set<Key> got;
            for (const auto& h : hits)
                got.insert({h.position, h.orientation == PrimerOrientation::Reverse ? 1 : 0, h.mismatch_mask, h.mismatches});
            const auto want = brute(ref, primer, o);
            if (got != want) {
                std::cerr << "terminal_max=" << tmax << " got " << got.size() << " want " << want.size() << '\n';
                return 1;
            }
            std::cout << "terminal_max=" << tmax << ": " << got.size() << " sites == brute force\n";
        }
        std::cout << "terminal_rule_search: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
