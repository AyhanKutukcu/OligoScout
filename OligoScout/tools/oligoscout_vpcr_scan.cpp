// oligoscout_vpcr_scan — exhaustive oligo site scan of any FASTA (virtual PCR core).
//
// Every position of every record is compared with every oligo on both strands
// (bit-parallel: four one-hot base planes, AND with the oligo's IUPAC planes,
// popcount). No index, no heuristic seed: the result is complete for the
// stated rule. Used by reference_py/oligoscout_thermo/vpcr.py for user
// templates (plasmids, amplicon sets, viral/bacterial assemblies, extra
// contigs) and for degenerate-primer genome scans.
//
// Usage:
//   oligoscout_vpcr_scan FASTA OLIGOS.tsv OUT.tsv [options]
//     --max-mismatch K        total mismatches per site (default 3)
//     --terminal-window W     3'-terminal window (default 7)
//     --terminal-max T        mismatches allowed in that window (default 1)
//     --circular all|NAME[,NAME..]   records treated as circular (sites may
//                             wrap the origin); default none
//     --template-n reject|mismatch|match  template N: site rejected (default, as the
//                             FM-index search), counted as a mismatch, or matching
//     --g-tolerant K          up to K G-containing mismatches (G.T, G.G, G.A between
//                             oligo and template strand) are reported but not
//                             counted against the limits (Kalendar et al. 2024)
//     --threads N             (default 2)
//     --primary-only          map NC_0000NN / chrN names, keep chr1..22,X,Y only
//
//   OLIGOS.tsv: oligo_id <TAB> core sequence (IUPAC, 8..512 nt); header allowed.
//   Template IUPAC codes match when the base sets intersect.
//   OUT.tsv: oligo_id chrom start end_exclusive strand mismatches terminal_mismatches
//            mask site_plus_strand wraps g_mismatches
//     mismatches counts every mismatch; with --g-tolerant the limits apply to
//     mismatches - g_mismatches (mask marks tolerated ones with 'g').
//     start is 0-based on the plus strand; on a circular record end_exclusive
//     may exceed the record length (wraps=1). mask is read from the oligo's
//     5' end ('X' = mismatch), as in oligoscout_panel's bindings.
//   A records file OUT.tsv.records lists name, length and topology.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

int iupac(char c) {
    switch (c) {
        case 'A': case 'a': return 1; case 'C': case 'c': return 2; case 'G': case 'g': return 4;
        case 'T': case 't': case 'U': case 'u': return 8;
        case 'R': case 'r': return 5; case 'Y': case 'y': return 10; case 'S': case 's': return 6;
        case 'W': case 'w': return 9; case 'K': case 'k': return 12; case 'M': case 'm': return 3;
        case 'B': case 'b': return 14; case 'D': case 'd': return 13; case 'H': case 'h': return 11;
        case 'V': case 'v': return 7; case 'N': case 'n': return 15;
        default: return 0;
    }
}
int comp_mask(int m) { return ((m & 1) << 3) | ((m & 8) >> 3) | ((m & 2) << 1) | ((m & 4) >> 1); }

constexpr int kMaxWords = 8;  // 512 nt

struct Oligo {
    std::string id, seq;
    int len = 0, words = 0;
    // planes[strand][base][word]: bit set if the oligo (strand 0) or its
    // reverse complement (strand 1) accepts `base` at that window position.
    // Window position j from the left maps to bit (len-1-j) (newest base = bit 0).
    std::uint64_t planes[2][4][kMaxWords] = {};
    std::uint64_t lenmask[kMaxWords] = {};
    std::uint64_t term[2][kMaxWords] = {};   // 3'-terminal window bits per strand
    // tol[strand][base][word]: mismatch at this position with template `base`
    // forms a G-containing pair (G.T, G.G, G.A)
    std::uint64_t tol[2][4][kMaxWords] = {};
};

struct Opts {
    int maxmm = 3, twin = 7, tmax = 1, threads = 2;
    int template_n = 0;   // 0 reject, 1 mismatch, 2 match
    bool primary_only = false, all_circular = false;
    int gmax = -1;        // -1: G-tolerance off
    std::set<std::string> circular;
};

void set_bit(std::uint64_t* w, int bit) { w[bit / 64] |= std::uint64_t{1} << (bit % 64); }

Oligo make_oligo(const std::string& id, const std::string& s, const Opts& o) {
    Oligo g; g.id = id; g.seq = s; g.len = static_cast<int>(s.size()); g.words = (g.len + 63) / 64;
    for (int b = 0; b < g.len; ++b) set_bit(g.lenmask, b);
    for (int j = 0; j < g.len; ++j) {
        const int m = iupac(s[j]);
        const int rcm = comp_mask(iupac(s[g.len - 1 - j]));   // revcomp at window position j
        const int bit = g.len - 1 - j;
        for (int b = 0; b < 4; ++b) {
            if (m & (1 << b)) set_bit(g.planes[0][b], bit);
            if (rcm & (1 << b)) set_bit(g.planes[1][b], bit);
        }
    }
    // G-pair tolerance. Window base b is the template's plus strand. The oligo
    // pairs with the opposite strand: plus-strand oligo base p (window j) pairs
    // with comp(b); minus-strand oligo (revcomp written in the window as q)
    // pairs with b through its own base comp(q).
    auto stable = [](int x, int y) {      // x, y in 0..3 = A C G T, a mismatch pair x.y
        const bool gx = x == 2, gy = y == 2;
        auto ok = [](int z) { return z == 0 || z == 2 || z == 3; };   // A, G, T
        return (gx && ok(y)) || (gy && ok(x));
    };
    const int cmp[4] = {3, 2, 1, 0};
    for (int j = 0; j < g.len; ++j) {
        const int bit = g.len - 1 - j;
        const int m = iupac(s[j]);
        const int rcm = comp_mask(iupac(s[g.len - 1 - j]));
        for (int b = 0; b < 4; ++b) {
            if (!(m & (1 << b))) {
                bool any = false;    // plus: oligo base x pairs with comp(b)
                for (int x = 0; x < 4; ++x) if ((m & (1 << x)) && stable(x, cmp[b])) any = true;
                if (any) set_bit(g.tol[0][b], bit);
            }
            if (!(rcm & (1 << b))) {
                bool any = false;    // minus: oligo base is comp(q) for window letter q, pairs with b
                for (int q = 0; q < 4; ++q) if ((rcm & (1 << q)) && stable(cmp[q], b)) any = true;
                if (any) set_bit(g.tol[1][b], bit);
            }
        }
    }
    const int w = std::min(o.twin, g.len);
    for (int k = 0; k < w; ++k) {
        set_bit(g.term[0], k);                 // plus: 3' end = window right = low bits
        set_bit(g.term[1], g.len - 1 - k);     // minus: 3' end = window left = high bits
    }
    return g;
}

struct Hit { int o; std::uint64_t start; char strand; int mm, tm, gm; std::uint64_t bits[kMaxWords]; std::uint64_t gbits[kMaxWords]; };

// shift a multi-word register left by one and insert `in` at bit 0
inline void push(std::uint64_t* r, int words, std::uint64_t in) {
    for (int k = words - 1; k > 0; --k) r[k] = (r[k] << 1) | (r[k - 1] >> 63);
    r[0] = (r[0] << 1) | in;
}

void scan(const std::string& seq, std::size_t n_real, std::size_t from, std::size_t to,
          const std::vector<Oligo>& oligos, int words, const Opts& o, std::vector<Hit>& out) {
    // planes of the template window: tp[base][word]
    std::uint64_t tp[4][kMaxWords] = {};
    std::uint64_t np[kMaxWords] = {};       // N positions in the window
    int maxlen = 0;
    for (const auto& g : oligos) maxlen = std::max(maxlen, g.len);
    const std::size_t begin = from >= static_cast<std::size_t>(maxlen) ? from - maxlen : 0;
    std::uint64_t filled = 0;
    for (std::size_t i = begin; i < to; ++i) {
        const int m = iupac(seq[i]);
        const bool is_n = (m == 15 && (seq[i] == 'N' || seq[i] == 'n'));
        const int eff = is_n ? (o.template_n == 2 ? 15 : 0) : m;
        for (int b = 0; b < 4; ++b) push(tp[b], words, (eff >> b) & 1);
        push(np, words, (is_n && o.template_n == 0) ? 1 : 0);
        ++filled;
        if (i < from) continue;
        for (std::size_t k = 0; k < oligos.size(); ++k) {
            const Oligo& g = oligos[k];
            if (filled < static_cast<std::uint64_t>(g.len)) continue;
            const std::uint64_t start = i + 1 - g.len;
            if (start >= n_real) continue;      // circular padding: only starts inside the record
            bool has_n = false;
            for (int w = 0; w < g.words; ++w) if (np[w] & g.lenmask[w]) { has_n = true; break; }
            if (has_n) continue;
            for (int s = 0; s < 2; ++s) {
                int mm = 0, tm = 0, gm = 0, gt = 0;
                std::uint64_t bits[kMaxWords] = {}, gbits[kMaxWords] = {};
                bool over = false;
                for (int w = 0; w < g.words; ++w) {
                    std::uint64_t match = 0, tol = 0;
                    for (int b = 0; b < 4; ++b) {
                        match |= tp[b][w] & g.planes[s][b][w];
                        tol |= tp[b][w] & g.tol[s][b][w];
                    }
                    const std::uint64_t mis = ~match & g.lenmask[w];
                    const std::uint64_t gmis = o.gmax >= 0 ? (mis & tol) : 0;
                    const std::uint64_t counted = mis & ~gmis;
                    bits[w] = mis;
                    gbits[w] = gmis;
                    gm += __builtin_popcountll(gmis);
                    gt += __builtin_popcountll(gmis & g.term[s][w]);
                    mm += __builtin_popcountll(counted);
                    tm += __builtin_popcountll(counted & g.term[s][w]);
                    if (mm > o.maxmm || tm > o.tmax) { over = true; break; }
                }
                if (!over && o.gmax >= 0 && gm > o.gmax) {
                    // more G-pairs than tolerated: tolerate terminal ones first
                    // (most favourable), the excess counts as ordinary mismatches
                    const int tol_t = std::min(gt, o.gmax);
                    const int tol_o = std::min(gm - gt, o.gmax - tol_t);
                    mm += (gt - tol_t) + (gm - gt - tol_o);
                    tm += gt - tol_t;
                    if (mm > o.maxmm || tm > o.tmax) over = true;
                }
                if (over) continue;
                Hit h{static_cast<int>(k), start, s ? '-' : '+', 0, tm, gm, {}, {}};
                h.mm = 0;
                for (int w = 0; w < g.words; ++w) h.mm += __builtin_popcountll(bits[w]);
                std::copy(bits, bits + kMaxWords, h.bits);
                std::copy(gbits, gbits + kMaxWords, h.gbits);
                out.push_back(h);
            }
        }
    }
}

std::string primary_name(const std::string& tok) {
    if (tok.rfind("NC_0000", 0) == 0) {
        const int n = std::atoi(tok.substr(7, 2).c_str());
        if (n >= 1 && n <= 22) return "chr" + std::to_string(n);
        if (n == 23) return "chrX";
        if (n == 24) return "chrY";
        return "";
    }
    const std::string t = tok.rfind("chr", 0) == 0 ? tok : "chr" + tok;
    for (int i = 1; i <= 22; ++i) if (t == "chr" + std::to_string(i)) return t;
    return (t == "chrX" || t == "chrY") ? t : "";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: oligoscout_vpcr_scan FASTA OLIGOS.tsv OUT.tsv [--max-mismatch K] "
                     "[--terminal-window W] [--terminal-max T] [--circular all|N1,N2] "
                     "[--template-n reject|mismatch|match] [--g-tolerant K] [--threads N] [--primary-only]\n";
        return 2;
    }
    Opts o;
    for (int i = 4; i < argc; ++i) {
        const std::string f = argv[i];
        auto val = [&]() -> std::string {
            if (i + 1 >= argc) { std::cerr << f << " needs a value\n"; std::exit(2); }
            return argv[++i];
        };
        if (f == "--max-mismatch") o.maxmm = std::atoi(val().c_str());
        else if (f == "--terminal-window") o.twin = std::atoi(val().c_str());
        else if (f == "--terminal-max") o.tmax = std::atoi(val().c_str());
        else if (f == "--threads") o.threads = std::max(1, std::atoi(val().c_str()));
        else if (f == "--primary-only") o.primary_only = true;
        else if (f == "--g-tolerant") o.gmax = std::max(0, std::atoi(val().c_str()));
        else if (f == "--template-n") {
            const std::string v = val();
            o.template_n = v == "match" ? 2 : v == "mismatch" ? 1 : 0;
        }
        else if (f == "--circular") {
            const std::string v = val();
            if (v == "all") o.all_circular = true;
            else { std::stringstream ss(v); std::string t; while (std::getline(ss, t, ',')) o.circular.insert(t); }
        } else { std::cerr << "unknown option " << f << "\n"; return 2; }
    }
    std::vector<Oligo> oligos;
    int words = 1, maxlen = 0;
    {
        std::ifstream in(argv[2]);
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            std::istringstream ss(line);
            std::string id, s;
            ss >> id >> s;
            if (s.empty()) continue;
            bool ok = true;
            for (char& c : s) { if (!iupac(c)) ok = false; c = static_cast<char>(std::toupper(c)); if (c == 'U') c = 'T'; }
            if (!ok) { if (id != "oligo_id" && id != "id") std::cerr << "skip " << id << " (non-IUPAC)\n"; continue; }
            if (s.size() < 8 || s.size() > 64 * kMaxWords) { std::cerr << "skip " << id << " (length 8..512)\n"; continue; }
            oligos.push_back(make_oligo(id, s, o));
            words = std::max(words, oligos.back().words);
            maxlen = std::max(maxlen, oligos.back().len);
        }
    }
    if (oligos.empty()) { std::cerr << "no oligos\n"; return 2; }
    std::ifstream fa(argv[1]);
    if (!fa) { std::cerr << "cannot open " << argv[1] << "\n"; return 2; }
    std::ofstream out(argv[3]);
    std::ofstream rec(std::string(argv[3]) + ".records");
    out << "oligo_id\tchrom\tstart\tend_exclusive\tstrand\tmismatches\tterminal_mismatches\tmask\tsite_plus_strand\twraps\tg_mismatches\n";
    rec << "chrom\tlength\ttopology\n";
    std::string line, name, seq;
    std::size_t total = 0, scanned = 0;
    auto flush = [&]() {
        if (name.empty()) { seq.clear(); return; }
        const std::size_t n = seq.size();
        const bool circ = o.all_circular || o.circular.count(name);
        rec << name << '\t' << n << '\t' << (circ ? "circular" : "linear") << '\n';
        std::string work = seq;
        if (circ && n > 0) work += seq.substr(0, std::min<std::size_t>(n, static_cast<std::size_t>(maxlen - 1)));
        const std::size_t span = work.size();
        std::vector<std::vector<Hit>> parts(o.threads);
        std::vector<std::thread> th;
        const std::size_t chunk = span / o.threads + 1;
        for (int t = 0; t < o.threads; ++t) {
            const std::size_t a = t * chunk, b = std::min(span, (t + 1) * chunk);
            th.emplace_back([&, t, a, b] { if (a < b) scan(work, n, a, b, oligos, words, o, parts[t]); });
        }
        for (auto& t : th) t.join();
        for (auto& part : parts) for (const Hit& h : part) {
            const Oligo& g = oligos[h.o];
            std::string mask(g.len, '.');
            for (int j = 0; j < g.len; ++j) {
                const int bit = g.len - 1 - j;
                if ((h.bits[bit / 64] >> (bit % 64)) & 1ULL)
                    mask[h.strand == '+' ? j : g.len - 1 - j] = ((h.gbits[bit / 64] >> (bit % 64)) & 1ULL) ? 'g' : 'X';
            }
            std::string site = work.substr(h.start, g.len);
            std::transform(site.begin(), site.end(), site.begin(), ::toupper);
            const bool wraps = h.start + g.len > n;
            out << g.id << '\t' << name << '\t' << h.start << '\t' << h.start + g.len << '\t' << h.strand << '\t'
                << h.mm << '\t' << h.tm << '\t' << mask << '\t' << site << '\t' << (wraps ? 1 : 0) << '\t'
                << h.gm << '\n';
            ++total;
        }
        scanned += n;
        seq.clear();
    };
    while (std::getline(fa, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') {
            flush();
            std::string tok = line.substr(1, line.find_first_of(" \t") - 1);
            name = o.primary_only ? primary_name(tok) : tok;
            continue;
        }
        if (!name.empty()) seq += line;
    }
    flush();
    std::cerr << "done: " << scanned << " bp, " << oligos.size() << " oligos, " << total << " sites\n";
    return 0;
}
