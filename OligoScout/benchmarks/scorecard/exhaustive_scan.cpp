// exhaustive_scan.cpp — brute-force, UNANCHORED <=K-mismatch primer scan of a FASTA genome.
//
// Independent of OligoScout's FM-index: every reference position on both strands
// is compared with every primer (2-bit packed XOR + popcount). Used to measure the
// "anchor blind spot" (sites whose mismatches fall inside the 3' 12-nt anchor),
// which an anchor-seeded search cannot enumerate by design.
//
// Build:  g++ -std=c++17 -O3 -march=native -pthread exhaustive_scan.cpp -o exhaustive_scan
// Usage:  exhaustive_scan GENOME.fa PRIMERS.tsv OUT.tsv [max_mm=3] [threads=2]
//   PRIMERS.tsv: primer_id <TAB> sequence (A/C/G/T, 8..32 nt), header line allowed
//   OUT.tsv: primer_id chrom start end_exclusive strand mismatches mask site_plus_strand
//   mask: position i from the primer's 5' end, 'X' = mismatch (same convention as
//   oligoscout_panel's panel.bindings.tsv).
// Only the 24 primary chromosomes are scanned (chr1..22, X, Y; NC_0000NN accessions
// are mapped to chrN).
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

struct Pat { std::string id, seq; int len; uint64_t fwd, rc, lenmask; };

static int code(char c) {
  switch (c) { case 'A': case 'a': return 0; case 'C': case 'c': return 1;
               case 'G': case 'g': return 2; case 'T': case 't': return 3; default: return -1; }
}
static std::string revcomp(const std::string& s) {
  std::string r(s.rbegin(), s.rend());
  for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
  return r;
}
static uint64_t pack(const std::string& s) { uint64_t w = 0; for (char c : s) w = (w << 2) | code(c); return w; }

static std::string chrom_name(const std::string& header) {
  std::string tok = header.substr(1, header.find_first_of(" \t") - 1);
  if (tok.rfind("NC_0000", 0) == 0) {
    int n = std::stoi(tok.substr(7, 2));
    if (n >= 1 && n <= 22) return "chr" + std::to_string(n);
    if (n == 23) return "chrX";
    if (n == 24) return "chrY";
    return "";
  }
  static const char* ok[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
    "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"};
  std::string t = tok.rfind("chr", 0) == 0 ? tok : "chr" + tok;
  for (auto* o : ok) if (t == o) return t;
  return "";
}

struct Hit { int p; uint32_t start; char strand; int mm; uint64_t bits; };

static void scan_range(const std::string& seq, size_t from, size_t to, const std::vector<Pat>& pats, int maxmm,
                       std::vector<Hit>& out) {
  const uint64_t M1 = 0x5555555555555555ULL;
  uint64_t w = 0; int valid = 0;  // number of consecutive valid bases ending here
  size_t begin = from >= 32 ? from - 32 : 0;
  for (size_t i = begin; i < to; ++i) {
    int c = code(seq[i]);
    if (c < 0) { valid = 0; w = 0; continue; }
    w = (w << 2) | (uint64_t)c; if (valid < 64) ++valid;
    if (i < from) continue;
    for (size_t k = 0; k < pats.size(); ++k) {
      const Pat& p = pats[k];
      if (valid < p.len) continue;
      uint64_t win = w & p.lenmask;
      uint64_t x = win ^ p.fwd; uint64_t m = (x | (x >> 1)) & M1 & p.lenmask;
      int mm = __builtin_popcountll(m);
      if (mm <= maxmm) out.push_back({(int)k, (uint32_t)(i + 1 - p.len), '+', mm, m});
      x = win ^ p.rc; m = (x | (x >> 1)) & M1 & p.lenmask;
      mm = __builtin_popcountll(m);
      if (mm <= maxmm) out.push_back({(int)k, (uint32_t)(i + 1 - p.len), '-', mm, m});
    }
  }
}

int main(int argc, char** argv) {
  if (argc < 4) { std::cerr << "usage: exhaustive_scan GENOME.fa PRIMERS.tsv OUT.tsv [max_mm] [threads]\n"; return 2; }
  int maxmm = argc > 4 ? std::atoi(argv[4]) : 3;
  int nthreads = argc > 5 ? std::atoi(argv[5]) : 2;
  std::vector<Pat> pats;
  { std::ifstream in(argv[2]); std::string line;
    while (std::getline(in, line)) {
      std::istringstream ss(line); std::string id, s; ss >> id >> s;
      if (s.empty() || s.find_first_not_of("ACGT") != std::string::npos) continue;
      if (s.size() < 8 || s.size() > 32) { std::cerr << "skip " << id << " (length)\n"; continue; }
      int L = (int)s.size(); uint64_t lm = L == 32 ? ~0ULL : ((1ULL << (2 * L)) - 1);
      pats.push_back({id, s, L, pack(s), pack(revcomp(s)), lm});
    } }
  std::cerr << pats.size() << " primers, max_mm " << maxmm << ", threads " << nthreads << "\n";
  std::ifstream fa(argv[1]); std::ofstream out(argv[3]);
  out << "primer_id\tchrom\tstart\tend_exclusive\tstrand\tmismatches\tmask\tsite_plus_strand\n";
  std::string line, name, seq; seq.reserve(260000000);
  size_t total = 0, scanned = 0;
  auto flush = [&]() {
    if (name.empty()) { seq.clear(); return; }
    std::vector<std::vector<Hit>> parts(nthreads); std::vector<std::thread> th;
    size_t chunk = seq.size() / nthreads + 1;
    for (int t = 0; t < nthreads; ++t) {
      size_t a = t * chunk, b = std::min(seq.size(), (t + 1) * chunk);
      th.emplace_back([&, t, a, b] { if (a < b) scan_range(seq, a, b, pats, maxmm, parts[t]); });
    }
    for (auto& t : th) t.join();
    size_t n = 0;
    for (auto& part : parts) for (const Hit& h : part) {
      const Pat& p = pats[h.p]; std::string mask(p.len, '.');
      for (int j = 0; j < p.len; ++j) {           // j = window position from the left (plus strand)
        bool bad = (h.bits >> (2 * (p.len - 1 - j))) & 1ULL;
        if (bad) mask[h.strand == '+' ? j : p.len - 1 - j] = 'X';
      }
      std::string site = seq.substr(h.start, p.len);
      std::transform(site.begin(), site.end(), site.begin(), ::toupper);
      out << p.id << '\t' << name << '\t' << h.start << '\t' << h.start + p.len << '\t' << h.strand << '\t'
          << h.mm << '\t' << mask << '\t' << site << '\n';
      ++n;
    }
    total += n; scanned += seq.size();
    std::cerr << name << " " << seq.size() << " bp, hits " << n << "\n";
    seq.clear();
  };
  while (std::getline(fa, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty() && line[0] == '>') { flush(); name = chrom_name(line); continue; }
    if (!name.empty()) seq += line;
  }
  flush();
  std::cerr << "done: " << scanned << " bp scanned, " << total << " hits\n";
  return 0;
}
