#include "primerpair/native_thermo.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <limits>
#include <map>
#include <stdexcept>

namespace primerpair::thermo {
namespace {

/*
 * SantaLucia 1998 birlesik NN seti (1 M NaCl, pH 7).
 * Degerler reference_py/oligoscout_thermo/params.py ile AYNI kaynaktan ve
 * ayni sayilardir; ikisi `cli/diff_cpp.py` ile sayisal olarak karsilastirilir.
 * Durum: transcribed_unverified (birincil PDF kapali erisim).
 */
const std::map<std::string, std::pair<double, double>>& wc_table() {
    static const std::map<std::string, std::pair<double, double>> t = {
        {"AA", {-7.9, -22.2}}, {"AT", {-7.2, -20.4}}, {"TA", {-7.2, -21.3}},
        {"CA", {-8.5, -22.7}}, {"GT", {-8.4, -22.4}}, {"CT", {-7.8, -21.0}},
        {"GA", {-8.2, -22.2}}, {"CG", {-10.6, -27.2}}, {"GC", {-9.8, -24.4}},
        {"GG", {-8.0, -19.9}},
    };
    return t;
}

constexpr std::pair<double, double> kInitAT{2.3, 4.1};
constexpr std::pair<double, double> kInitGC{0.1, -2.8};
constexpr std::pair<double, double> kSymmetry{0.0, -1.4};

// Owczarzy 2008, denklem 16 katsayilari
constexpr double kA = 3.92e-5;
constexpr double kB = -9.11e-6;
constexpr double kC = 6.26e-5;
constexpr double kD = 1.42e-5;
constexpr double kE = -4.82e-4;
constexpr double kF = 5.25e-4;
constexpr double kG = 8.31e-5;

std::pair<double, double> stack(const std::string& step) {
    const auto& t = wc_table();
    auto it = t.find(step);
    if (it != t.end()) {
        return it->second;
    }
    // Ters tamamlayici esdegeri: 5'XY3'/3'X'Y'5' == 5'Y'X'3'/3'YX5'
    it = t.find(reverse_complement(step));
    if (it == t.end()) {
        throw std::runtime_error("NN basamagi tabloda yok: " + step);
    }
    return it->second;
}

bool pairs(char a, char b) {
    return (a == 'A' && b == 'T') || (a == 'T' && b == 'A') ||
           (a == 'G' && b == 'C') || (a == 'C' && b == 'G');
}

std::pair<double, double> stack_energy(const std::string& top) {
    double dH = 0.0;
    double dS = 0.0;
    for (std::size_t i = 0; i + 1 < top.size(); ++i) {
        const auto hs = stack(top.substr(i, 2));
        dH += hs.first;
        dS += hs.second;
    }
    return {dH, dS};
}

std::pair<double, double> init_energy(char base) {
    return (base == 'A' || base == 'T') ? kInitAT : kInitGC;
}

std::pair<double, double> helix_energy(const std::string& top) {
    auto hs = stack_energy(top);
    for (char base : {top.front(), top.back()}) {
        const auto init = init_energy(base);
        hs.first += init.first;
        hs.second += init.second;
    }
    return hs;
}

#include "primerpair/generated/thermo_tables.inc"

/*
 * Tek ic uyumsuzlugun iki komsu yigini icin Watson-Crick'e gore fark.
 * reference_py/oligoscout_thermo/nn.py::_mismatch_stack_adjustment ile ayni
 * sozlesme: anahtar, yigin ile 180 derece donmus esdeginin leksikografik
 * kucugudur.
 */
std::pair<double, double> mismatch_adjustment(const std::string& top,
                                              const std::string& bottom,
                                              int position) {
    double dH = 0.0;
    double dS = 0.0;
    for (int i : {position - 1, position}) {
        const std::string ts = top.substr(static_cast<std::size_t>(i), 2);
        const std::string bs = bottom.substr(static_cast<std::size_t>(i), 2);
        const std::string direct = ts + "/" + bs;
        const std::string rotated =
            std::string{bs[1], bs[0]} + "/" + std::string{ts[1], ts[0]};
        const std::string key = std::min(direct, rotated);
        const auto& table = mismatch_table();
        const auto it = table.find(key);
        if (it == table.end()) {
            throw std::runtime_error("uyumsuzluk yigini tabloda yok: " + key);
        }
        const auto wc = stack(ts);
        dH += it->second.first - wc.first;
        dS += it->second.second - wc.second;
    }
    return {dH, dS};
}

// Tek bazlik bulge uzunluk terimi (dH = 0 kaynak varsayimi).
std::pair<double, double> bulge1_energy() {
    return {0.0, -1000.0 * kBulge1_dG37 / 310.15};
}

}  // namespace

std::string reverse_complement(std::string_view seq) {
    std::string out;
    out.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (*it) {
            case 'A': out.push_back('T'); break;
            case 'T': out.push_back('A'); break;
            case 'C': out.push_back('G'); break;
            case 'G': out.push_back('C'); break;
            default: throw std::invalid_argument("tanimsiz karakter");
        }
    }
    return out;
}

std::string validate(std::string_view seq) {
    std::string s;
    s.reserve(seq.size());
    for (char ch : seq) {
        const char up = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        if (up == ' ' || up == '\t' || up == '\n' || up == '\r') {
            continue;
        }
        if (up != 'A' && up != 'C' && up != 'G' && up != 'T') {
            throw std::invalid_argument("tanimsiz karakter (yalniz A/C/G/T)");
        }
        s.push_back(up);
    }
    if (s.empty()) {
        throw std::invalid_argument("bos dizi");
    }
    return s;
}

double gc_fraction(std::string_view seq) {
    const std::string s = validate(seq);
    const auto gc = std::count_if(s.begin(), s.end(),
                                  [](char c) { return c == 'G' || c == 'C'; });
    return static_cast<double>(gc) / static_cast<double>(s.size());
}

bool is_self_complementary(std::string_view seq) {
    const std::string s = validate(seq);
    return s == reverse_complement(s);
}

DuplexThermo watson_crick(std::string_view seq) {
    const std::string s = validate(seq);
    if (s.size() < 2) {
        throw std::invalid_argument("NN modeli icin en az 2 nt gerekir");
    }
    double dH = 0.0;
    double dS = 0.0;
    for (std::size_t i = 0; i + 1 < s.size(); ++i) {
        const auto hs = stack(s.substr(i, 2));
        dH += hs.first;
        dS += hs.second;
    }
    for (char base : {s.front(), s.back()}) {
        const auto init = (base == 'A' || base == 'T') ? kInitAT : kInitGC;
        dH += init.first;
        dS += init.second;
    }
    const bool selfc = (s == reverse_complement(s));
    if (selfc) {
        dH += kSymmetry.first;
        dS += kSymmetry.second;
    }
    return DuplexThermo{dH, dS, selfc, static_cast<int>(s.size())};
}

double owczarzy_2004_monovalent(double tm_1m_K, double na_M, double f_gc) {
    if (na_M <= 0.0) {
        throw std::invalid_argument("[Na+] pozitif olmali");
    }
    const double ln_na = std::log(na_M);
    const double inv = 1.0 / tm_1m_K + (4.29 * f_gc - 3.95) * 1e-5 * ln_na +
                       9.40e-6 * ln_na * ln_na;
    return 1.0 / inv;
}

double von_ahsen_2001_na_eq(double mon_M, double mg_M, double dntp_M) {
    const double free_mg_mM = std::max(0.0, (mg_M - dntp_M) * 1e3);
    return mon_M + 120.0 * std::sqrt(free_mg_mM) / 1e3;
}

SaltResult owczarzy_2008_divalent(double tm_1m_K, double mg_M, double mon_M,
                                  double f_gc, int n_bp) {
    if (mg_M <= 0.0) {
        if (mon_M <= 0.0) {
            throw std::invalid_argument("ne Mg2+ ne tek degerlikli katyon var");
        }
        return {owczarzy_2004_monovalent(tm_1m_K, mon_M, f_gc), "monovalent_only"};
    }
    const double ratio =
        mon_M > 0.0 ? std::sqrt(mg_M) / mon_M : std::numeric_limits<double>::infinity();
    if (ratio < 0.22) {
        return {owczarzy_2004_monovalent(tm_1m_K, mon_M, f_gc),
                "monovalent_dominant"};
    }

    double a = kA;
    double d = kD;
    double g = kG;
    std::string regime = "divalent_dominant";
    if (ratio < 6.0) {
        const double ln_mon = std::log(mon_M);
        a = kA * (0.843 - 0.352 * std::sqrt(mon_M) * ln_mon);
        d = kD * (1.279 - 4.03e-3 * ln_mon - 8.03e-3 * ln_mon * ln_mon);
        g = kG * (0.486 - 0.258 * ln_mon + 5.25e-3 * ln_mon * ln_mon * ln_mon);
        regime = "mixed_mono_divalent";
    }

    const double ln_mg = std::log(mg_M);
    const double inv = 1.0 / tm_1m_K + a + kB * ln_mg +
                       f_gc * (kC + d * ln_mg) +
                       (1.0 / (2.0 * (n_bp - 1))) *
                           (kE + kF * ln_mg + g * ln_mg * ln_mg);
    return {1.0 / inv, regime};
}

SaltResult correct_salt(double tm_1m_K, double f_gc, int n_bp,
                        const Conditions& cond, std::string_view model) {
    if (model == "owczarzy_2004") {
        return {owczarzy_2004_monovalent(tm_1m_K, cond.monovalent_M, f_gc),
                "monovalent"};
    }
    if (model == "owczarzy_2008") {
        return owczarzy_2008_divalent(tm_1m_K, cond.free_magnesium_M(),
                                      cond.monovalent_M, f_gc, n_bp);
    }
    if (model == "von_ahsen_2001") {
        const double na_eq =
            von_ahsen_2001_na_eq(cond.monovalent_M, cond.magnesium_M, cond.dntp_M);
        return {owczarzy_2004_monovalent(tm_1m_K, na_eq, f_gc), "na_equivalent"};
    }
    throw std::invalid_argument("bilinmeyen tuz modeli");
}

double effective_concentration(std::string_view conc_model,
                               const Conditions& cond, double total_strand_M) {
    if (conc_model == "equal") {
        return total_strand_M / 4.0;
    }
    if (conc_model == "self_complementary") {
        return total_strand_M;
    }
    if (conc_model == "primer_excess") {
        return cond.primer_M;
    }
    throw std::invalid_argument("bilinmeyen konsantrasyon modeli");
}

double MeltingModel::dG_kcal(double temp_c) const noexcept {
    const double t = temp_c + kKelvin;
    return (duplex.dH_kcal * 1000.0 - t * dS_eff_cal) / 1000.0;
}

double MeltingModel::theta(double temp_c) const {
    const double t = temp_c + kKelvin;
    const double dG_cal = duplex.dH_kcal * 1000.0 - t * dS_eff_cal;
    const double k = std::exp(-dG_cal / (kR_cal * t));
    if (conc_model == "primer_excess") {
        const double kp = k * cond.primer_M;
        return kp / (1.0 + kp);
    }
    if (conc_model == "equal") {
        const double a = c_eff_M * 2.0;
        const double ka = k * a;
        const double disc =
            std::sqrt(std::max(0.0, (2 * ka + 1) * (2 * ka + 1) - 4 * ka * ka));
        const double x = 2 * ka * a / ((2 * ka + 1) + disc);
        return x / a;
    }
    if (conc_model == "self_complementary") {
        const double s = c_eff_M;
        const double ks = k * s;
        const double disc =
            std::sqrt(std::max(0.0, (4 * ks + 1) * (4 * ks + 1) - 16 * ks * ks));
        const double x = 2 * ks * s / ((4 * ks + 1) + disc);
        return 2 * x / s;
    }
    throw std::invalid_argument("bilinmeyen konsantrasyon modeli");
}

MeltingModel build(std::string_view seq, const Conditions& cond,
                   std::string_view conc_model, std::string_view salt_model,
                   double total_strand_M) {
    const std::string s = validate(seq);
    const DuplexThermo duplex = watson_crick(s);
    const double f_gc = gc_fraction(s);
    std::string model{conc_model};
    if (model.empty()) {
        model = duplex.self_complementary ? "self_complementary" : "primer_excess";
    }
    if (total_strand_M < 0.0) {
        total_strand_M = cond.primer_M;
    }
    const double c_eff = effective_concentration(model, cond, total_strand_M);
    if (c_eff <= 0.0) {
        throw std::invalid_argument("etkin konsantrasyon pozitif olmali");
    }
    const double denom = duplex.dS_cal + kR_cal * std::log(c_eff);
    if (denom >= 0.0) {
        throw std::invalid_argument("gecersiz termodinamik: dS + R ln C >= 0");
    }
    const double tm_1m_K = duplex.dH_kcal * 1000.0 / denom;
    const SaltResult salt = correct_salt(tm_1m_K, f_gc, duplex.n_bp, cond, salt_model);
    const double dS_eff =
        duplex.dH_kcal * 1000.0 / salt.tm_K - kR_cal * std::log(c_eff);

    MeltingModel m;
    m.duplex = duplex;
    m.cond = cond;
    m.f_gc = f_gc;
    m.conc_model = model;
    m.salt_model = std::string{salt_model};
    m.salt_regime = salt.regime;
    m.c_eff_M = c_eff;
    m.tm_1m_K = tm_1m_K;
    m.tm_K = salt.tm_K;
    m.dS_eff_cal = dS_eff;
    return m;
}

double terminal_dG37(std::string_view seq, int window) {
    const std::string s = validate(seq);
    const std::size_t w = static_cast<std::size_t>(std::max(1, window));
    const std::string tail = s.size() <= w ? s : s.substr(s.size() - w);
    double dH = 0.0;
    double dS = 0.0;
    for (std::size_t i = 0; i + 1 < tail.size(); ++i) {
        const auto hs = stack(tail.substr(i, 2));
        dH += hs.first;
        dS += hs.second;
    }
    // Baslatma terimi BILEREK disarida: literaturdeki "3' end stability"
    // tanimi yalniz yigin toplamidir (Python tarafiyla ayni).
    return dH - 310.15 * dS / 1000.0;
}

DimerResult best_dimer(std::string_view a, std::string_view b,
                       std::string_view mode) {
    if (mode != "any" && mode != "3prime") {
        throw std::invalid_argument("bilinmeyen mod");
    }
    const std::string sa = validate(a);
    const std::string sb = validate(b);
    std::string rb{sb.rbegin(), sb.rend()};          // b, 3'->5'
    const int na = static_cast<int>(sa.size());
    const int nb = static_cast<int>(rb.size());

    DimerResult best{};

    auto consider = [&](double dH, double dS, int offset, int helix_len,
                        bool terminal, int mismatches, int bulges) {
        const double dG = dH - 310.15 * dS / 1000.0;
        if (best.found && dG >= best.dG37_kcal) {
            return;
        }
        best.found = true;
        best.dH_kcal = dH;
        best.dS_cal = dS;
        best.dG37_kcal = dG;
        best.offset = offset;
        best.helix_len = helix_len;
        best.three_prime_paired = terminal;
        best.mismatch_count = mismatches;
        best.bulge_count = bulges;
    };

    // --- kesintisiz helis ve tek ic uyumsuzluk -----------------------------
    for (int offset = -(nb - 1); offset < na; ++offset) {
        const int lo = std::max(0, offset);
        const int hi = std::min(na, offset + nb);
        if (hi - lo < 2) {
            continue;
        }
        int start = -1;
        int prev_lo = -1;
        int prev_hi = -1;
        for (int idx = lo; idx <= hi; ++idx) {
            const bool paired =
                idx < hi && pairs(sa[static_cast<std::size_t>(idx)],
                                  rb[static_cast<std::size_t>(idx - offset)]);
            if (paired) {
                if (start < 0) {
                    start = idx;
                }
                continue;
            }
            if (start < 0) {
                continue;
            }
            const int run_lo = start;
            const int run_hi = idx;
            start = -1;
            if (run_hi - run_lo < 2) {
                prev_lo = run_lo;
                prev_hi = run_hi;
                continue;
            }
            // Aday 1: kosunun kendisi. Aday 2: tek eslesmemis sutunla ayrilan
            // onceki kosuyla birlesimi (tek ic uyumsuzluk).
            for (int variant = 0; variant < 2; ++variant) {
                int left = run_lo;
                int mismatch = -1;
                if (variant == 1) {
                    if (prev_hi < 0 || run_lo != prev_hi + 1 ||
                        prev_hi - prev_lo < 2) {
                        continue;
                    }
                    left = prev_lo;
                    mismatch = prev_hi;
                }
                const bool terminal = (run_hi == na) || (left == offset);
                if (mode == "3prime" && !terminal) {
                    continue;
                }
                const std::string top =
                    sa.substr(static_cast<std::size_t>(left),
                              static_cast<std::size_t>(run_hi - left));
                auto hs = helix_energy(top);
                if (mismatch >= 0) {
                    const std::string bottom =
                        rb.substr(static_cast<std::size_t>(left - offset),
                                  static_cast<std::size_t>(run_hi - left));
                    const auto adj =
                        mismatch_adjustment(top, bottom, mismatch - left);
                    hs.first += adj.first;
                    hs.second += adj.second;
                }
                const int count = (mismatch >= 0) ? 1 : 0;
                consider(hs.first, hs.second, offset, run_hi - left - count,
                         terminal, count, 0);
            }
            prev_lo = run_lo;
            prev_hi = run_hi;
        }
    }

    // --- tek bazlik ic bulge -----------------------------------------------
    // Sol segment `offset`, sag segment `offset -+ 1` kaymasindadir; aradaki
    // tek baz eslesmeden kalir. Iki segment de bulge'dan UZAGA buyutulur;
    // bu sozlesme aynali simetriyi garanti eder (align.py ile ayni).
    const auto bulge = bulge1_energy();
    for (int offset = -(nb - 1); offset < na; ++offset) {
        const int lo = std::max(0, offset);
        const int hi = std::min(na, offset + nb);
        for (int side = 0; side < 2; ++side) {
            const bool on_b = (side == 0);
            const int o2 = on_b ? offset - 1 : offset + 1;
            const int lo2 = std::max(0, o2);
            const int hi2 = std::min(na, o2 + nb);
            for (int h1 = lo + 2; h1 <= hi; ++h1) {
                if (!pairs(sa[static_cast<std::size_t>(h1 - 1)],
                           rb[static_cast<std::size_t>(h1 - 1 - offset)])) {
                    continue;
                }
                int left = h1 - 1;
                while (left - 1 >= lo &&
                       pairs(sa[static_cast<std::size_t>(left - 1)],
                             rb[static_cast<std::size_t>(left - 1 - offset)])) {
                    --left;
                }
                if (h1 - left < 2) {
                    continue;
                }
                const int l2 = on_b ? h1 : h1 + 1;
                if (l2 < lo2 || l2 >= hi2) {
                    continue;
                }
                if (!pairs(sa[static_cast<std::size_t>(l2)],
                           rb[static_cast<std::size_t>(l2 - o2)])) {
                    continue;
                }
                int h2 = l2;
                while (h2 < hi2 && pairs(sa[static_cast<std::size_t>(h2)],
                                         rb[static_cast<std::size_t>(h2 - o2)])) {
                    ++h2;
                }
                if (h2 - l2 < 2) {
                    continue;
                }
                const bool terminal = (h2 == na) || (left == offset);
                if (mode == "3prime" && !terminal) {
                    continue;
                }

                double dH = 0.0;
                double dS = 0.0;
                for (const auto& segment :
                     {sa.substr(static_cast<std::size_t>(left),
                                static_cast<std::size_t>(h1 - left)),
                      sa.substr(static_cast<std::size_t>(l2),
                                static_cast<std::size_t>(h2 - l2))}) {
                    const auto hs = stack_energy(segment);
                    dH += hs.first;
                    dS += hs.second;
                }
                for (char base : {sa[static_cast<std::size_t>(left)],
                                  sa[static_cast<std::size_t>(h2 - 1)]}) {
                    const auto init = init_energy(base);
                    dH += init.first;
                    dS += init.second;
                }
                dH += bulge.first;
                dS += bulge.second;
                // Bulge'u kapatan her A/T cifti icin ek ceza
                // (align.py::_AT_BULGE_CLOSURE_DG37 ile ayni sozlesme).
                {
                    const char c1 = sa[static_cast<std::size_t>(h1 - 1)];
                    const char c2 = sa[static_cast<std::size_t>(l2)];
                    const int at_closures = ((c1 == 'A' || c1 == 'T') ? 1 : 0) +
                                            ((c2 == 'A' || c2 == 'T') ? 1 : 0);
                    dS -= 1000.0 * kAtBulgeClosure_dG37 * at_closures / 310.15;
                }
                // Denklem 13: bulge uzerinden gecen komsu WC yigini,
                // kesintisiz iplikte okunur.
                std::string step;
                if (on_b) {
                    step = {sa[static_cast<std::size_t>(h1 - 1)],
                            sa[static_cast<std::size_t>(h1)]};
                } else {
                    const int j = h1 - 1 - offset;
                    step = {rb[static_cast<std::size_t>(j + 1)],
                            rb[static_cast<std::size_t>(j)]};
                }
                const auto across = stack(step);
                dH += across.first;
                dS += across.second;

                consider(dH, dS, offset, (h1 - left) + (h2 - l2), terminal,
                         0, 1);
            }
        }
    }
    return best;
}

}  // namespace primerpair::thermo
