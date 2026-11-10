/*
 * oligoscout_thermo_dump
 *
 * Yerel termodinamik cekirdegin (native_thermo) ciktisini yuksek hassasiyetle
 * TSV olarak basar. Tek amaci DIFERANSIYEL TEST: `cli/diff_cpp.py` ayni
 * dizileri Python referans katmanindan gecirip iki uygulamayi karsilastirir.
 *
 *   oligoscout_thermo_dump --sequences girdi.tsv [--na 50] [--mg 1.5]
 *       [--dntp 0.6] [--primer-nm 200] [--ta 60] [--salt-model owczarzy_2008]
 *
 * Girdi TSV basligi: id<TAB>sequence[<TAB>partner]
 */

#include "primerpair/native_thermo.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream ss(line);
    while (std::getline(ss, cur, '\t')) {
        out.push_back(cur);
    }
    return out;
}

double arg_double(int argc, char** argv, const std::string& flag, double fallback) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (flag == argv[i]) {
            return std::stod(argv[i + 1]);
        }
    }
    return fallback;
}

std::string arg_string(int argc, char** argv, const std::string& flag,
                       const std::string& fallback) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (flag == argv[i]) {
            return argv[i + 1];
        }
    }
    return fallback;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string path = arg_string(argc, argv, "--sequences", "");
    if (path.empty()) {
        std::cerr << "kullanim: oligoscout_thermo_dump --sequences girdi.tsv\n";
        return 2;
    }

    primerpair::thermo::Conditions cond;
    cond.monovalent_M = arg_double(argc, argv, "--na", 50.0) / 1e3;
    cond.magnesium_M = arg_double(argc, argv, "--mg", 1.5) / 1e3;
    cond.dntp_M = arg_double(argc, argv, "--dntp", 0.6) / 1e3;
    cond.primer_M = arg_double(argc, argv, "--primer-nm", 200.0) / 1e9;
    cond.anneal_c = arg_double(argc, argv, "--ta", 60.0);
    const std::string salt_model =
        arg_string(argc, argv, "--salt-model", "owczarzy_2008");

    std::ifstream in(path);
    if (!in) {
        std::cerr << "dosya acilamadi: " << path << "\n";
        return 2;
    }

    std::cout << std::setprecision(17);
    std::cout << "id\tsequence\tdH_kcal\tdS_cal\tn_bp\tself_complementary\t"
                 "gc_fraction\ttm_1M_c\ttm_c\tsalt_regime\tconc_model\t"
                 "theta_ta\tdG37_kcal\tterminal_dG37_5\tselfdimer_any_dG37\t"
                 "selfdimer_3p_dG37\thetero_any_dG37\thetero_3p_dG37\t"
                 "hetero_any_helix\n";

    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (header) {
            header = false;
            if (line.rfind("id\t", 0) == 0) {
                continue;
            }
        }
        const auto cols = split_tab(line);
        if (cols.size() < 2) {
            continue;
        }
        const std::string id = cols[0];
        const std::string seq = cols[1];
        const std::string partner = cols.size() > 2 ? cols[2] : std::string{};

        const auto model = primerpair::thermo::build(seq, cond, "", salt_model);
        const auto d_any = primerpair::thermo::self_dimer(seq, "any");
        const auto d_3p = primerpair::thermo::self_dimer(seq, "3prime");

        double het_any = 0.0;
        double het_3p = 0.0;
        int het_len = 0;
        if (!partner.empty()) {
            const auto h1 = primerpair::thermo::best_dimer(seq, partner, "any");
            const auto h2 = primerpair::thermo::best_dimer(seq, partner, "3prime");
            het_any = h1.found ? h1.dG37_kcal : 0.0;
            het_3p = h2.found ? h2.dG37_kcal : 0.0;
            het_len = h1.helix_len;
        }

        std::cout << id << '\t' << seq << '\t' << model.duplex.dH_kcal << '\t'
                  << model.duplex.dS_cal << '\t' << model.duplex.n_bp << '\t'
                  << (model.duplex.self_complementary ? 1 : 0) << '\t'
                  << model.f_gc << '\t'
                  << model.tm_1m_K - primerpair::thermo::kKelvin << '\t'
                  << model.tm_c() << '\t' << model.salt_regime << '\t'
                  << model.conc_model << '\t' << model.theta(cond.anneal_c)
                  << '\t' << model.dG_kcal(37.0) << '\t'
                  << primerpair::thermo::terminal_dG37(seq, 5) << '\t'
                  << (d_any.found ? d_any.dG37_kcal : 0.0) << '\t'
                  << (d_3p.found ? d_3p.dG37_kcal : 0.0) << '\t' << het_any
                  << '\t' << het_3p << '\t' << het_len << '\n';
    }
    return 0;
}
