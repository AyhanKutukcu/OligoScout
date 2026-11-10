#pragma once

/*
 * NATIVE_THERMO_V1
 *
 * Primer3'ten BAGIMSIZ, yerel nearest-neighbor termodinamigi.
 *
 * Bu baslik `thermodynamic_backend.hpp`in (Primer3 ntthal/oligotm adaptoru)
 * yerine gecmek uzere yazilmistir: ayni buyuklukleri alt surec cagirmadan,
 * tablolari depo icinde tutarak hesaplar.
 *
 * Kapsam
 *   - Watson-Crick NN duplex dH/dS (SantaLucia 1998 birlesik set)
 *   - Tuz duzeltmeleri: Owczarzy 2004 / Owczarzy 2008 / von Ahsen 2001
 *   - Konsantrasyon modelleri: equal | self_complementary | primer_excess
 *   - Tm, theta(T), dG(T)
 *   - 3' uc kararliligi (son N bazin yigin dG37 toplami)
 *   - Dimer: iki oligo arasindaki en kararli yapi. Puanlanan sinif,
 *     Python tarafiyla AYNIDIR: kesintisiz Watson-Crick helisi ya da en
 *     fazla bir kusur (tek ic uyumsuzluk VEYA tek bazlik ic bulge; her
 *     iki yaninda >= 2 WC cifti). Uyumsuzluk ve bulge tablolari
 *     `include/primerpair/generated/thermo_tables.inc` dosyasindan gelir
 *     ve `tools/gen_thermo_tables.py` tarafindan dogrudan params.py'den
 *     uretilir -- C++ tarafinda ikinci bir transkripsiyon yoktur.
 *
 * Bu kapsam, `reference_py/oligoscout_thermo` Python referans katmaniyla
 * BIREBIR ayni modeldir ve `cli/diff_cpp.py` ikisini sayisal olarak
 * karsilastirir (diferansiyel test).
 *
 * Kapsam DISI (Python tarafinda da bekliyor):
 *   - uyumsuzluk ve sarkan uc NN tablolari (kapali erisim, `pending`)
 *   - dongu (hairpin/bulge/ic dongu) parametreleri
 */

#include <string>
#include <string_view>
#include <utility>

namespace primerpair::thermo {

inline constexpr double kR_cal = 1.9872;     // cal / (mol K)
inline constexpr double kKelvin = 273.15;

struct Conditions {
    double monovalent_M{0.050};
    double magnesium_M{0.0015};
    double dntp_M{0.0006};
    double primer_M{2.0e-7};
    double anneal_c{60.0};

    [[nodiscard]] double free_magnesium_M() const noexcept {
        const double v = magnesium_M - dntp_M;
        return v > 0.0 ? v : 0.0;
    }
};

struct DuplexThermo {
    double dH_kcal{0.0};
    double dS_cal{0.0};
    bool self_complementary{false};
    int n_bp{0};

    [[nodiscard]] double dG_kcal(double temp_c = 37.0) const noexcept {
        return dH_kcal - (temp_c + kKelvin) * dS_cal / 1000.0;
    }
};

struct SaltResult {
    double tm_K{0.0};
    std::string regime;
};

struct MeltingModel {
    DuplexThermo duplex{};
    Conditions cond{};
    double f_gc{0.0};
    std::string conc_model;
    std::string salt_model;
    std::string salt_regime;
    double c_eff_M{0.0};
    double tm_1m_K{0.0};
    double tm_K{0.0};
    double dS_eff_cal{0.0};

    [[nodiscard]] double tm_c() const noexcept { return tm_K - kKelvin; }
    [[nodiscard]] double dG_kcal(double temp_c = 37.0) const noexcept;
    [[nodiscard]] double theta(double temp_c) const;
};

struct DimerResult {
    bool found{false};
    double dH_kcal{0.0};
    double dS_cal{0.0};
    double dG37_kcal{0.0};
    int offset{0};
    int helix_len{0};
    bool three_prime_paired{false};
    int mismatch_count{0};   // 0 veya 1 -- tek ic uyumsuzluk
    int bulge_count{0};      // 0 veya 1 -- tek bazlik ic bulge
};

// -- dizi yardimcilari -------------------------------------------------------
[[nodiscard]] std::string reverse_complement(std::string_view seq);
[[nodiscard]] std::string validate(std::string_view seq);     // ACGT disi -> throw
[[nodiscard]] double gc_fraction(std::string_view seq);
[[nodiscard]] bool is_self_complementary(std::string_view seq);

// -- NN ve tuz ---------------------------------------------------------------
[[nodiscard]] DuplexThermo watson_crick(std::string_view seq);
[[nodiscard]] double owczarzy_2004_monovalent(double tm_1m_K, double na_M,
                                              double f_gc);
[[nodiscard]] double von_ahsen_2001_na_eq(double mon_M, double mg_M,
                                          double dntp_M);
[[nodiscard]] SaltResult owczarzy_2008_divalent(double tm_1m_K, double mg_M,
                                                double mon_M, double f_gc,
                                                int n_bp);
[[nodiscard]] SaltResult correct_salt(double tm_1m_K, double f_gc, int n_bp,
                                      const Conditions& cond,
                                      std::string_view model);

// -- erime modeli ------------------------------------------------------------
[[nodiscard]] double effective_concentration(std::string_view conc_model,
                                             const Conditions& cond,
                                             double total_strand_M);
[[nodiscard]] MeltingModel build(std::string_view seq,
                                 const Conditions& cond = {},
                                 std::string_view conc_model = "",
                                 std::string_view salt_model = "owczarzy_2008",
                                 double total_strand_M = -1.0);

// -- QC yardimcilari ---------------------------------------------------------
[[nodiscard]] double terminal_dG37(std::string_view seq, int window = 5);
[[nodiscard]] DimerResult best_dimer(std::string_view a, std::string_view b,
                                     std::string_view mode = "any");
[[nodiscard]] inline DimerResult self_dimer(std::string_view a,
                                            std::string_view mode = "any") {
    return best_dimer(a, a, mode);
}

}  // namespace primerpair::thermo
