// OTOMATIK URETILDI -- ELLE DUZENLEMEYIN.
// Ureten : tools/gen_thermo_tables.py
// Kaynak : reference_py/oligoscout_thermo/params.py
// mismatch_nn_allawi_peyret
//   durum   : transcribed_unverified
//   sha256  : 730a53555fe1d24c90c80bb595b2a7884be5d3e6d3f74e10dc79d6ec08848c7a
// dna_loop_length_santalucia_hicks_2004
//   durum   : transcribed_unverified
//   sha256  : 1b557f6ed7a8bbf938c1dd16bce9486cf36d8ad9aeece63091996efa796d6b81
//
// Anahtar sozlesmesi: "XY/ZW" -> ust iplik 5'-XY-3',
// alt iplik 3'-ZW-5'. Ikinci sutun uyumsuz ciftir. Anahtarlar
// donme esdegerinin leksikografik kucugudur.

const std::map<std::string, std::pair<double, double>>& mismatch_table() {
    static const std::map<std::string, std::pair<double, double>> t = {
        {"AA/AT", {4.7, 12.9}}, {"AA/CT", {7.6, 20.2}}, {"AA/GT", {3.0, 7.4}},
        {"AA/TA", {1.2, 1.7}}, {"AA/TC", {2.3, 4.6}}, {"AA/TG", {-0.6, -2.3}},
        {"AC/AG", {-2.9, -9.8}}, {"AC/CG", {-0.7, -3.8}}, {"AC/GG", {0.5, 3.2}},
        {"AC/TA", {5.3, 14.6}}, {"AC/TC", {0.0, -4.4}}, {"AC/TT", {0.7, 0.2}},
        {"AG/AC", {-0.9, -4.2}}, {"AG/CC", {0.6, -0.6}}, {"AG/GC", {-4.0, -13.2}},
        {"AG/TA", {-0.7, -2.3}}, {"AG/TG", {-3.1, -9.5}}, {"AG/TT", {1.0, 0.9}},
        {"AT/TC", {-1.2, -6.2}}, {"AT/TG", {-2.5, -8.3}}, {"AT/TT", {-2.7, -10.8}},
        {"CA/AT", {3.4, 8.0}}, {"CA/CT", {6.1, 16.4}}, {"CA/GC", {1.9, 3.7}},
        {"CA/GG", {-0.7, -2.3}}, {"CA/TT", {1.0, 0.7}}, {"CC/AG", {5.2, 14.2}},
        {"CC/CG", {3.6, 8.9}}, {"CC/GC", {-1.5, -7.2}}, {"CC/GT", {-0.8, -4.5}},
        {"CC/TG", {5.2, 13.5}}, {"CG/GG", {-4.9, -15.3}}, {"CG/GT", {-4.1, -11.7}},
        {"CG/TC", {-1.5, -6.1}}, {"CT/GG", {-2.8, -8.0}}, {"CT/GT", {-5.0, -15.8}},
        {"GA/AT", {0.7, 0.7}}, {"GA/CG", {-0.6, -1.0}}, {"GA/GT", {1.6, 3.6}},
        {"GA/TT", {-1.3, -5.3}}, {"GC/CT", {2.3, 5.4}}, {"GC/GG", {-6.0, -15.8}},
        {"GC/TG", {-4.4, -12.3}}, {"GG/CT", {3.3, 10.4}}, {"GT/CT", {-2.2, -8.4}},
        {"TA/CT", {1.2, 0.7}}, {"TA/GT", {-0.1, -1.7}}, {"TA/TT", {0.2, -1.5}},
    };
    return t;
}

// Tek bazlik bulge uzunluk terimi (dG37, kcal/mol).
// Tablo yalniz dG37 tasir; kaynak varsayimi dH = 0'dir, bu yuzden
// terim tamamen entropiye yuklenir (align.py ile ayni sozlesme).
constexpr double kBulge1_dG37 = 4.0;

// Bulge'u kapatan A/T baz cifti basina ek ceza (dG37, kcal/mol).
// Kaynak: params.MODEL_CONSTANTS['at_bulge_closure_dg37']
// Literatur tablosu DEGIL, beyan edilmis proje modelleme sabiti;
// registry_digest()e girer.
constexpr double kAtBulgeClosure_dG37 = 0.5;

