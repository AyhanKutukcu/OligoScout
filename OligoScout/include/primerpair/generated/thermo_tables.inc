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

// dangling_ends_bommarito_2000
//   durum   : transcribed_unverified
//   sha256  : ec8243ce524eba1ca5848182bcd55d683fefd8a0414575f4254eb719ea949474
// Anahtar: "D3:K:P-Q" / "D5:K:P-Q": sarkan K, uc cifti P-Q (P sarkan iplikte)
const std::map<std::string, std::pair<double, double>>& dangling_table() {
    static const std::map<std::string, std::pair<double, double>> t = {
        {"D3:A:A-T", {-0.5, -1.1}}, {"D3:A:C-G", {-5.9, -16.5}}, {"D3:A:G-C", {-2.1, -3.9}},
        {"D3:A:T-A", {-0.7, -0.8}}, {"D3:C:A-T", {4.7, 14.2}}, {"D3:C:C-G", {-2.6, -7.4}},
        {"D3:C:G-C", {-0.2, -0.1}}, {"D3:C:T-A", {4.4, 14.9}}, {"D3:G:A-T", {-4.1, -13.1}},
        {"D3:G:C-G", {-3.2, -10.4}}, {"D3:G:G-C", {-3.9, -11.2}}, {"D3:G:T-A", {-1.6, -3.6}},
        {"D3:T:A-T", {-3.8, -12.6}}, {"D3:T:C-G", {-5.2, -15.0}}, {"D3:T:G-C", {-4.4, -13.1}},
        {"D3:T:T-A", {2.9, 10.4}}, {"D5:A:A-T", {0.2, 2.3}}, {"D5:A:C-G", {-6.3, -17.1}},
        {"D5:A:G-C", {-3.7, -10.0}}, {"D5:A:T-A", {-2.9, -7.6}}, {"D5:C:A-T", {0.6, 3.3}},
        {"D5:C:C-G", {-4.4, -12.6}}, {"D5:C:G-C", {-4.0, -11.9}}, {"D5:C:T-A", {-4.1, -13.0}},
        {"D5:G:A-T", {-1.1, -1.6}}, {"D5:G:C-G", {-5.1, -14.0}}, {"D5:G:G-C", {-3.9, -10.9}},
        {"D5:G:T-A", {-4.2, -15.0}}, {"D5:T:A-T", {-6.9, -20.0}}, {"D5:T:C-G", {-4.0, -10.9}},
        {"D5:T:G-C", {-4.9, -13.8}}, {"D5:T:T-A", {-0.2, -0.5}},
    };
    return t;
}

// terminal_mismatch_santalucia_peyret_2001
//   durum   : transcribed_unverified
//   sha256  : acee72a6f38f16d526952250a37bf43b2e979c529e0ec899c8002ca58fb7c929
// Anahtar: "XY/ZW": ust 5'-XY-3', alt 3'-ZW-5', Y-W uyumsuz (heliks ucu)
const std::map<std::string, std::pair<double, double>>& terminal_mismatch_table() {
    static const std::map<std::string, std::pair<double, double>> t = {
        {"AA/TA", {-3.1, -7.8}}, {"AA/TC", {-1.6, -4.0}}, {"AA/TG", {-1.9, -4.4}},
        {"AC/TA", {-1.8, -3.8}}, {"AC/TC", {-0.1, 0.5}}, {"AC/TT", {-0.9, -1.7}},
        {"AG/TA", {-2.5, -5.9}}, {"AG/TG", {-1.1, -2.1}}, {"AG/TT", {-3.2, -8.7}},
        {"AT/TC", {-2.3, -6.3}}, {"AT/TG", {-3.5, -9.4}}, {"AT/TT", {-2.4, -6.5}},
        {"CA/GA", {-4.3, -10.7}}, {"CA/GC", {-2.6, -5.9}}, {"CA/GG", {-3.9, -9.6}},
        {"CC/GA", {-2.7, -6.0}}, {"CC/GC", {-2.1, -5.1}}, {"CC/GT", {-3.2, -8.0}},
        {"CG/GA", {-6.0, -15.5}}, {"CG/GG", {-3.8, -9.5}}, {"CG/GT", {-3.8, -9.0}},
        {"CT/GC", {-3.9, -10.6}}, {"CT/GG", {-6.6, -18.7}}, {"CT/GT", {-6.1, -16.9}},
        {"GA/CA", {-8.0, -22.5}}, {"GA/CC", {-5.0, -13.8}}, {"GA/CG", {-4.3, -11.1}},
        {"GC/CA", {-3.2, -7.1}}, {"GC/CC", {-3.9, -10.6}}, {"GC/CT", {-4.9, -13.5}},
        {"GG/CA", {-4.6, -11.4}}, {"GG/CG", {-6.914879999999999, -19.2}}, {"GG/CT", {-5.7, -15.9}},
        {"GT/CC", {-3.0, -7.8}}, {"GT/CG", {-5.9, -16.1}}, {"GT/CT", {-7.4, -21.2}},
        {"TA/AA", {-2.5, -6.3}}, {"TA/AC", {-2.3, -5.9}}, {"TA/AG", {-2.0, -4.7}},
        {"TC/AA", {-2.7, -7.0}}, {"TC/AC", {-0.7, -1.3}}, {"TC/AT", {-2.5, -6.3}},
        {"TG/AA", {-2.4, -5.8}}, {"TG/AG", {-1.1, -2.7}}, {"TG/AT", {-3.9, -10.5}},
        {"TT/AC", {-0.7, -1.2}}, {"TT/AG", {-3.6, -9.8}}, {"TT/AT", {-3.2, -8.9}},
    };
    return t;
}

