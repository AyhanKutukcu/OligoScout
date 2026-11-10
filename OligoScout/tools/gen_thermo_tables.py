#!/usr/bin/env python3
"""params.py tablolarini C++ basligina cevirir (TEK DOGRULUK KAYNAGI).

NEDEN URETILIYOR. Parametreler deneysel veridir ve depoda tek bir yerde
durur: `reference_py/oligoscout_thermo/params.py`. C++ cekirdegine elle
ikinci kez yazmak yeni bir transkripsiyon hatasi riskidir. Bu betik sayilari
oradan okuyup `include/primerpair/generated/thermo_tables.inc` dosyasina
yazar; boylece iki uygulama ayni sayilari paylasir ve `cli/diff_cpp.py`
farkin modelden geldigini gosterebilir.

Kullanim (depo kokunden):
    python3 tools/gen_thermo_tables.py
    python3 tools/gen_thermo_tables.py --check    # fark varsa cikis 1

`--check` surumu CI/regresyon icindir: uretilmis dosya params.py ile
tutmuyorsa hata verir.
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
ROOT = os.path.dirname(REPO)
sys.path.insert(0, os.path.join(ROOT, "reference_py"))

from oligoscout_thermo import params  # noqa: E402

OUT_PATH = os.path.join(REPO, "include", "primerpair", "generated",
                        "thermo_tables.inc")


def rotate(key: str) -> str:
    top, bottom = key.split("/")
    return f"{bottom[1]}{bottom[0]}/{top[1]}{top[0]}"


def render() -> str:
    mismatch = params.MISMATCH_ALLAWI_PEYRET
    loops = params.DNA_LOOP_SANTALUCIA_HICKS_2004
    data = mismatch.require()

    # params.py anahtarlari kanonik (donme esdegerinin leksikografik kucugu)
    # olmali; `nn._mismatch_stack_adjustment` bu varsayimla arama yapar.
    for key in data:
        canonical = min(key, rotate(key))
        if key != canonical:
            raise SystemExit(
                f"{key}: params.py anahtari kanonik degil (beklenen {canonical})"
            )

    bulge1 = loops.require()["BULGE:1"][0]
    # align.py'deki proje modelleme sabiti. Parametre kayit defterinde
    # DEGILDIR (kunye/SHA-256 tasimaz); buraya tek kaynaktan tasinir ki
    # C++ ile Python ayrismasin.
    at_closure = params.AT_BULGE_CLOSURE.value

    lines = [
        "// OTOMATIK URETILDI -- ELLE DUZENLEMEYIN.",
        "// Ureten : tools/gen_thermo_tables.py",
        "// Kaynak : reference_py/oligoscout_thermo/params.py",
        f"// {mismatch.provenance.table_id}",
        f"//   durum   : {mismatch.provenance.verification_status}",
        f"//   sha256  : {mismatch.sha256}",
        f"// {loops.provenance.table_id}",
        f"//   durum   : {loops.provenance.verification_status}",
        f"//   sha256  : {loops.sha256}",
        "//",
        "// Anahtar sozlesmesi: \"XY/ZW\" -> ust iplik 5'-XY-3',",
        "// alt iplik 3'-ZW-5'. Ikinci sutun uyumsuz ciftir. Anahtarlar",
        "// donme esdegerinin leksikografik kucugudur.",
        "",
        "const std::map<std::string, std::pair<double, double>>& "
        "mismatch_table() {",
        "    static const std::map<std::string, std::pair<double, double>> t = {",
    ]
    items = sorted(data.items())
    for index in range(0, len(items), 3):
        chunk = items[index:index + 3]
        row = " ".join(
            '{"%s", {%s, %s}},' % (k, repr(v[0]), repr(v[1])) for k, v in chunk
        )
        lines.append("        " + row)
    lines += [
        "    };",
        "    return t;",
        "}",
        "",
        "// Tek bazlik bulge uzunluk terimi (dG37, kcal/mol).",
        "// Tablo yalniz dG37 tasir; kaynak varsayimi dH = 0'dir, bu yuzden",
        "// terim tamamen entropiye yuklenir (align.py ile ayni sozlesme).",
        f"constexpr double kBulge1_dG37 = {bulge1!r};",
        "",
        "// Bulge'u kapatan A/T baz cifti basina ek ceza (dG37, kcal/mol).",
        "// Kaynak: params.MODEL_CONSTANTS['at_bulge_closure_dg37']",
        "// Literatur tablosu DEGIL, beyan edilmis proje modelleme sabiti;",
        "// registry_digest()e girer.",
        f"constexpr double kAtBulgeClosure_dG37 = {at_closure!r};",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="dosyayi yazma, params.py ile tuttugunu dogrula")
    ap.add_argument("--out", default=OUT_PATH)
    args = ap.parse_args()

    text = render()
    if args.check:
        if not os.path.exists(args.out):
            print(f"uretilmis dosya yok: {args.out}", file=sys.stderr)
            return 1
        with open(args.out, encoding="utf-8") as fh:
            current = fh.read()
        if current != text:
            print(f"{args.out} params.py ile tutmuyor; "
                  "tools/gen_thermo_tables.py calistirin", file=sys.stderr)
            return 1
        print(f"guncel: {args.out}")
        return 0

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print(f"yazildi: {args.out} ({os.path.getsize(args.out)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
