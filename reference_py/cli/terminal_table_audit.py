#!/usr/bin/env python3
"""Terminal uyumsuzluk tablosundaki kaynak farklarinin CT-04'e etkisi.

NEDEN. `params.TERMINAL_MISMATCH_SANTALUCIA_PEYRET` iki girdide bagimsiz
kaynakla (primer3 tstack2) ayrilir ve bu ayriliklar `SOURCE_NOTE` olarak
kayda alinmistir:

  GG/CG dH : tablo -6.9149 kcal/mol  (basili dG37 = -0.96 ve dS = -19.2'den
             turetilmis), primer3 -0.7000 kcal/mol (= basili dH).
  AC/TC dS : tablo +0.5 cal/(mol K), primer3 -0.5 cal/(mol K) (isaret).

Bu betik ikisinin de CT-04 ciktisina etkisini olcer: her terminal anahtari
icin sentetik bir ARMS baglami kurar (son baz uyumsuz, sondan ikinci baz
Watson-Crick), CT-04'u tablo degeriyle ve oracle degeriyle hesaplar, farki
ve karar degisimini yazar.

Kullanim:
    python3 cli/terminal_table_audit.py [--out out/terminal_audit]
                                        [--ta 60] [--prefix ACGTTGCAAGGTCCATTC]

Cikti: terminal_audit.tsv, terminal_audit_summary.md
Cikis kodu her zaman 0; bu bir olcum betigidir, kapi degildir.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import engine, params, profiles, salt  # noqa: E402

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# (anahtar, buyukluk, tablo degeri, oracle degeri) -- crosscheck SOURCE_NOTE'lari
SOURCE_NOTES = (
    ("GG/CG", "dH", -6.9149, -0.7000),
    ("AC/TC", "dS", 0.5000, -0.5000),
)


def contexts(prefix: str) -> list[tuple[str, str, str]]:
    """Her terminal anahtari icin (anahtar, primer, uyumsuz alel) uret."""
    out = []
    for key in sorted(params.TERMINAL_MISMATCH_SANTALUCIA_PEYRET.require()):
        top, bottom = key.split("/")
        x, y = top[0], top[1]
        z, w = bottom[0], bottom[1]
        if COMPLEMENT[x] != z:          # sondan ikinci cift WC olmali
            continue
        if COMPLEMENT[y] == w:          # son cift uyumsuz olmali
            continue
        primer = prefix + x + y
        mismatched = prefix + x + COMPLEMENT[w]
        out.append((key, primer, mismatched))
    return out


def ct04_values(eng, key: str, primer: str, mismatched: str) -> dict:
    result = eng.evaluate_ct04(f"AUDIT_{key.replace('/', '_')}", primer,
                               primer, mismatched, "+").as_dict()
    return result


def pick(result: dict, *names):
    """Kayit semasindan ilk bulunan alani dondur (sema esnekligi icin)."""
    for name in names:
        if name in result:
            return result[name]
    for value in result.values():
        if isinstance(value, dict):
            for name in names:
                if name in value:
                    return value[name]
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "out", "terminal_audit"))
    ap.add_argument("--prefix", default="ACGTTGCAAGGTCCATTC",
                    help="sentetik primerin sabit 5' bolumu")
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    args = ap.parse_args()

    cond = salt.Conditions(
        monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
        dntp_M=args.dntp / 1e3, primer_M=args.primer_nm / 1e9,
        anneal_c=args.ta,
    )
    eng = engine.ThermoEngine(cond, profiles.get("ctdna").qc, "owczarzy_2008")
    table = params.TERMINAL_MISMATCH_SANTALUCIA_PEYRET.data

    cases = contexts(args.prefix.upper())
    rows = []
    changed = []
    for key, primer, mismatched in cases:
        base = ct04_values(eng, key, primer, mismatched)
        note = next((n for n in SOURCE_NOTES if n[0] == key), None)
        if note is None:
            rows.append({
                "key": key, "primer_3p": primer[-2:], "allele_3p": mismatched[-1],
                "degisen": "hayir", "buyukluk": "", "tablo": "", "oracle": "",
                "d_theta_tablo": pick(base, "theta_margin"),
                "d_theta_oracle": "", "ddG37_tablo": pick(base, "ddG37_kcal_mol",
                                                          "delta_dG37_kcal_mol",
                                                          "ddG37"),
                "ddG37_oracle": "", "karar_tablo": base.get("severity",
                                                            base.get("verdict", "")),
                "karar_oracle": "", "fark_ddG37": "",
            })
            continue

        _, quantity, ours, theirs = note
        index = 0 if quantity == "dH" else 1
        original = table[key][index]
        table[key][index] = theirs
        try:
            alt = ct04_values(eng, key, primer, mismatched)
        finally:
            table[key][index] = original

        dd_base = pick(base, "delta_delta_G37_kcal_mol")
        dd_alt = pick(alt, "delta_delta_G37_kcal_mol")
        th_base = pick(base, "theta_margin")
        th_alt = pick(alt, "theta_margin")
        sev_base = base.get("severity", base.get("verdict", ""))
        sev_alt = alt.get("severity", alt.get("verdict", ""))
        rows.append({
            "key": key, "primer_3p": primer[-2:], "allele_3p": mismatched[-1],
            "degisen": "evet", "buyukluk": quantity,
            "tablo": f"{ours:+.4f}", "oracle": f"{theirs:+.4f}",
            "d_theta_tablo": th_base, "d_theta_oracle": th_alt,
            "ddG37_tablo": dd_base, "ddG37_oracle": dd_alt,
            "karar_tablo": sev_base, "karar_oracle": sev_alt,
            "fark_ddG37": (None if dd_base is None or dd_alt is None
                           else round(dd_base - dd_alt, 4)),
        })
        changed.append(rows[-1])

    os.makedirs(args.out, exist_ok=True)
    tsv_path = os.path.join(args.out, "terminal_audit.tsv")
    fields = ["key", "primer_3p", "allele_3p", "degisen", "buyukluk", "tablo",
              "oracle", "d_theta_tablo", "d_theta_oracle", "ddG37_tablo",
              "ddG37_oracle", "karar_tablo", "karar_oracle", "fark_ddG37"]
    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    md_path = os.path.join(args.out, "terminal_audit_summary.md")
    with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# Terminal uyumsuzluk tablosu - kaynak farkinin CT-04'e etkisi\n\n")
        fh.write(f"- sentetik baglam: {len(rows)}\n")
        fh.write(f"- kaynak farki tasiyan baglam: {len(changed)}\n")
        fh.write(f"- kosul: Ta {args.ta} C, Na {args.na} mM, Mg {args.mg} mM, "
                 f"primer {args.primer_nm} nM\n\n")
        fh.write("| Anahtar | Buyukluk | Tablo | Oracle | ddG37 (tablo) | "
                 "ddG37 (oracle) | Fark | Karar (tablo -> oracle) |\n")
        fh.write("|---|---|---:|---:|---:|---:|---:|---|\n")
        for row in changed:
            fh.write(f"| {row['key']} | {row['buyukluk']} | {row['tablo']} | "
                     f"{row['oracle']} | {row['ddG37_tablo']} | "
                     f"{row['ddG37_oracle']} | {row['fark_ddG37']} | "
                     f"{row['karar_tablo']} -> {row['karar_oracle']} |\n")

    print(f"sentetik baglam          : {len(rows)}")
    print(f"kaynak farki tasiyan     : {len(changed)}")
    for row in changed:
        print(f"  {row['key']} {row['buyukluk']}: tablo={row['tablo']} "
              f"oracle={row['oracle']} | ddG37 {row['ddG37_tablo']} -> "
              f"{row['ddG37_oracle']} (fark {row['fark_ddG37']}) | "
              f"karar {row['karar_tablo']} -> {row['karar_oracle']}")
    print("ciktilar:")
    for path in (tsv_path, md_path):
        print(f"  {path}  ({os.path.getsize(path)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
