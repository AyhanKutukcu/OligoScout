#!/usr/bin/env python3
"""params.py tablolarini bagimsiz kaynaklara karsi denetler (PARAM-XCHK).

Zincirdeki yeri:
    param_transcribe.py  -> tablolari birincil makaleden aktarir
    param_crosscheck.py  -> aktarimi bagimsiz kaynakla sinar   <-- burasi
    thermo_report.py     -> kunyeyi her rapor satirina yazar

Kaynaklar (yalnizca OKUNUR; hicbir deger urune kopyalanmaz):
    primer3_thal      primer3_config/{stack,stackmm,dangle,loops}.{dh,ds}
    primer3_oligotm   src/oligotm.c SantaLucia_1998_dH/dS + baslatma terimleri

Cikti : crosscheck.tsv, crosscheck_summary.md, crosscheck_records.jsonl
Cikis kodu 1: en az bir girdi toleransi asiyor ya da eslesmiyor.

Bu denetim bir tabloyu "verified" YAPMAZ. Projenin tanimi geregerli kalir:
"verified" icin birincil PDF'ten ikinci bir kisinin bagimsiz aktarimi ve
diff'i gerekir. Buradaki sonuc kunyeye ayri bir `cross_check` alani olarak
yazilir.

Cikti dosyalarinda zaman damgasi YOKTUR: ayni girdi ayni dosyayi uretir
(REP-01).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import crosscheck as xc, params  # noqa: E402


def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def oracle_files(sources: dict[str, str]) -> list[tuple[str, str]]:
    """(oracle, dosya yolu) listesi; kunyeye checksum yazmak icin."""
    out: list[tuple[str, str]] = []
    cfg = sources.get("config_dir")
    if cfg:
        for name in ("stack.dh", "stack.ds", "stackmm.dh", "stackmm.ds",
                     "dangle.dh", "dangle.ds", "loops.dh", "loops.ds"):
            out.append((xc.ORACLE_THAL, os.path.join(cfg, name)))
    if sources.get("oligotm_c"):
        out.append((xc.ORACLE_OLIGOTM, sources["oligotm_c"]))
    if sources.get("terminal_dir"):
        for name in ("tstack2.dh", "tstack2.ds"):
            out.append((xc.ORACLE_TERMINAL, os.path.join(sources["terminal_dir"], name)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primer3-src",
                    help="primer3 'src' klasoru (varsayilan: depo icindeki "
                         "OligoScout/third_party/primer3/src)")
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "out", "param_crosscheck"),
        help="cikti klasoru")
    ap.add_argument("--all-rows", action="store_true",
                    help="TSV'ye eslesen girdileri de yaz (varsayilan: yazar)")
    args = ap.parse_args()

    sources = xc.oracle_sources(args.primer3_src)
    if not sources:
        print("primer3 kaynagi bulunamadi.", file=sys.stderr)
        print("  --primer3-src ile yol verin ya da OLIGOSCOUT_PRIMER3_SRC "
              "ortam degiskenini ayarlayin.", file=sys.stderr)
        return 2

    entries = xc.run_crosscheck(args.primer3_src)
    summary = xc.summarise(entries)
    os.makedirs(args.out, exist_ok=True)

    tsv_path = os.path.join(args.out, "crosscheck.tsv")
    fields = ["table_id", "oracle", "key", "quantity", "ours", "oracle_value",
              "delta", "tolerance", "status"]
    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for e in entries:
            w.writerow(e.as_row())

    checksums = [
        {"oracle": oracle, "file": os.path.basename(path),
         "sha256": file_sha256(path)}
        for oracle, path in oracle_files(sources)
    ]

    records = []
    for row in summary:
        table = params.REGISTRY[row["table_id"]]
        records.append({
            "check": "PARAM-XCHK",
            "table_id": row["table_id"],
            "table_sha256": table.sha256,
            "verification_status": table.provenance.verification_status,
            "oracle": row["oracle"],
            "oracle_files": [c for c in checksums if c["oracle"] == row["oracle"]],
            "compared": int(row["compared"]),
            "match": int(row["match"]),
            "diff": int(row["diff"]),
            "source_notes": int(row["source_notes"]),
            "missing_ours": int(row["missing_ours"]),
            "missing_oracle": int(row["missing_oracle"]),
            "max_abs_delta": float(row["max_abs_delta"] or 0.0),
            "tolerances": {"dH_kcal": xc.TOL_DH_KCAL, "dS_cal": xc.TOL_DS_CAL,
                           "dG37_kcal": xc.TOL_DG37_KCAL},
            "severity": ("FAIL" if any(int(row[k]) for k in ("diff", "missing_ours", "missing_oracle"))
                         else "WARN" if int(row["source_notes"]) else "INFO"),
            "note": ("diferansiyel denetim; 'verified' durumu icin ikinci "
                     "insan aktarimi hala gereklidir"),
        })

    rec_path = os.path.join(args.out, "crosscheck_records.jsonl")
    with open(rec_path, "w", encoding="utf-8", newline="\n") as fh:
        for r in records:
            fh.write(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n")

    differences = [e for e in entries if e.status != "MATCH"]
    problems = [e for e in differences if e.status != "SOURCE_NOTE"]
    md_path = os.path.join(args.out, "crosscheck_summary.md")
    with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# Parametre tablolari - diferansiyel denetim (PARAM-XCHK)\n\n")
        fh.write("Kaynaklar yalnizca okunur; hicbir deger `params.py`'ye "
                 "kopyalanmaz. Bu denetim transkripsiyon hatasini yakalar, "
                 "birincil PDF'ten ikinci insan aktariminin yerine GECMEZ.\n\n")
        fh.write("| Tablo | Oracle | Karsilastirilan | Eslesen | Sapan | "
                 "Bizde yok | Oracle'da yok | En buyuk mutlak fark |\n")
        fh.write("|---|---|---:|---:|---:|---:|---:|---:|\n")
        for row in summary:
            fh.write("| {table_id} | {oracle} | {compared} | {match} | "
                     "{diff} | {missing_ours} | {missing_oracle} | "
                     "{max_abs_delta} |\n".format(**row))
        fh.write("\n## Toleranslar\n\n")
        fh.write(f"- dH: {xc.TOL_DH_KCAL} kcal/mol\n")
        fh.write(f"- dS: {xc.TOL_DS_CAL} cal/(mol K)\n")
        fh.write(f"- dG37: {xc.TOL_DG37_KCAL} kcal/mol "
                 "(dongu tablosu dH/dS'ten turetilir)\n")
        fh.write("\n## Kaynak dosya checksum'lari\n\n")
        fh.write("| Oracle | Dosya | SHA-256 |\n|---|---|---|\n")
        for c in checksums:
            fh.write(f"| {c['oracle']} | {c['file']} | {c['sha256']} |\n")
        fh.write("\n## Tablo checksum'lari\n\n")
        fh.write("| Tablo | Durum | SHA-256 |\n|---|---|---|\n")
        for key in sorted(params.REGISTRY):
            t = params.REGISTRY[key]
            fh.write(f"| {key} | {t.provenance.verification_status} | "
                     f"{t.sha256} |\n")
        if differences:
            fh.write("\n## Sapan girdiler\n\n")
            fh.write("| Tablo | Oracle | Anahtar | Buyukluk | Bizim | "
                     "Oracle | Fark |\n|---|---|---|---|---:|---:|---:|\n")
            for e in differences:
                r = e.as_row()
                fh.write("| {table_id} | {oracle} | {key} | {quantity} | "
                         "{ours} | {oracle_value} | {delta} |\n".format(**r))
            if any(e.status == "SOURCE_NOTE" for e in differences):
                fh.write("\nSOURCE_NOTE kayıtları. AC/TC dS işareti birincil PDF'deki +0,5 değeridir. "
                         "GG/CG dH, birincil PDF'deki dG37 ve dS'den yeniden türetildi.\n")
        else:
            fh.write("\n## Sapan girdiler\n\nYok.\n")

    for row in summary:
        print(f"{row['table_id']:38s} {row['oracle']:16s} "
              f"karsilastirilan={row['compared']:>4s} eslesen={row['match']:>4s} "
              f"sapan={row['diff']:>3s} eksik={row['missing_ours']}/"
              f"{row['missing_oracle']} en_buyuk_fark={row['max_abs_delta']}")
    if problems:
        print(f"\nSAPMA: {len(problems)} girdi toleransi asiyor ya da eksik.")
        for e in problems[:20]:
            print(f"  {e.status} {e.table_id} {e.oracle} {e.key} "
                  f"{e.quantity} bizim={e.ours} oracle={e.theirs}")
    elif differences:
        print(f"\n{len(differences)} belgelenmiş kaynak farkı (SOURCE_NOTE).")
    else:
        print("\nTum girdiler tolerans icinde.")
    print("bu denetim 'verified' durumu vermez; ikinci insan aktarimi gerekir")
    print("ciktilar:")
    for p in (tsv_path, md_path, rec_path):
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
