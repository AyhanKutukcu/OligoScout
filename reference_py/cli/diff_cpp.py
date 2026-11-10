#!/usr/bin/env python3
"""Python referans katmani ile C++ `native_thermo` cekirdeginin diferansiyel testi.

Zincirdeki yeri:
    oligoscout_thermo_dump   -> C++ cekirdegin TSV dokumu (src/thermo_dump_main.cpp)
    diff_cpp.py              -> ayni dizileri Python'dan gecirir, karsilastirir  <-- burasi

Kullanim:
    python3 cli/diff_cpp.py --sequences girdi.tsv --dump cpp_dump.tsv
    python3 cli/diff_cpp.py --sequences girdi.tsv \\
        --binary ../OligoScout/build/release/oligoscout_thermo_dump

Girdi TSV basligi: id<TAB>sequence[<TAB>partner]

Cikti: diff_cpp.tsv (satir basina bir id x sutun), diff_cpp_summary.md,
       diff_cpp_records.jsonl

SUTUN SINIFLARI
  cekirdek  : dH, dS, n_bp, self_complementary, gc_fraction, tm_1M_c, tm_c,
              theta_ta, dG37, terminal_dG37_5, salt_regime, conc_model.
              Iki uygulama ayni modeli kodluyor; fark = HATA.
  dimer     : selfdimer_*, hetero_*. Python tarafinda tek bazlik bulge
              uretim yolundadir (align.MODEL_NAME). C++ cekirdek bu
              genislemeyi tasimiyorsa fark MODEL BOSLUGUDUR, hata degil.
              Bu yuzden her dimer sutunu Python'un IKI varyantiyla
              (bulge acik / kapali) karsilastirilir ve hangisinin tuttugu
              yazilir.

Cikis kodu 1: en az bir CEKIRDEK sutun toleransi asiyor.
Cikis kodu 0: cekirdek sutunlar tutuyor (dimer boslugu raporlanir).

Zaman damgasi yazilmaz; ayni girdi ayni dosyayi uretir (REP-01).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import align, concentration, engine, salt  # noqa: E402

CORE_COLUMNS = (
    "dH_kcal", "dS_cal", "n_bp", "self_complementary", "gc_fraction",
    "tm_1M_c", "tm_c", "salt_regime", "conc_model", "theta_ta",
    "dG37_kcal", "terminal_dG37_5",
)
DIMER_COLUMNS = (
    "selfdimer_any_dG37", "selfdimer_3p_dG37",
    "hetero_any_dG37", "hetero_3p_dG37", "hetero_any_helix",
)
TEXT_COLUMNS = ("salt_regime", "conc_model")

DEFAULT_TOL = {
    "dH_kcal": 1e-6, "dS_cal": 1e-6, "n_bp": 0.0, "self_complementary": 0.0,
    "gc_fraction": 1e-9, "tm_1M_c": 1e-6, "tm_c": 1e-6, "theta_ta": 1e-9,
    "dG37_kcal": 1e-6, "terminal_dG37_5": 1e-6,
    "salt_regime": 0.0, "conc_model": 0.0,
    "selfdimer_any_dG37": 1e-6, "selfdimer_3p_dG37": 1e-6,
    "hetero_any_dG37": 1e-6, "hetero_3p_dG37": 1e-6, "hetero_any_helix": 0.0,
}


def read_sequences(path: str) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for cols in reader:
            if not cols or not cols[0].strip():
                continue
            if cols[0].strip() == "id":
                continue
            if len(cols) < 2:
                raise SystemExit(f"{path}: 'id<TAB>sequence' bekleniyor: {cols}")
            partner = cols[2].strip() if len(cols) > 2 else ""
            rows.append((cols[0].strip(), cols[1].strip().upper(),
                         partner.upper()))
    if not rows:
        raise SystemExit(f"{path}: dizi yok")
    return rows


def read_dump(path: str) -> dict[str, dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(CORE_COLUMNS + DIMER_COLUMNS) - set(reader.fieldnames or ())
        if missing:
            raise SystemExit(f"{path}: eksik sutunlar: {sorted(missing)}")
        return {row["id"]: row for row in reader}


def run_binary(binary: str, sequences: str, args) -> str:
    cmd = [binary, "--sequences", sequences,
           "--na", str(args.na), "--mg", str(args.mg), "--dntp", str(args.dntp),
           "--primer-nm", str(args.primer_nm), "--ta", str(args.ta),
           "--salt-model", args.salt_model]
    done = subprocess.run(cmd, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"{binary} cikis kodu {done.returncode}: {done.stderr[:500]}")
    return done.stdout


def python_side(pid: str, seq: str, partner: str, eng, cond,
                allow_bulge: bool) -> dict[str, float | str]:
    """Bir dizinin Python tarafindaki degerleri, C++ dump sutun adlariyla.

    YUVARLAMA. `eng.evaluate()` kaydi rapor icin yuvarlanmis degerler tasir
    (Tm 3 ondalik, theta 6 ondalik...). Diferansiyel testte yuvarlama
    gercek bir fark gibi gorunur: C++ dokumu 17 anlamli basamak basar.
    Bu yuzden buyuklukler kaydin degil, `concentration.build` modelinin
    HAM alanlarindan okunur -- C++ tarafindaki `primerpair::thermo::build`
    ile birebir ayni nesne.
    """
    model = concentration.build(seq, cond=cond, salt_model=eng.salt_model)
    duplex = model.duplex
    th04 = engine._terminal_dG37(seq, 5)

    self_any = align.self_dimer(seq, "any", allow_bulge=allow_bulge)
    self_3p = align.self_dimer(seq, "3prime", allow_bulge=allow_bulge)
    if partner:
        het_any = align.best_dimer(seq, partner, "any", allow_bulge=allow_bulge)
        het_3p = align.best_dimer(seq, partner, "3prime", allow_bulge=allow_bulge)
    else:
        het_any = het_3p = None

    def dg(result):
        if result is None or not result.found:
            return 0.0
        return result.dG37_kcal

    def kusur(result):
        """Yapinin kusur sinifi: C++ cekirdegi yalniz kesintisiz helis tarar."""
        if result is None or not result.found:
            return ""
        if result.bulge_count:
            return "bulge"
        if result.mismatch_count:
            return "uyumsuzluk"
        return ""

    flags = {
        "selfdimer_any_dG37": kusur(self_any),
        "selfdimer_3p_dG37": kusur(self_3p),
        "hetero_any_dG37": kusur(het_any),
        "hetero_3p_dG37": kusur(het_3p),
        "hetero_any_helix": kusur(het_any),
    }
    return flags, {
        "dH_kcal": duplex.dH_kcal,
        "dS_cal": duplex.dS_cal,
        "n_bp": float(duplex.n_bp),
        "self_complementary": 1.0 if duplex.self_complementary else 0.0,
        "gc_fraction": model.f_gc,
        "tm_1M_c": model.tm_1m_K - concentration.KELVIN,
        "tm_c": model.tm_c,
        "salt_regime": model.salt_detail.get("regime", ""),
        "conc_model": model.conc_model,
        "theta_ta": model.theta(cond.anneal_c),
        "dG37_kcal": model.dG_kcal(37.0),
        "terminal_dG37_5": th04,
        "selfdimer_any_dG37": dg(self_any),
        "selfdimer_3p_dG37": dg(self_3p),
        "hetero_any_dG37": dg(het_any),
        "hetero_3p_dG37": dg(het_3p),
        "hetero_any_helix": float(het_any.helix_len if het_any
                                  and het_any.found else 0),
    }


def compare(column: str, ours, theirs_text: str, tol: float) -> tuple[str, str]:
    """(durum, fark metni). Metin sutunlari birebir karsilastirilir."""
    if column in TEXT_COLUMNS:
        return ("MATCH" if str(ours) == theirs_text.strip() else "DIFF", "")
    try:
        theirs = float(theirs_text)
    except ValueError:
        return "UNPARSED", ""
    delta = float(ours) - theirs
    return ("MATCH" if abs(delta) <= tol else "DIFF", f"{delta:+.6g}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sequences", required=True,
                    help="id<TAB>sequence[<TAB>partner] TSV'si")
    ap.add_argument("--dump", help="C++ dokumunun TSV'si")
    ap.add_argument("--binary", help="oligoscout_thermo_dump yolu (dump yerine)")
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "out", "diff_cpp"))
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--salt-model", default="owczarzy_2008",
                    choices=list(salt.MODELS))
    ap.add_argument("--tol", type=float, default=None,
                    help="tum sayisal sutunlar icin tek tolerans")
    args = ap.parse_args()

    if not args.dump and not args.binary:
        ap.error("--dump ya da --binary verin")

    sequences = read_sequences(args.sequences)
    os.makedirs(args.out, exist_ok=True)

    if args.binary:
        text = run_binary(args.binary, args.sequences, args)
        dump_path = os.path.join(args.out, "cpp_dump.tsv")
        with open(dump_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    else:
        dump_path = args.dump
    cpp = read_dump(dump_path)

    cond = salt.Conditions(
        monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
        dntp_M=args.dntp / 1e3, primer_M=args.primer_nm / 1e9,
        anneal_c=args.ta,
    )
    eng = engine.ThermoEngine(cond=cond, profile=engine.Profile(),
                              salt_model=args.salt_model)

    rows: list[dict[str, str]] = []
    records: list[dict] = []
    core_diffs = 0
    dimer_diffs = 0
    bulge_off_rescues = 0
    missing: list[str] = []

    for pid, seq, partner in sequences:
        if pid not in cpp:
            missing.append(pid)
            continue
        theirs = cpp[pid]
        if theirs["sequence"].strip().upper() != seq:
            raise SystemExit(f"{pid}: dump dizisi girdi dizisiyle ayni degil")
        flags, with_bulge = python_side(pid, seq, partner, eng, cond, True)
        _, without_bulge = python_side(pid, seq, partner, eng, cond, False)

        record = {"id": pid, "sequence": seq, "columns": {}}
        for column in CORE_COLUMNS + DIMER_COLUMNS:
            tol = args.tol if args.tol is not None else DEFAULT_TOL[column]
            if column in TEXT_COLUMNS:
                tol = 0.0
            status, delta = compare(column, with_bulge[column],
                                    theirs[column], tol)
            variant = "bulge_acik"
            if status == "DIFF" and column in DIMER_COLUMNS:
                alt_status, alt_delta = compare(column, without_bulge[column],
                                                theirs[column], tol)
                # NOT: C++ cekirdek artik Python ile AYNI yapi sinifini
                # tariyor (tek ic uyumsuzluk + tek bazlik bulge). Bu yuzden
                # dimer farki da HATADIR; "model boslugu" diye gecistirilmez.
                # Asagidaki bilgi yalniz teshis icindir: Python'un bulge'u
                # kapaliyken tutuyorsa, fark bulge dalindadir.
                if alt_status == "MATCH":
                    variant = "yalniz_bulge_kapaliyken_tutuyor"
                    bulge_off_rescues += 1
                elif flags.get(column):
                    variant = "python_yapisi: " + flags[column]
            if status == "DIFF":
                if column in CORE_COLUMNS:
                    core_diffs += 1
                else:
                    dimer_diffs += 1
            rows.append({
                "id": pid, "column": column, "sinif":
                    "cekirdek" if column in CORE_COLUMNS else "dimer",
                "python": f"{with_bulge[column]}",
                "python_bulge_kapali": f"{without_bulge[column]}",
                "cpp": theirs[column].strip(),
                "delta": delta, "tolerance": f"{tol:g}",
                "eslesen_varyant": variant if status != "DIFF" else "",
                "status": status,
            })
            record["columns"][column] = status
        records.append(record)

    tsv_path = os.path.join(args.out, "diff_cpp.tsv")
    fields = ["id", "column", "sinif", "python", "python_bulge_kapali", "cpp",
              "delta", "tolerance", "eslesen_varyant", "status"]
    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    rec_path = os.path.join(args.out, "diff_cpp_records.jsonl")
    with open(rec_path, "w", encoding="utf-8", newline="\n") as fh:
        for record in records:
            fh.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")

    by_column: dict[str, dict[str, int]] = {}
    for row in rows:
        bucket = by_column.setdefault(row["column"], {})
        bucket[row["status"]] = bucket.get(row["status"], 0) + 1

    md_path = os.path.join(args.out, "diff_cpp_summary.md")
    with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# Python <-> C++ diferansiyel testi\n\n")
        fh.write(f"- dizi sayisi: {len(records)}\n")
        fh.write(f"- cekirdek sutun farki: {core_diffs}\n")
        fh.write(f"- dimer sutun farki: {dimer_diffs}\n")
        fh.write(f"- bulge kapaliyken tutan dimer degeri: {bulge_off_rescues}\n")
        fh.write(f"- kosul: Na {args.na} mM, Mg {args.mg} mM, dNTP {args.dntp} mM, "
                 f"primer {args.primer_nm} nM, Ta {args.ta} C, "
                 f"tuz modeli {args.salt_model}\n")
        if missing:
            fh.write(f"- dump'ta bulunmayan id: {', '.join(missing)}\n")
        fh.write("\n| Sutun | Sinif | MATCH | DIFF |\n")
        fh.write("|---|---|---:|---:|\n")
        for column in CORE_COLUMNS + DIMER_COLUMNS:
            bucket = by_column.get(column, {})
            fh.write(f"| {column} | "
                     f"{'cekirdek' if column in CORE_COLUMNS else 'dimer'} | "
                     f"{bucket.get('MATCH', 0)} | {bucket.get('DIFF', 0)} |\n")
        worst = [r for r in rows if r["status"] == "DIFF"]
        fh.write("\n## Farklar\n\n")
        if worst:
            fh.write("| id | sutun | python | cpp | fark |\n|---|---|---:|---:|---:|\n")
            for row in worst[:50]:
                fh.write(f"| {row['id']} | {row['column']} | {row['python']} | "
                         f"{row['cpp']} | {row['delta']} |\n")
        else:
            fh.write("Yok.\n")

    print(f"dizi sayisi              : {len(records)}")
    print(f"cekirdek sutun farki     : {core_diffs}")
    print(f"dimer sutun farki        : {dimer_diffs}")
    print(f"bulge kapaliyken tutan   : {bulge_off_rescues}")
    if missing:
        print(f"dump'ta yok              : {', '.join(missing)}")
    for row in rows:
        if row["status"] == "DIFF":
            print(f"  FARK [{row['sinif']}] {row['id']} {row['column']}: "
                  f"python={row['python']} cpp={row['cpp']} "
                  f"fark={row['delta']} {row['eslesen_varyant']}")
    print("ciktilar:")
    for path in (tsv_path, md_path, rec_path):
        print(f"  {path}  ({os.path.getsize(path)} bayt)")
    return 1 if (core_diffs or dimer_diffs) else 0


if __name__ == "__main__":
    raise SystemExit(main())
