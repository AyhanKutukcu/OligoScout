#!/usr/bin/env python3
"""Primer paneli -> termodinamik QC raporu + erime egrileri.

Girdi TSV'si iki bicimi de kabul eder:
  (a) primer_id <tab> sequence
  (b) pair_id <tab> forward_primer <tab> reverse_primer [... diger sutunlar]
      -> her satir <pair_id>_F ve <pair_id>_R olarak acilir

Cikti (varsayilan out/ klasoru):
  thermo_report.tsv      her primer icin duzlestirilmis ozet
  thermo_records.jsonl   her primer icin tam REP-02 kanit kaydi
  melting_curves.tsv     primer_id, temp_c, theta
  param_registry.tsv     parametre tablolari, durumlari ve SHA-256'lari

Ornek:
  python3 cli/thermo_report.py --panel ../OligoScout/data/benchmarks/grch38_chr22_panel_32.tsv \
      --ta 60 --mg 3.0 --primer-nm 400 --out out
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import engine, nn, params, profiles, salt  # noqa: E402

DEMO_PANEL = [
    ("DEMO_01_F", "GTCACTGGCATTCAAGAGCTACAT"),
    ("DEMO_01_R", "ACCTGTGACATTGCCTGGAGATTA"),
    ("DEMO_02_F", "TGCAAGGTCTTACCGAATCGATCA"),
    ("DEMO_02_R", "GGCTTAACGTTGCACAAGTCTTCA"),
    ("DEMO_03_F", "CAAGTGCCTAAGGTCATCGGATTC"),
    ("DEMO_03_R", "TTGCACCAGTGACTTAGCCATGAA"),
    ("GC_RICH", "GCGCGGCGGCGGCGCGGCGC"),
    ("AT_RICH", "ATTATTAATTATTAATTATT"),
    ("SELF_COMP", "GGGGCCCC"),
    ("SHORT", "ACGTACGTACGT"),
]


def read_panel(path: str):
    """(primer listesi, cift listesi) dondurur."""
    rows: list[tuple[str, str]] = []
    pairs: list[tuple[str, str, str]] = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        fwd = next((cols[k] for k in ("forward_primer", "primer1", "left_primer",
                                      "fwd", "forward") if k in cols), None)
        rev = next((cols[k] for k in ("reverse_primer", "primer2", "right_primer",
                                      "rev", "reverse") if k in cols), None)
        if fwd and rev:
            idc = cols.get("pair_id") or (reader.fieldnames or ["id"])[0]
            for r in reader:
                pid = r[idc]
                rows.append((f"{pid}_F", r[fwd]))
                rows.append((f"{pid}_R", r[rev]))
                pairs.append((pid, r[fwd], r[rev]))
        elif "sequence" in cols:
            idc = cols.get("primer_id") or (reader.fieldnames or ["id"])[0]
            for r in reader:
                rows.append((r[idc], r[cols["sequence"]]))
        else:
            raise SystemExit(
                "TSV basligi taninmadi; 'primer_id/sequence' veya "
                "'pair_id/forward_primer/reverse_primer' bekleniyor"
            )
    return rows, pairs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", help="primer TSV'si; verilmezse gomulu demo panel")
    ap.add_argument("--out", default="out")
    ap.add_argument("--ct04-alleles",
                    help="TSV: primer_id, primer, strand, matched_plus, mismatched_plus")
    ap.add_argument("--ta", type=float, default=60.0, help="annealing sicakligi, C")
    ap.add_argument("--na", type=float, default=50.0, help="tek degerlikli katyon, mM")
    ap.add_argument("--mg", type=float, default=1.5, help="Mg2+, mM")
    ap.add_argument("--dntp", type=float, default=0.6, help="toplam dNTP, mM")
    ap.add_argument("--primer-nm", type=float, default=200.0, help="primer, nM")
    ap.add_argument("--salt-model", default="owczarzy_2008",
                    choices=list(salt.MODELS))
    ap.add_argument("--tm-window", default="58:62", help="profil Tm penceresi, C")
    ap.add_argument("--curve-min", type=float, default=20.0)
    ap.add_argument("--curve-max", type=float, default=95.0)
    ap.add_argument("--curve-step", type=float, default=0.25)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    lo, hi = (float(x) for x in args.tm_window.split(":"))

    cond = salt.Conditions(
        monovalent_M=args.na / 1e3,
        magnesium_M=args.mg / 1e3,
        dntp_M=args.dntp / 1e3,
        primer_M=args.primer_nm / 1e9,
        anneal_c=args.ta,
    )
    prof = engine.Profile(tm_window_c=(lo, hi))
    eng = engine.ThermoEngine(cond=cond, profile=prof, salt_model=args.salt_model)
    ct04 = []
    if args.ct04_alleles:
        ct_engine = engine.ThermoEngine(cond, profiles.get("ctdna").qc, args.salt_model)
        with open(args.ct04_alleles, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            keys = ("primer_id", "primer", "matched_plus", "mismatched_plus", "strand")
            if not set(keys).issubset(reader.fieldnames or ()):
                ap.error(f"CT-04 TSV sütunları: {', '.join(keys)}")
            for row in reader:
                try:
                    ct04.append(ct_engine.evaluate_ct04(*(row[k] for k in keys)).as_dict())
                except (ValueError, KeyError, AttributeError) as exc:
                    ap.error(f"CT-04 TSV satır {reader.line_num}: {exc}")
        if not ct04:
            ap.error("CT-04 TSV dosyası boş")

    if args.panel:
        panel, pairs = read_panel(args.panel)
    else:
        panel, pairs = DEMO_PANEL, []

    # --- parametre kayit defteri ------------------------------------------
    reg_path = os.path.join(args.out, "param_registry.tsv")
    with open(reg_path, "w", newline="", encoding="utf-8") as fh:
        rows = params.registry_summary()
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # --- raporlar ----------------------------------------------------------
    rep_path = os.path.join(args.out, "thermo_report.tsv")
    rec_path = os.path.join(args.out, "thermo_records.jsonl")
    crv_path = os.path.join(args.out, "melting_curves.tsv")

    header = [
        "primer_id", "sequence", "length_nt", "gc_pct", "tm_c", "tm_1M_c",
        "dH_kcal_mol", "dS_cal_mol_K", "dG37_kcal_mol", "theta_at_ta",
        "dG37_last5_kcal_mol", "salt_spread_c", "selfdimer_dG37",
        "selfdimer3p_dG37", "conc_model", "salt_regime",
        "worst_severity", "flags",
    ]

    def _dg(check):
        v = check["value"]
        return v.get("dG37_kcal_mol", "") if v.get("found", True) else ""
    n_fail = n_warn = 0
    skipped: list[tuple[str, str]] = []

    with open(rep_path, "w", newline="", encoding="utf-8") as fr, \
         open(rec_path, "w", encoding="utf-8") as fj, \
         open(crv_path, "w", newline="", encoding="utf-8") as fc:
        rw = csv.writer(fr, delimiter="\t", lineterminator="\n")
        rw.writerow(header)
        cw = csv.writer(fc, delimiter="\t", lineterminator="\n")
        cw.writerow(["primer_id", "temp_c", "theta"])

        for pid, seq in panel:
            try:
                nn.validate(seq)
            except ValueError as exc:
                skipped.append((pid, str(exc)))
                continue
            r = eng.evaluate(pid, seq)
            by_id = {c["id"]: c for c in r["checks"]}
            flags = [f"{c['id']}:{c['severity']}" for c in r["checks"]
                     if c["severity"] != "INFO"]
            if r["worst_severity"] == "FAIL":
                n_fail += 1
            elif r["worst_severity"] == "WARN":
                n_warn += 1

            rw.writerow([
                pid, r["sequence"], r["length_nt"],
                round(100 * r["gc_fraction"], 2),
                r["thermo"]["tm_c"], r["thermo"]["tm_1M_c"],
                r["thermo"]["dH_kcal_mol"], r["thermo"]["dS_1M_cal_mol_K"],
                r["thermo"]["dG37_kcal_mol"],
                by_id["TH-03"]["value"]["theta_at_ta"],
                next(v for k, v in by_id["TH-04"]["value"].items()
                     if k.startswith("dG37_last")),
                by_id["TH-02"]["value"]["spread_c"],
                _dg(by_id["DIM-01"]), _dg(by_id["DIM-02"]),
                r["thermo"]["conc_model"], r["thermo"]["salt_regime"],
                r["worst_severity"], ";".join(flags) or "-",
            ])
            fj.write(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n")
            for t, th in eng.melting_curve(seq, args.curve_min, args.curve_max,
                                           args.curve_step):
                cw.writerow([pid, t, round(th, 6)])

    # --- cift duzeyi rapor (DIM-03, DIM-04, PAIR-01) -----------------------
    pair_path = os.path.join(args.out, "pair_report.tsv")
    n_pair_fail = n_pair_warn = 0
    if pairs:
        with open(pair_path, "w", newline="", encoding="utf-8") as fp, \
             open(rec_path, "a", encoding="utf-8") as fj:
            pw = csv.writer(fp, delimiter="\t", lineterminator="\n")
            pw.writerow(["pair_id", "delta_tm_c", "tm_f_c", "tm_r_c",
                         "heterodimer_dG37", "heterodimer3p_dG37",
                         "worst_severity", "flags"])
            for pid, f_seq, r_seq in pairs:
                try:
                    pr = eng.evaluate_pair(pid, f_seq, r_seq)
                except ValueError:
                    continue
                bid = {c["id"]: c for c in pr["checks"]}
                if pr["worst_severity"] == "FAIL":
                    n_pair_fail += 1
                elif pr["worst_severity"] == "WARN":
                    n_pair_warn += 1
                pw.writerow([
                    pid, bid["PAIR-01"]["value"]["delta_tm_c"],
                    bid["PAIR-01"]["value"]["tm_f_c"],
                    bid["PAIR-01"]["value"]["tm_r_c"],
                    _dg(bid["DIM-03"]), _dg(bid["DIM-04"]),
                    pr["worst_severity"],
                    ";".join(f"{c['id']}:{c['severity']}" for c in pr["checks"]
                             if c["severity"] != "INFO") or "-",
                ])
                fj.write(json.dumps(pr, ensure_ascii=False, sort_keys=True) + "\n")

    if ct04:
        with open(rec_path, "a", encoding="utf-8") as fj, \
             open(os.path.join(args.out, "ct04_records.jsonl"), "w", encoding="utf-8") as fc, \
             open(os.path.join(args.out, "ct04.tsv"), "w", newline="", encoding="utf-8") as ft:
            writer = csv.DictWriter(ft, fieldnames=list(ct04[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in ct04:
                writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, dict) else v
                                 for k, v in row.items()})
                fc.write(json.dumps(row, ensure_ascii=False) + "\n")
                fj.write(json.dumps({"target": row["target"], "checks": [row],
                                     "worst_severity": row["severity"], "engine_version": row["engine_version"]},
                                    ensure_ascii=False) + "\n")
        print(f"CT-04 kayıt sayısı   : {len(ct04)}")

    n = len(panel) - len(skipped)
    print(f"primer sayisi        : {n}")
    if pairs:
        print(f"cift sayisi          : {len(pairs)}")
        print(f"cift FAIL / WARN     : {n_pair_fail} / {n_pair_warn}")
    print(f"FAIL / WARN          : {n_fail} / {n_warn}")
    for pid, why in skipped:
        print(f"  atlandi {pid}: {why}")
    print(f"parametre digest     : {params.registry_digest()[:16]}...")
    print(f"motor surumu         : {engine.ENGINE_VERSION}")
    print("ciktilar:")
    outputs = [rep_path, rec_path, crv_path, reg_path]
    if pairs:
        outputs.insert(1, pair_path)
    for p in outputs:
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
