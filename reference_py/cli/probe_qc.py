#!/usr/bin/env python3
"""Tasarlanan ciftler icin hidroliz probu secer ve PRB-01..PRB-05 uygular.

Zincirdeki yeri:
    design_primers.py  -> cift adaylari (bu betigin girdisi)
    probe_qc.py        -> her cift icin prob + PRB kanit kayitlari   <-- burasi
    oligoscout_panel   -> ozgulluk taramasi (proba gerek duymaz)

Girdi : --design klasoru (design_pairs.tsv + design_panel.tsv) ve --target
        FASTA (amplikon dizisi buradan kesilir).
        Hazir problari denetlemek icin: --probes pair_id<TAB>probe TSV'si.
Cikti : probe_panel.tsv, probe_qc.tsv, probe_qc_records.jsonl

Esikler literatur kilavuz degeridir, KALIBRE DEGILDIR (REP-03).
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import nn, probe as prb, profiles, salt  # noqa: E402

SEVERITY_ORDER = {"INFO": 0, "WARN": 1, "FAIL": 2}


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path: str) -> str:
    parts = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if parts:
                    break
                continue
            parts.append(line.strip())
    return nn.validate("".join(parts))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True,
                    help="design_pairs.tsv + design_panel.tsv iceren klasor")
    ap.add_argument("--target", required=True, help="tasarimda kullanilan FASTA")
    ap.add_argument("--out", help="cikti klasoru (varsayilan: --design)")
    ap.add_argument("--probes", help="hazir problar: pair_id<TAB>probe TSV'si")
    ap.add_argument("--profile", default="ddpcr", choices=sorted(profiles.PROFILES))
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    ap.add_argument("--probe-length", default="13:30", help="min:max nt")
    ap.add_argument("--tm-delta", default="6:10", help="PRB-01 penceresi, C")
    args = ap.parse_args()

    out_dir = args.out or args.design
    os.makedirs(out_dir, exist_ok=True)

    target = read_fasta(args.target)
    pairs = read_tsv(os.path.join(args.design, "design_pairs.tsv"))
    panel = read_tsv(os.path.join(args.design, "design_panel.tsv"))

    lmin, lmax = (int(x) for x in args.probe_length.split(":"))
    dmin, dmax = (float(x) for x in args.tm_delta.split(":"))
    cons = prb.ProbeConstraints(length_min=lmin, length_max=lmax,
                                tm_delta_min=dmin, tm_delta_max=dmax)
    cond = salt.Conditions(
        monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
        dntp_M=args.dntp / 1e3, primer_M=args.primer_nm / 1e9,
        anneal_c=args.ta,
    )

    given: dict[str, str] = {}
    if args.probes:
        for row in read_tsv(args.probes):
            given[row["pair_id"]] = row["probe"]

    panel_rows: list[list] = []
    qc_rows: list[list] = []
    records: list[dict] = []
    rejections: collections.Counter = collections.Counter()
    verdicts: collections.Counter = collections.Counter()
    no_probe: list[str] = []

    # design_panel.tsv, design_pairs.tsv'nin ilk N satiriyla ayni siradadir.
    for row, src in zip(panel, pairs):
        if row["primer1"] != src["left_primer"]:
            raise SystemExit("design_panel ve design_pairs eslesmiyor; "
                             "ikisini ayni kosudan alin")
        pair_id = row["pair_id"]
        fwd, rev = src["left_primer"], src["right_primer"]
        amp = target[int(src["product_start"]):int(src["product_end"])]

        if pair_id in given:
            seq = nn.validate(given[pair_id])
            cand = None
            source = "verilen"
        else:
            cands, rej = prb.enumerate_probes(amp, fwd, rev, cond, cons)
            for rule, _ in rej:
                rejections[rule] += 1
            if not cands:
                no_probe.append(pair_id)
                verdicts["PROB YOK"] += 1
                panel_rows.append([pair_id, "", "", "", "", "", "", "",
                                   "", "", "PROB YOK", ""])
                continue
            cand = cands[0]
            seq = cand.sequence
            source = "tasarlandi"

        res = prb.evaluate_probe(f"{pair_id}_P", seq, fwd, rev, amplicon=amp,
                                 cond=cond, cons=cons,
                                 profile_name=args.profile)
        verdicts[res["worst_severity"]] += 1
        records.append(res)

        start = end = strand = ""
        for ch in res["checks"]:
            if ch["id"] == "PRB-04" and ch["value"].get("found_in_amplicon"):
                start = ch["value"]["start"]
                end = ch["value"]["end"]
                strand = ch["value"]["strand"]
        g, c = prb.gc_counts(seq)
        panel_rows.append([
            pair_id, seq, strand, start, end, len(seq),
            round(100 * nn.gc_fraction(seq), 2), res["tm_c"],
            res["delta_tm_c"], f"G{g}/C{c}", res["worst_severity"], source,
        ])

        for ch in res["checks"]:
            qc_rows.append([
                pair_id, ch["id"], ch["target"], ch["severity"],
                json.dumps(ch["value"], ensure_ascii=False,
                           separators=(",", ":")),
                json.dumps(ch["threshold"], ensure_ascii=False,
                           separators=(",", ":")),
                ch["note"],
            ])

    panel_path = os.path.join(out_dir, "probe_panel.tsv")
    with open(panel_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pair_id", "probe", "strand", "start_in_amplicon",
                    "end_in_amplicon", "length_nt", "gc_pct", "tm_c",
                    "delta_tm_to_primers_c", "g_c_counts", "verdict",
                    "source"])
        w.writerows(panel_rows)

    qc_path = os.path.join(out_dir, "probe_qc.tsv")
    with open(qc_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pair_id", "check", "target", "severity", "value",
                    "threshold", "note"])
        w.writerows(qc_rows)

    rec_path = os.path.join(out_dir, "probe_qc_records.jsonl")
    with open(rec_path, "w", encoding="utf-8", newline="\n") as fh:
        for r in records:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")

    worst_by_check = collections.Counter(
        (r[1], r[3]) for r in qc_rows if r[3] != "INFO")
    print(f"cift sayisi         : {len(panel)}")
    print(f"prob secildi        : {len(records)}"
          + (f", prob bulunamadi: {len(no_probe)}" if no_probe else ""))
    print("karar               : " +
          ", ".join(f"{k}={v}" for k, v in sorted(verdicts.items())))
    if records:
        dts = [r["delta_tm_c"] for r in records]
        print(f"PRB-01 dTm araligi  : {min(dts):.2f} .. {max(dts):.2f} C "
              f"(pencere {dmin:g}-{dmax:g} C)")
        worst = min(records,
                    key=lambda r: (r["worst_dimer"] or {}).get(
                        "dG37_kcal_mol", 0.0))
        if worst["worst_dimer"]:
            print(f"PRB-05 en kotu dimer: "
                  f"{worst['worst_dimer']['dG37_kcal_mol']:.3f} kcal/mol "
                  f"({worst['probe_id']} x {worst['worst_dimer']['partner']})")
    if worst_by_check:
        print("WARN/FAIL dagilimi  : " +
              ", ".join(f"{c}:{s}={n}" for (c, s), n in
                        sorted(worst_by_check.items())))
    if rejections:
        print("eleme gerekceleri   : " +
              ", ".join(f"{k}={v}" for k, v in rejections.most_common()))
    if no_probe:
        print("prob bulunamayan    : " + ", ".join(no_probe))
    print("esikler kalibre degildir (REP-03)")
    print("ciktilar:")
    for p in (panel_path, qc_path, rec_path):
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
