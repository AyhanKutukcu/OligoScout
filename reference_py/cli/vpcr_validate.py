#!/usr/bin/env python3
"""Kalendar et al. 2024 (virtualPCR) insan genomu sonuclarinin OligoScout ile tekrari.

Girdi: validation/kalendar2024_virtualpcr_results.tsv, validation/kalendar2024_human_genes.tsv
Calistirir: cli/virtual_pcr.py --genome grch38 (iki deney) ve karsilastirir.
Cikti: out/vpcr_validation/summary.json
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAL = os.path.join(HERE, "validation")
OUT = os.path.join(HERE, "out", "vpcr_validation")


def read(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))


def run(name, args):
    out = os.path.join(OUT, name)
    subprocess.run([sys.executable, os.path.join(HERE, "cli", "virtual_pcr.py"), "--genome", "grch38",
                    "--out", out, *args], check=True)
    return read(os.path.join(out, "products.tsv")), json.load(open(os.path.join(out, "summary.json")))


def match(ref, prods):
    found = []
    for r in ref:
        hit = [p for p in prods if p["record"] == r["chrom"] and int(p["start"]) + 1 == int(r["start"])
               and int(p["end_exclusive"]) == int(r["end"])]
        found.append({"id": r["id"], "region": f"{r['chrom']}:{r['start']}-{r['end']}", "found": bool(hit),
                      "orientation": hit[0]["orientation"] if hit else None,
                      "mismatches": [int(hit[0]["left_mismatches"]), int(hit[0]["right_mismatches"])] if hit else None,
                      "ta_virtualpcr": float(r["ta_c"]), "ta_oligoscout": float(hit[0]["ta_c"]) if hit and hit[0]["ta_c"] else None})
    return found


def main():
    ref = read(os.path.join(VAL, "kalendar2024_virtualpcr_results.tsv"))
    os.makedirs(OUT, exist_ok=True)
    genes, gs = run("six_genes", ["--oligos", os.path.join(VAL, "kalendar2024_human_genes.tsv"),
                                  "--min-product", "100", "--max-product", "1000",
                                  "--terminal-window", "7", "--terminal-max", "0"])
    six = match([r for r in ref if r["experiment"] == "six_genes"], genes)
    rt, rs = run("rt_copia", ["--primer", "RTp=CARATGGAYGTNAARAC", "--primer", "RTm=CATRTCRTCNACRTA",
                              "--mode", "all", "--min-product", "200", "--max-product", "500",
                              "--terminal-window", "7", "--terminal-max", "1", "--max-mismatch", "3", "--no-thermo"])
    rtm = match([r for r in ref if r["experiment"] == "rt_copia"], rt)
    summary = {
        "source": "Kalendar et al. 2024 Front. Bioinform. 4:1464197; Zenodo 13625500 virtualPCR outputs",
        "six_genes": {"virtualpcr_products": 6, "oligoscout_products": len(genes),
                      "identical_coordinates": sum(x["found"] for x in six), "pairs": six,
                      "rule": "<=3 mismatches/primer, 0 in the 3'-terminal 7 nt (number3errors=0), 100-1000 bp"},
        "rt_copia_degenerate": {"virtualpcr_products": 11, "found_by_oligoscout": sum(x["found"] for x in rtm),
                                "oligoscout_products_same_rule": len(rt), "products": rtm,
                                "rule": "IUPAC primers, <=3 mismatches/primer, <=1 in the 3'-terminal 7 nt, 200-500 bp, all primer combinations",
                                "note": ("virtualPCR does not document a total-mismatch cap per primer; its 11 products "
                                         "carry 0-3 mismatches per site and are a subset of the OligoScout set under the "
                                         "stated rule")},
        "ta_note": ("Ta formula identical (min primer Tm + ln length); absolute values differ because the primer Tm "
                    "models differ (virtualPCR/FastPCR vs OligoScout NN at 50 mM Na+, 1.5 mM Mg2+, 200 nM)."),
    }
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "source"}, indent=1)[:1500])


if __name__ == "__main__":
    main()
