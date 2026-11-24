#!/usr/bin/env python3
"""LNA modelinin McTigue 2004 Tablo 1 verisine karsi dogrulamasi.

Iki yontem (makalenin "Comparison Among Different Methods" bolumu):
  M3  olculen referans DNA dH/dS + LNA tablosu  -> dTm   (tablonun kendisini sinar;
      makale 32NN combined icin RMSD(dTm)=0.88 C bildirir)
  M4  OligoScout DNA birlesik NN + LNA tablosu, 1.011 M Na+ -> mutlak Tm
      (uygulamadaki yol; makale 16(RY)XL(RY) ile RMSD(Tm)=1.47 C bildirir)

Girdi : reference_py/validation/mctigue2004_table1.tsv
Cikti : out/lna_validation/{per_oligo.tsv, summary.json}
"""
from __future__ import annotations

import csv
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

from oligoscout_thermo import concentration, lna, oligomod, salt  # noqa: E402
from oligoscout_thermo.params import R_CAL  # noqa: E402

DATA = os.path.join(HERE, "validation", "mctigue2004_table1.tsv")
COND = salt.Conditions(monovalent_M=1.011, magnesium_M=0.0, dntp_M=0.0)


def load():
    with open(DATA, encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def tm_eq5(dH, dS, ct):
    return dH * 1000.0 / (dS + R_CAL * math.log(ct / 4.0)) - 273.15


def rmsd(xs):
    return math.sqrt(sum(x * x for x in xs) / len(xs)) if xs else float("nan")


def run():
    rows = load()
    ref = {r["name"]: r for r in rows if r["reference"] == "-"}
    out = []
    for r in rows:
        if r["reference"] == "-":
            continue
        d = ref[r["reference"]]
        p = oligomod.parse(r["sequence"])
        ct = float(r["ct_M"])
        inc = lna.increments(p.dna, p.lna_positions)
        m3 = tm_eq5(float(d["dH_kcal"]) + inc.ddH_kcal, float(d["dS_cal"]) + inc.ddS_cal, ct)
        m3_ref = tm_eq5(float(d["dH_kcal"]), float(d["dS_cal"]), ct)
        m4 = concentration.build(p.dna, COND, "equal", "owczarzy_2004", ct,
                                 lna_positions=p.lna_positions).tm_c
        m4_ref = concentration.build(p.dna, COND, "equal", "owczarzy_2004", ct).tm_c
        out.append({
            "name": r["name"], "set": r["set"], "sequence": r["sequence"], "flag": r["flag"],
            "tm_exp": float(r["tm_c"]), "dtm_exp": float(r["tm_c"]) - float(d["tm_c"]),
            "dtm_m3": m3 - m3_ref, "tm_m4": m4, "dtm_m4": m4 - m4_ref, "ref_tm_exp": float(d["tm_c"]),
            "ref_tm_m4": m4_ref, "in_domain": inc.in_domain,
        })
    return out


def summarize(out):
    use = [o for o in out if not o["flag"]]
    e3 = [o["dtm_m3"] - o["dtm_exp"] for o in use]
    e4 = [o["tm_m4"] - o["tm_exp"] for o in use]
    e4d = [o["dtm_m4"] - o["dtm_exp"] for o in use]
    refs = {o["sequence"].replace("+", ""): 0 for o in use}
    eref = []
    seen = set()
    for o in use:
        key = (o["ref_tm_exp"], round(o["ref_tm_m4"], 6))
        if key not in seen:
            seen.add(key)
            eref.append(o["ref_tm_m4"] - o["ref_tm_exp"])
    by_set = {}
    for s in ("training", "test"):
        u = [o for o in use if o["set"] == s]
        by_set[s] = {"n": len(u),
                     "rmsd_dtm_m3": round(rmsd([o["dtm_m3"] - o["dtm_exp"] for o in u]), 3),
                     "rmsd_tm_m4": round(rmsd([o["tm_m4"] - o["tm_exp"] for o in u]), 3)}
    return {
        "source": "McTigue, Peterson & Kahn, Biochemistry 2004, Table 1 (100 LNA oligos, 24 DNA references)",
        "excluded": [o["name"] for o in out if o["flag"]],
        "n": len(use),
        "M3_measured_reference_plus_table": {
            "rmsd_dtm_c": round(rmsd(e3), 3), "mean_c": round(statistics.mean(e3), 3),
            "max_abs_c": round(max(map(abs, e3)), 3), "paper_reported_rmsd_dtm_c": 0.88},
        "M4_oligoscout_unified_dna_plus_table": {
            "rmsd_tm_c": round(rmsd(e4), 3), "mean_c": round(statistics.mean(e4), 3),
            "max_abs_c": round(max(map(abs, e4)), 3), "within_2c": sum(abs(x) <= 2 for x in e4),
            "rmsd_dtm_c": round(rmsd(e4d), 3),
            "dna_reference_rmsd_tm_c": round(rmsd(eref), 3),
            "paper_reported_rmsd_tm_c_16RY_test_set": 1.47},
        "by_set": by_set,
        "conditions": {"na_M": COND.monovalent_M, "salt_model": "owczarzy_2004",
                       "concentration": "CT/4 (5 uM; 1 uM for LT4/LT5)"},
    }


def main():
    out = run()
    summ = summarize(out)
    dst = os.path.join(HERE, "out", "lna_validation")
    os.makedirs(dst, exist_ok=True)
    with open(os.path.join(dst, "per_oligo.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        for o in out:
            w.writerow({k: (round(v, 3) if isinstance(v, float) else v) for k, v in o.items()})
    with open(os.path.join(dst, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
