#!/usr/bin/env python3
"""OligoScout thermodynamics vs. primer3 executables (oligotm, ntthal).

Sequence-level benchmark on the Test #81 panel (2,048 primers, 1,024 pairs).
Complements reference_py/cli/param_crosscheck.py (table-level) and
salt_crosscheck.py (10 sequences).

Comparisons
  A  Tm, 1 M Na+, concentration convention ALIGNED
     OligoScout tm_1M_c (C = primer conc.)  vs  oligotm -tp 1 -sc 2 -d 4*C
     -> isolates the nearest-neighbour layer (should agree to print precision).
  B  Tm as each tool is used by default (Na 50, Mg 1.5, dNTP 0.6 mM, 200 nM)
     OligoScout owczarzy_2008 (C)  vs  oligotm defaults (-sc 1, C/4)
     -> the offset a user actually sees.
  C  Tm with the same salt model and aligned C
     OligoScout owczarzy_2008 (Owczarzy 2004/2008 with the 0.22 crossover)  vs  oligotm -sc 2 -d 4*C
     (primer3's -sc 2 is the same 2004/2008 combination, see oligotm.c)
  D  hetero-dimer dG37 (1 M Na+)   OligoScout align.best_dimer(any)  vs  ntthal -a ANY
  E  self-dimer  dG37 (1 M Na+)    same, a == b
  F  hairpin dG37 (1 M Na+)        OligoScout fold.zuker_mfe        vs  ntthal -a HAIRPIN

Usage (WSL / Linux, primer3 built from third_party/primer3/src):
  python3 thermo_vs_primer3.py --primer3 ~/p3/src --panel ../test81/inputs/test81_panel_1024.tsv \
      --out results/thermo_vs_primer3 [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics as st
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                      # AYHAN-OligoScout
sys.path.insert(0, str(ROOT / "reference_py"))
from oligoscout_thermo import align, engine, fold, salt  # noqa: E402

NA, MG, DNTP, PNM = 50.0, 1.5, 0.6, 200.0   # mM, mM, mM, nM


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout


def oligotm(p3, seq, mv, dv, n, d_nM, sc):
    out = run([str(p3 / "oligotm"), "-tp", "1", "-sc", str(sc), "-mv", str(mv), "-dv", str(dv),
               "-n", str(n), "-d", str(d_nM), seq]).strip()
    return float(out)


def ntthal(p3, a, b, mode):
    cmd = [str(p3 / "ntthal"), "-path", str(p3 / "primer3_config") + "/", "-mv", "1000", "-dv", "0",
           "-n", "0", "-d", str(PNM), "-t", "37", "-a", mode, "-s1", a]
    if b is not None:
        cmd += ["-s2", b]
    out = run(cmd)
    line = out.splitlines()[0] if out else ""
    if "No secondary structure" in out or "dG" not in line:
        return 0.0
    dg = float(line.split("dG = ")[1].split()[0])
    return dg / 1000.0


def oligoscout_tm(seq, model, cond):
    r = engine.ThermoEngine(cond=cond, salt_model=model).evaluate("x", seq)["thermo"]
    return r


def stats(xs, ys):
    d = [x - y for x, y in zip(xs, ys)]
    n = len(d)
    out = {"n": n, "mean_diff": st.fmean(d), "sd_diff": st.pstdev(d) if n > 1 else 0.0,
           "mae": st.fmean(abs(v) for v in d), "max_abs": max(abs(v) for v in d)}
    if n > 2 and st.pstdev(xs) > 0 and st.pstdev(ys) > 0:
        out["pearson_r"] = st.correlation(xs, ys)
        rx = _ranks(xs); ry = _ranks(ys)
        out["spearman_rho"] = st.correlation(rx, ry)
    return {k: round(v, 4) if isinstance(v, float) else v for k, v in out.items()}


def _ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2
        i = j + 1
    return r


def agreement(xs, ys, thr):
    """Flag agreement at the same dG threshold (dG < thr -> flagged)."""
    tp = sum(1 for x, y in zip(xs, ys) if x < thr and y < thr)
    fp = sum(1 for x, y in zip(xs, ys) if x < thr <= y)
    fn = sum(1 for x, y in zip(xs, ys) if y < thr <= x)
    tn = len(xs) - tp - fp - fn
    return {"threshold_kcal_mol": thr, "both_flag": tp, "oligoscout_only": fp, "primer3_only": fn,
            "neither": tn, "agreement": round((tp + tn) / len(xs), 4),
            "cohen_kappa": round(_kappa(tp, fp, fn, tn), 4)}


def _kappa(tp, fp, fn, tn):
    n = tp + fp + fn + tn
    po = (tp + tn) / n
    pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n / n
    return (po - pe) / (1 - pe) if pe < 1 else 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primer3", default=os.path.expanduser("~/p3/src"))
    ap.add_argument("--panel", default=str(HERE.parent / "test81" / "inputs" / "test81_panel_1024.tsv"))
    ap.add_argument("--out", default=str(HERE / "results" / "thermo_vs_primer3"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    p3 = Path(a.primer3)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    pairs = list(csv.DictReader(open(a.panel), delimiter="\t"))
    if a.limit:
        pairs = pairs[: a.limit]
    primers = [(p["primer1_id"], p["primer1"]) for p in pairs] + [(p["primer2_id"], p["primer2"]) for p in pairs]
    t0 = time.time()

    c1m = salt.Conditions(monovalent_M=1.0, magnesium_M=0.0, dntp_M=0.0, primer_M=PNM / 1e9, anneal_c=60.0)
    cpcr = salt.Conditions(monovalent_M=NA / 1e3, magnesium_M=MG / 1e3, dntp_M=DNTP / 1e3,
                           primer_M=PNM / 1e9, anneal_c=60.0)

    def per_primer(item):
        pid, s = item
        os_def = oligoscout_tm(s, "owczarzy_2008", cpcr)
        return {
            "id": pid, "seq": s, "len": len(s),
            "A_os": os_def["tm_1M_c"], "A_p3": oligotm(p3, s, 1000, 0, 0, 4 * PNM, 2),
            "B_os": os_def["tm_c"], "B_p3": oligotm(p3, s, NA, MG, DNTP, PNM, 1),
            "C_os": os_def["tm_c"], "C_p3": oligotm(p3, s, NA, MG, DNTP, 4 * PNM, 2),
            "E_os": align.best_dimer(s, s, "any").as_dict().get("dG37_kcal_mol", 0.0)
            if align.best_dimer(s, s, "any").found else 0.0,
            "E_p3": ntthal(p3, s, s, "ANY"),
            "F_os": fold.zuker_mfe(s)[0], "F_p3": ntthal(p3, s, None, "HAIRPIN"),
        }

    def per_pair(p):
        d = align.best_dimer(p["primer1"], p["primer2"], "any")
        return {"pair_id": p["pair_id"], "D_os": d.as_dict()["dG37_kcal_mol"] if d.found else 0.0,
                "D_p3": ntthal(p3, p["primer1"], p["primer2"], "ANY")}

    with ThreadPoolExecutor(a.jobs) as ex:
        prim = list(ex.map(per_primer, primers))
        pr = list(ex.map(per_pair, pairs))

    for name, rows in (("per_primer.tsv", prim), ("per_pair.tsv", pr)):
        with open(out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)

    col = lambda rows, k: [r[k] for r in rows]  # noqa: E731
    summary = {
        "primers": len(prim), "pairs": len(pr), "seconds": round(time.time() - t0, 1),
        "conditions": {"na_mM": NA, "mg_mM": MG, "dntp_mM": DNTP, "primer_nM": PNM},
        "primer3_version_hint": "third_party/primer3/src (built in WSL, gcc -O2)",
        "A_tm_1M_aligned_conc": stats(col(prim, "A_os"), col(prim, "A_p3")),
        "B_tm_default_each_tool": stats(col(prim, "B_os"), col(prim, "B_p3")),
        "C_tm_same_salt_model_aligned_conc": stats(col(prim, "C_os"), col(prim, "C_p3")),
        "D_heterodimer_dG37": {**stats(col(pr, "D_os"), col(pr, "D_p3")),
                               "flag_agreement": agreement(col(pr, "D_os"), col(pr, "D_p3"), -6.0)},
        "E_selfdimer_dG37": {**stats(col(prim, "E_os"), col(prim, "E_p3")),
                             "flag_agreement": agreement(col(prim, "E_os"), col(prim, "E_p3"), -6.0)},
        "F_hairpin_dG37": {**stats(col(prim, "F_os"), col(prim, "F_p3")),
                           "flag_agreement": agreement(col(prim, "F_os"), col(prim, "F_p3"), -2.0)},
        "notes": [
            "A isolates the NN layer: OligoScout uses C (primer excess), primer3 C/4; oligotm is given -d 4*C.",
            "B is what users see: different salt model (Owczarzy 2008 vs SantaLucia 1998) AND concentration convention.",
            "D-F: ntthal includes dangling ends, terminal mismatches and full loop tables; OligoScout dimers are "
            "WC stacks + <=1 internal mismatch or one 1-nt bulge, so OligoScout dG is expected to be less negative.",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
