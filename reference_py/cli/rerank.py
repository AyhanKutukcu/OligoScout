#!/usr/bin/env python3
"""Tasarim siralamasini ozgulluk taramasina gore yeniden duzenler.

Zincirin son halkasi:
    design_primers.py  -> termodinamik siralama (ozgullukten habersiz)
    oligoscout_panel   -> tum genomda eksiksiz tarama
    panel_qc.py        -> hedef / hedef disi siniflandirma
    variant_qc.py      -> VAR-01..VAR-08 (istege bagli, --variants)
    rerank.py          -> nihai karar ve siralama          <-- bu dosya

Karar kurali (hepsi beyan edilmis, kalibre edilmemis):
    FAIL  hedef tam olarak bir kez bulunmuyor (SPC-01), ya da
          tam eslesmeli hedef disi urun var (SPC-03)
    WARN  yalniz uyumsuzluklu hedef disi urun var
    PASS  genom genelinde tek urun

--variants verilirse varyant katmani da karara girer:
    FAIL  herhangi bir VAR kontrolu FAIL (orn. 3' ucta yaygin varyant)
    WARN  herhangi bir VAR kontrolu WARN
    PASS  temiz
Nihai karar iki katmanin EN KOTUSUDUR; ozgullukte PASS olup varyantta FAIL
alan bir cift nihai FAIL'dir. Iki katman birbirinin yerine gecmez: ozgulluk
"baska yerde urun var mi", varyant "hedefte baglanma her bireyde tutuyor mu"
sorusunu yanitlar.

Siralama: once nihai karar (PASS < WARN < FAIL), sonra hedef disi urun sayisi,
sonra varyant WARN/FAIL sayisi, sonra termodinamik skor.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo.variants import (  # noqa: E402
    VERDICT_ORDER, combine_verdicts, verdict_from_severities,
)


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True, help="design_* ciktilarinin klasoru")
    ap.add_argument("--qc", required=True, help="panel_qc.tsv iceren klasor")
    ap.add_argument("--out", help="cikti klasoru (varsayilan: --qc)")
    ap.add_argument("--variants", help="variant_qc.tsv iceren klasor "
                                       "(verilmezse varyant katmani atlanir)")
    args = ap.parse_args()

    out_dir = args.out or args.qc
    os.makedirs(out_dir, exist_ok=True)

    pairs = read_tsv(os.path.join(args.design, "design_pairs.tsv"))
    panel = read_tsv(os.path.join(args.design, "design_panel.tsv"))
    qc = {r["pair_id"]: r for r in read_tsv(os.path.join(args.qc, "panel_qc.tsv"))}

    # --- varyant katmani (istege bagli) ---------------------------------
    var_verdict: dict[str, str] = {}
    var_flags: dict[str, int] = {}
    var_driver: dict[str, str] = {}
    if args.variants:
        import json
        vpath = os.path.join(args.variants, "variant_qc.tsv")
        if not os.path.exists(vpath):
            raise SystemExit(f"varyant tablosu yok: {vpath}\n"
                             "Once: RUN_VARIANT.cmd <tarama> <tasarim>")
        best_af: dict[str, float] = {}
        for r in read_tsv(vpath):
            pid = r["pair_id"]
            cur = var_verdict.get(pid, "PASS")
            sev = r["severity"]
            if sev != "INFO":
                var_flags[pid] = var_flags.get(pid, 0) + 1
                value = json.loads(r["value"])
                af = float(value.get("af", value.get("max_af", 0)) or 0)
                if af > best_af.get(pid, 0.0):
                    best_af[pid] = af
                    var_driver[pid] = (f'{value.get("variant", "")} '
                                       f'{r["check"]} AF={af:g}').strip()
            var_verdict[pid] = verdict_from_severities([
                "WARN" if cur == "WARN" else "FAIL" if cur == "FAIL" else "INFO",
                sev])

    # design_panel.tsv, design_pairs.tsv'nin ilk N satiriyla ayni siradadir.
    rows = []
    for row, src in zip(panel, pairs):
        if row["primer1"] != src["left_primer"]:
            raise SystemExit("design_panel ve design_pairs eslesmiyor; "
                             "ikisini ayni kosudan alin")
        q = qc.get(row["pair_id"])
        if q is None:
            continue

        total = int(q["products_total"])
        intended = int(q["products_intended"])
        off = int(q["offtarget_total"])
        off_exact = int(q["offtarget_exact"])

        if intended != 1 or off_exact > 0:
            verdict = "FAIL"
        elif off > 0:
            verdict = "WARN"
        else:
            verdict = "PASS"

        vv = var_verdict.get(row["pair_id"], "")
        combined = combine_verdicts(verdict, vv)

        rows.append({
            "pair_id": row["pair_id"],
            "verdict": combined,
            "specificity_verdict": verdict,
            "variant_verdict": vv or "degerlendirilmedi",
            "variant_flags": var_flags.get(row["pair_id"], 0),
            "variant_driver": var_driver.get(row["pair_id"], ""),
            "left_primer": row["primer1"],
            "right_primer": row["primer2"],
            "product_length": src["product_length"],
            "tm_left_c": src["tm_left_c"],
            "tm_right_c": src["tm_right_c"],
            "delta_tm_c": src["delta_tm_c"],
            "theta_f_at_ta": q["theta_f_at_ta"],
            "theta_r_at_ta": q["theta_r_at_ta"],
            "products_total": total,
            "products_intended": intended,
            "offtarget_total": off,
            "offtarget_exact": off_exact,
            "thermo_rank": src["pareto_rank"],
            "thermo_score": src["uncalibrated_score"],
        })

    rows.sort(key=lambda r: (VERDICT_ORDER[r["verdict"]], r["offtarget_total"],
                             r["variant_flags"], float(r["thermo_score"])))
    for i, r in enumerate(rows, 1):
        r["final_rank"] = i

    path = os.path.join(out_dir, "final_ranked.tsv")
    cols = ["final_rank", "verdict", "specificity_verdict", "variant_verdict",
            "variant_flags", "variant_driver",
            "pair_id", "left_primer", "right_primer",
            "product_length", "tm_left_c", "tm_right_c", "delta_tm_c",
            "theta_f_at_ta", "theta_r_at_ta", "products_total",
            "products_intended", "offtarget_total", "offtarget_exact",
            "thermo_rank", "thermo_score"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    import collections
    counts = collections.Counter(r["verdict"] for r in rows)
    print(f"cift sayisi : {len(rows)}")
    print(f"nihai karar : " +
          ", ".join(f"{k}={counts.get(k, 0)}" for k in ("PASS", "WARN", "FAIL")))
    if args.variants:
        spec = collections.Counter(r["specificity_verdict"] for r in rows)
        varc = collections.Counter(r["variant_verdict"] for r in rows)
        print("  ozgulluk  : " +
              ", ".join(f"{k}={spec.get(k, 0)}" for k in ("PASS", "WARN", "FAIL")))
        print("  varyant   : " +
              ", ".join(f"{k}={varc.get(k, 0)}" for k in ("PASS", "WARN", "FAIL")))
        flipped = [r for r in rows
                   if r["specificity_verdict"] == "PASS" and r["verdict"] != "PASS"]
        print(f"  ozgullukte PASS, varyantta dusen: {len(flipped)} cift"
              + (f" ({', '.join(r['pair_id'] for r in flipped[:5])})" if flipped else ""))

    # Termodinamik siralamanin ozgullugu ne kadar yanlis tahmin ettigi
    moved = [r for r in rows if abs(r["final_rank"] - int(r["thermo_rank"] or 0)) > 0]
    best_thermo = min(rows, key=lambda r: float(r["thermo_score"]))
    print(f"termodinamik en iyi : {best_thermo['pair_id']} "
          f"({best_thermo['verdict']}, {best_thermo['offtarget_total']} hedef disi)")
    if rows:
        print(f"nihai en iyi        : {rows[0]['pair_id']} "
              f"({rows[0]['verdict']}, {rows[0]['offtarget_total']} hedef disi, "
              f"{rows[0]['product_length']} bp)")
    print(f"cikti: {path}  ({os.path.getsize(path)} bayt)")
    print("skor kalibre degildir (REP-03)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
