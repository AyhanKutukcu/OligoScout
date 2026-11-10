#!/usr/bin/env python3
"""Hedef dizi -> aday primer ciftleri (ozgulluk taramasi ONCESI).

Girdi: FASTA dosyasi ya da dogrudan dizi.
Cikti:
  design_pairs.tsv        siralanmis cift adaylari (Pareto katmani + skor)
  design_candidates.tsv   filtreleri gecen tekil primerler
  design_rejections.tsv   elenen her aday, QC madde kimligi ve degeriyle
  design_panel.tsv        ilk N cift, oligoscout_panel'in okuyacagi bicimde

Sonraki adim: design_panel.tsv dosyasini RUN_PANEL.cmd ile tum genomda
tarayip panel_qc.py ile yeniden siralamak. Skor kalibre degildir.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import design as dsg  # noqa: E402
from oligoscout_thermo import profiles  # noqa: E402
from oligoscout_thermo import nn, salt  # noqa: E402
from oligoscout_thermo.engine import Profile  # noqa: E402


def read_locus(value: str) -> tuple[str, int] | None:
    """FASTA basligindan 'chrN:start-end' okur; bulamazsa None.

    Boylece tasarlanan panel, beklenen genom koordinatlarini tasiyabilir ve
    ozgulluk taramasinda hedef/hedef disi ayrimi kesinlesir.
    """
    if not os.path.exists(value):
        return None
    with open(value, encoding="utf-8") as fh:
        header = fh.readline().strip()
    m = re.match(r">?\s*([\w.]+):(\d+)-(\d+)", header)
    if not m:
        return None
    return m.group(1), int(m.group(2))


def read_target(value: str) -> str:
    looks_like_path = (
        os.sep in value or "/" in value
        or value.lower().endswith((".fa", ".fasta", ".fna", ".txt", ".seq"))
    )
    if looks_like_path and not os.path.exists(value):
        raise SystemExit(
            f"hedef dosyasi bulunamadi: {value}\n"
            "  Tam yol verin ya da diziyi dogrudan yapistirin.\n"
            "  Ornek: RUN_DESIGN.cmd reference_py\\out\\_design\\target.fa"
        )
    if os.path.exists(value):
        parts = []
        with open(value, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    if parts:
                        break
                    continue
                parts.append(line.strip())
        return "".join(parts)
    return value


def parse_region(value: str | None, length: int) -> tuple[int, int] | None:
    if not value:
        return None
    a, b = value.split(":")
    return (int(a), min(int(b), length))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, help="FASTA yolu ya da dizi")
    ap.add_argument("--out", default="out/design")
    ap.add_argument("--left-region", help="a:b (0-tabanli, dislayici)")
    ap.add_argument("--right-region", help="a:b")
    ap.add_argument("--product", default="70:200", help="min:max bp")
    ap.add_argument("--tm", default="58:62", help="min:max C")
    ap.add_argument("--tm-opt", type=float, default=60.0)
    ap.add_argument("--length", default="18:27", help="min:max nt")
    ap.add_argument("--gc", default="40:60", help="min:max yuzde (orn. 40:60)")
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    ap.add_argument("--require-th03", action="store_true",
                    help="Ta'da baglanma kesri TH-03 FAIL esiginin altindaki primerleri ele")
    ap.add_argument("--max-pairs", type=int, default=2000)
    ap.add_argument("--top", type=int, default=20, help="panelde kac cift")
    ap.add_argument("--profile", choices=sorted(profiles.PROFILES),
                    help="uygulama profili; kisitlar buradan gelir, "
                         "komut satirinda acikca verilen bayraklar yine onceliklidir")

    # Profil YALNIZCA varsayilanlari degistirir; kullanici bir bayragi acikca
    # verdiyse o kazanir.
    pre, _ = ap.parse_known_args()
    prof = None
    if pre.profile:
        prof = profiles.get(pre.profile)
        c = prof.constraints
        ap.set_defaults(
            product=f"{c.product_min}:{c.product_max}",
            tm=f"{c.tm_min}:{c.tm_max}", tm_opt=c.tm_opt,
            length=f"{c.length_min}:{c.length_max}",
            gc=f"{c.gc_min * 100:g}:{c.gc_max * 100:g}",
        )
    args = ap.parse_args()

    tmin, tmax = (float(x) for x in args.tm.split(":"))
    qc_profile = prof.qc if prof is not None else Profile()
    theta_min = qc_profile.theta_fail if args.require_th03 else 0.0
    target = nn.validate(read_target(args.target))
    os.makedirs(args.out, exist_ok=True)
    if args.require_th03:
        print(f"TH-03 filtresi     : theta >= {theta_min:g} (Ta {args.ta:g} C, {args.primer_nm:g} nM)")

    if prof is not None:
        print(f"profil: {prof.name} — {prof.description}")
        for note in prof.notes:
            print(f"  not     : {note}")
        for pc in prof.pending_checks:
            print(f"  BEKLEMEDE: {pc}")
        with open(os.path.join(args.out, "profile.tsv"), "w",
                  encoding="utf-8", newline="\n") as fh:
            fh.write("alan\tdeger\n")
            settings = prof.as_dict()
            settings["tm_window_c"] = [tmin, tmax]
            if args.require_th03:
                settings["design_theta_min"] = theta_min
            for k, v in settings.items():
                fh.write(f"{k}\t{v}\n")

    pmin, pmax = (int(x) for x in args.product.split(":"))
    lmin, lmax = (int(x) for x in args.length.split(":"))
    gmin, gmax = (float(x) / 100.0 for x in args.gc.split(":"))

    cons = dsg.DesignConstraints(
        length_min=lmin, length_max=lmax,
        tm_min=tmin, tm_max=tmax, tm_opt=args.tm_opt,
        gc_min=gmin, gc_max=gmax,
        product_min=pmin, product_max=pmax,
        max_pairs=args.max_pairs,
        theta_min=theta_min,
    )
    cond = salt.Conditions(
        monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
        dntp_M=args.dntp / 1e3, primer_M=args.primer_nm / 1e9,
        anneal_c=args.ta,
    )

    result = dsg.design(
        target, cons, cond,
        parse_region(args.left_region, len(target)),
        parse_region(args.right_region, len(target)),
    )

    pairs_path = os.path.join(args.out, "design_pairs.tsv")
    cand_path = os.path.join(args.out, "design_candidates.tsv")
    rej_path = os.path.join(args.out, "design_rejections.tsv")
    panel_path = os.path.join(args.out, "design_panel.tsv")

    with open(pairs_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["left_primer", "right_primer", "product_start",
                    "product_end", "product_length", "tm_left_c", "tm_right_c",
                    "delta_tm_c", "heterodimer_dG37", "heterodimer3p_dG37",
                    "uncalibrated_score", "pareto_rank"])
        for p in result["pairs"]:
            w.writerow(p.as_row())

    with open(cand_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["orientation", "sequence", "start", "end", "length",
                    "gc_pct", "tm_c", "theta_at_ta", "selfdimer_dG37",
                    "selfdimer3p_dG37"])
        for c in result["left_candidates"] + result["right_candidates"]:
            w.writerow(c.as_row())

    with open(rej_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["orientation", "sequence", "start", "rule", "reason",
                    "value"])
        for r in result["rejections"]:
            w.writerow(r.as_row())

    locus = read_locus(args.target)
    with open(panel_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if locus:
            w.writerow(["pair_id", "primer1", "primer2", "left_start",
                        "right_start", "expected_amplicon_length"])
        else:
            w.writerow(["pair_id", "primer1", "primer2"])
        for i, p in enumerate(result["pairs"][: args.top]):
            row = [f"DESIGN_{i:03d}", p.left.sequence, p.right.sequence]
            if locus:
                offset = locus[1]
                row += [offset + p.product_start, offset + p.product_end,
                        p.product_length]
            w.writerow(row)

    extra_paths = []
    if prof is not None and prof.name == "ctdna" and result["pairs"]:
        # CT-02 / CT-03: her aday amplikon icin butun-fragman ve Poisson
        # tespit olasiligi. cfDNA tepe uzunlugu 166 bp (Snyder ve ark. 2016).
        det_path = os.path.join(args.out, "ctdna_detectability.tsv")
        with open(det_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["pair_id", "amplicon_bp", "p_intact_fragment",
                        "genome_equivalents", "vaf", "lambda", "p_detect"])
            for i, p in enumerate(result["pairs"][: args.top]):
                for vaf in (0.0001, 0.001, 0.01):
                    d = profiles.detectability(p.product_length, vaf)
                    w.writerow([f"DESIGN_{i:03d}", d["amplicon_bp"],
                                d["p_intact_fragment"], d["genome_equivalents"],
                                d["vaf"], d["lambda"], d["p_detect"]])
        extra_paths.append(det_path)
        d0 = profiles.detectability(result["pairs"][0].product_length, 0.001)
        print(f"CT-02 butun fragman : {d0['p_intact_fragment']:.3f} "
              f"(amplikon {d0['amplicon_bp']} bp, cfDNA 166 bp)")
        print(f"CT-03 VAF %0.1'de   : lambda {d0['lambda']:.2f}, "
              f"P(tespit) {d0['p_detect']:.3f} (10 ng cfDNA)")

    import collections
    by_rule = collections.Counter(r.rule for r in result["rejections"])
    print(f"hedef uzunlugu      : {result['target_length']} nt")
    print(f"aday primer (F/R)   : {len(result['left_candidates'])} / "
          f"{len(result['right_candidates'])}")
    print(f"cift adayi          : {len(result['pairs'])}")
    if result["pairs"]:
        best = result["pairs"][0]
        print(f"en iyi cift         : {best.product_length} bp, "
              f"dTm {best.delta_tm:.2f} C, skor {best.score:.4f} "
              f"(Pareto {best.pareto_rank})")
    print("elenme gerekceleri  : " +
          ", ".join(f"{k}={v}" for k, v in by_rule.most_common()))
    print("skor kalibre degildir; ozgulluk taramasi yapilmadi")
    print("ciktilar:")
    for p in [pairs_path, cand_path, rej_path, panel_path] + extra_paths:
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
