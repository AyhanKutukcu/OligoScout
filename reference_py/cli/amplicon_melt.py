#!/usr/bin/env python3
"""Amplikon erime egrisi ve QPCR-02 (tek erime piki) kontrolu.

Girdi : --design klasoru (design_pairs.tsv + design_panel.tsv) ve tasarimda
        kullanilan --target FASTA; amplikon dizileri buradan kesilir.
        Ek olarak --sequences id<TAB>sequence dosyasi verilebilir (ornegin
        hedef disi urunler icin).
Cikti : amplicon_melt.tsv        cift basina Tm, pik sayisi, QPCR-02 karari
        amplicon_melt_curves.tsv uzun bicimde egriler (T, theta, -dtheta/dT)
        amplicon_melt.png        matplotlib varsa egri + turev grafigi

Model Zimm-Bragg tipi isbirlikli helis-kangaldir; sigma KALIBRE EDILMEMIS bir
model sabitidir ve dongu entropisi henuz yoktur (tablolar `pending`). Karar
pik SAYISI ve pikler arasi ayrim uzerinden verilir, mutlak Tm uzerinden degil.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import amplicon as am  # noqa: E402
from oligoscout_thermo import nn, salt  # noqa: E402


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path):
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
    ap.add_argument("--design", help="design_pairs.tsv + design_panel.tsv klasoru")
    ap.add_argument("--target", help="tasarimda kullanilan FASTA")
    ap.add_argument("--sequences", help="id<TAB>sequence TSV'si (istege bagli)")
    ap.add_argument("--out", help="varsayilan: --design")
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--t-min", type=float, default=50.0)
    ap.add_argument("--t-max", type=float, default=98.0)
    ap.add_argument("--step", type=float, default=0.25)
    ap.add_argument("--min-separation", type=float, default=2.0,
                    help="QPCR-02: bu kadar ayrik iki pik WARN")
    ap.add_argument("--tm-resolution", type=float, default=1.0,
                    help="iki urun Tm'i bu kadar yakinsa ayirt edilemez")
    args = ap.parse_args()

    items: list[tuple[str, str]] = []
    if args.design and args.target:
        target = read_fasta(args.target)
        pairs = read_tsv(os.path.join(args.design, "design_pairs.tsv"))
        panel = read_tsv(os.path.join(args.design, "design_panel.tsv"))
        for row, src in zip(panel, pairs):
            if row["primer1"] != src["left_primer"]:
                raise SystemExit("design_panel ve design_pairs eslesmiyor")
            items.append((row["pair_id"],
                          target[int(src["product_start"]):int(src["product_end"])]))
    if args.sequences:
        for r in read_tsv(args.sequences):
            items.append((r["id"], nn.validate(r["sequence"])))
    if not items:
        raise SystemExit("girdi yok: --design ile --target, ya da --sequences verin")

    out_dir = args.out or args.design or "."
    os.makedirs(out_dir, exist_ok=True)
    cond = salt.Conditions(monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
                           dntp_M=args.dntp / 1e3, anneal_c=args.ta)

    rows, curves, profiles = [], [], {}
    for pid, seq in items:
        p = am.melt_profile(seq, cond=cond, t_min=args.t_min, t_max=args.t_max,
                            step=args.step)
        v = am.qpcr02(p, args.min_separation)
        profiles[pid] = p
        rows.append([pid, len(seq), round(100 * nn.gc_fraction(seq), 2),
                     round(p.tm_c, 3), v["n_peaks"],
                     ";".join(f"{t:.2f}" for t in v["peaks_c"]),
                     v["severity"], v["severity"], v["note"]])
        for t, th, d in zip(p.temperatures, p.helicity, p.derivative):
            curves.append([pid, t, round(th, 6), round(d, 6)])

    # Urunler arasi Tm ayrimi: ayni tupte iki urun ancak Tm'leri ayrilirsa
    # erime egrisinde ayirt edilebilir.
    tms = sorted((r[3], r[0]) for r in rows if r[3] == r[3])
    gaps = [(b[0] - a[0], a[1], b[1]) for a, b in zip(tms, tms[1:])]
    unresolved = [g for g in gaps if g[0] < args.tm_resolution]

    melt_path = os.path.join(out_dir, "amplicon_melt.tsv")
    with open(melt_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        # amp02 == qpcr02: AMP-02 (amplikon erime profili, dtheta/dT pik
        # sayisi) ile QPCR-02 AYNI hesaptir; iki katalog kimligi de ayni
        # karari tasir, farkli sutunlarda acikca raporlanir.
        w.writerow(["id", "length_bp", "gc_pct", "tm_c", "n_peaks", "peaks_c",
                    "qpcr02", "amp02", "note"])
        w.writerows(rows)

    curve_path = os.path.join(out_dir, "amplicon_melt_curves.tsv")
    with open(curve_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "temperature_c", "theta", "minus_dtheta_dt"])
        w.writerows(curves)

    fig_path = None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
        for pid, p in profiles.items():
            axes[0].plot(p.temperatures, p.helicity, linewidth=0.8, alpha=0.75)
            axes[1].plot(p.temperatures, p.derivative, linewidth=0.8, alpha=0.75)
        axes[0].set_xlabel("Temperature (C)"); axes[0].set_ylabel("Paired fraction")
        axes[1].set_xlabel("Temperature (C)"); axes[1].set_ylabel("-d(theta)/dT")
        axes[0].set_title("Melting curves", loc="left")
        axes[1].set_title("Derivative (melt peaks)", loc="left")
        fig.suptitle(f"{len(profiles)} amplicon, uncalibrated cooperativity model",
                     x=0.005, ha="left", fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        fig_path = os.path.join(out_dir, "amplicon_melt.png")
        fig.savefig(fig_path, dpi=300)
    except ImportError:
        pass

    multi = [r for r in rows if r[6] == "WARN"]
    print(f"amplikon sayisi      : {len(rows)}")
    print(f"Tm araligi           : "
          f"{min(r[3] for r in rows):.2f} .. {max(r[3] for r in rows):.2f} C")
    print(f"QPCR-02 WARN         : {len(multi)} urun birden cok pik veriyor"
          + (f" ({', '.join(r[0] for r in multi[:6])})" if multi else ""))
    print(f"Tm ile ayrilamayan    : {len(unresolved)} cift"
          + (f" (ornek: {unresolved[0][1]}/{unresolved[0][2]}, "
             f"{unresolved[0][0]:.2f} C)" if unresolved else ""))
    print("model sabiti sigma kalibre degildir; dongu entropisi yok (REP-03)")
    print("ciktilar:")
    for p_ in (melt_path, curve_path, fig_path):
        if p_:
            print(f"  {p_}  ({os.path.getsize(p_)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
