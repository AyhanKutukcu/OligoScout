#!/usr/bin/env python3
"""Iki uygulama profilinin ozgulluk/multipleks karsilastirma sekli ve tablosu.

Makaledeki iddia: cift duzeyinde olculen ozgulluk, ayni primerler bir arada
calistirildiginda ortaya cikan davranisi ONGORMUYOR.

Girdi : her profil icin RUN_PANEL cikti klasoru (panel_qc.tsv, final_ranked.tsv,
        multiplex/multiplex_summary.tsv, multiplex/cross_products.tsv)
Cikti : profile_comparison.png + profile_comparison.tsv

DIKKAT: Bu betik matplotlib gerektirir ve referans katmanin "yalniz standart
kutuphane" kuralinin DISINDADIR; WSL boru hattina (run_panel.sh) bagli degildir.
Sekil uretimi icin ayri calistirilir:

    python3 cli/compare_profiles.py \
        --panel qpcr_sybr=out/panel_design_qpcr_sybr \
        --panel ctdna=out/panel_design_ctdna \
        --label qpcr_sybr="qPCR profile\n(70-150 bp)" \
        --out out
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

BLUE, VERM, GREY = "#0072B2", "#D55E00", "#9E9E9E"
CYCLE = [BLUE, VERM, "#009E73", "#CC79A7"]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def collect(folder):
    qc = read_tsv(os.path.join(folder, "panel_qc.tsv"))
    fin = read_tsv(os.path.join(folder, "final_ranked.tsv"))
    mpx = {r["check"]: r for r in
           read_tsv(os.path.join(folder, "multiplex", "multiplex_summary.tsv"))}
    cross = [r for r in read_tsv(os.path.join(folder, "multiplex",
                                              "cross_products.tsv"))
             if r["kind"] != "designed"]
    note = mpx["MPX-02"]["note"]
    designed = int(note.split("tasarlanmis ")[1].split(",")[0])
    exact = sum(1 for r in cross if int(r["total_mismatches"]) == 0)
    cross_load = designed_load = None
    if "MPX-02f" in mpx:
        cross_load, designed_load = (float(x) for x in
                                     mpx["MPX-02f"]["value"].split(" / "))
    verd = collections.Counter(r["verdict"] for r in fin)
    return {
        "per_pair_offtarget": [int(r["offtarget_total"]) for r in qc],
        "offtarget_total": sum(int(r["offtarget_total"]) for r in qc),
        "offtarget_exact": sum(int(r["offtarget_exact"]) for r in qc),
        "pairs": len(qc),
        "pairs_with_offtarget": sum(1 for r in qc
                                    if int(r["offtarget_total"]) > 0),
        "binding_sites": sum(int(r["sites_f"]) + int(r["sites_r"]) for r in qc),
        "verdict": verd,
        "designed": designed, "cross": len(cross), "cross_exact": exact,
        "cross_mismatched": len(cross) - exact,
        "cross_combinations": note.split("capraz kombinasyon ")[1],
        "mpx01": mpx["MPX-01"]["value"],
        "designed_load": designed_load, "cross_load": cross_load,
        "load_ratio": cross_load / designed_load if designed_load else None,
        "worst_dimer": mpx["DIM-06"]["value"].split()[0],
    }


def figure(names, d, labels, colors, path):
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8,
                         "axes.labelsize": 8, "xtick.labelsize": 6,
                         "ytick.labelsize": 6, "legend.fontsize": 6,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 300, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.0))
    rng = np.random.default_rng(0)

    ax = axes[0]
    for i, n in enumerate(names):
        v = np.array(d[n]["per_pair_offtarget"], dtype=float)
        ax.scatter(i + rng.uniform(-0.13, 0.13, v.size), v, s=16,
                   facecolor=colors[n], edgecolor="white", linewidth=0.4, zorder=3)
        ax.hlines(np.median(v), i - 0.22, i + 0.22, color=colors[n],
                  linewidth=2, zorder=4)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_yticks([0, 1, 10, 100]); ax.set_yticklabels(["0", "1", "10", "100"])
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([labels[n] for n in names])
    ax.set_xlim(-0.5, len(names) - 0.4); ax.set_ylim(-0.3, 900)
    ax.set_ylabel("Off-target products\nper pair")
    ax.set_title("Each pair screened alone", loc="left")
    for i, n in enumerate(names):
        ax.annotate(f"{d[n]['offtarget_total']} in total",
                    xy=(i + 0.06, max(1, max(d[n]['per_pair_offtarget']))),
                    xytext=(i + 0.40, 230), fontsize=7, color=colors[n],
                    ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", color=colors[n], lw=0.6))

    ax = axes[1]
    for i, n in enumerate(names):
        c = d[n]
        ax.bar(i, c["designed"], 0.55, color=GREY, edgecolor="white", linewidth=0.5)
        ax.bar(i, c["cross_exact"], 0.55, bottom=c["designed"], color=colors[n],
               edgecolor="white", linewidth=0.5)
        ax.bar(i, c["cross_mismatched"], 0.55,
               bottom=c["designed"] + c["cross_exact"], color=colors[n],
               alpha=0.38, edgecolor="white", linewidth=0.5)
        ax.text(i, c["designed"] + c["cross"] + 10, f"{c['cross']} cross",
                ha="center", fontsize=7, color=colors[n])
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([labels[n] for n in names])
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylim(0, 1.7 * max(d[n]["designed"] + d[n]["cross"] for n in names))
    ax.set_ylabel("Genome-wide products")
    ax.set_title("All primer combinations", loc="left")
    ax.legend(handles=[mpl.patches.Patch(facecolor=GREY, label="designed pairs"),
                       mpl.patches.Patch(facecolor="#555555", label="cross, perfect match"),
                       mpl.patches.Patch(facecolor="#555555", alpha=0.38,
                                         label="cross, with mismatches")],
              loc="upper center", handlelength=1.1, borderpad=0.2)

    ax = axes[2]
    measured = []
    for i, n in enumerate(names):
        c = d[n]
        if c["designed_load"] is None:
            ax.text(i, 0.04, "not measured", transform=ax.get_xaxis_transform(),
                    ha="center", fontsize=7, color=colors[n])
            continue
        measured.extend((c["designed_load"], c["cross_load"]))
        ax.bar(i - 0.17, c["designed_load"], 0.3, color=GREY,
               edgecolor="white", linewidth=0.5)
        ax.bar(i + 0.17, c["cross_load"], 0.3, color=colors[n],
               edgecolor="white", linewidth=0.5)
        if c["load_ratio"] is not None:
            ax.annotate(f"{c['load_ratio']:.2f}",
                        (i, max(c["designed_load"], c["cross_load"])),
                        xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=8, fontweight="bold", color=colors[n])
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([labels[n] for n in names])
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylim(0, 1.75 * max([1.0, *measured]))
    ax.set_ylabel("cfDNA-weighted load")
    ax.set_title("166 bp fragment weighting", loc="left")
    ax.legend(handles=[mpl.patches.Patch(facecolor=GREY, label="designed"),
                       mpl.patches.Patch(facecolor="#555555",
                                         label="perfect-match cross")],
              loc="upper center", handlelength=1.1, borderpad=0.2,
              title="bold number = cross / designed", title_fontsize=6)

    for ax, letter in zip(axes, "abc"):
        ax.text(-0.26, 1.16, letter, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top", ha="left")
    fig.suptitle("Pair-level specificity does not predict multiplex behaviour",
                 x=0.005, y=0.99, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path)
    return path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", action="append", required=True,
                    metavar="ad=klasor", help="en az iki kez verin")
    ap.add_argument("--label", action="append", default=[],
                    metavar="ad=etiket", help="sekildeki eksen etiketi")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    names, folders = [], {}
    for item in args.panel:
        name, _, folder = item.partition("=")
        names.append(name); folders[name] = folder
    if len(names) < 2:
        raise SystemExit("en az iki --panel gerekir")
    labels = {n: n for n in names}
    for item in args.label:
        name, _, text = item.partition("=")
        labels[name] = text.replace("\\n", "\n")
    colors = {n: CYCLE[i % len(CYCLE)] for i, n in enumerate(names)}

    d = {n: collect(folders[n]) for n in names}
    os.makedirs(args.out, exist_ok=True)
    png = figure(names, d, labels, colors,
                 os.path.join(args.out, "profile_comparison.png"))

    cols = ["profile", "pairs", "binding_sites_genomewide",
            "offtarget_products_total", "offtarget_products_exact",
            "pairs_with_offtarget", "verdict_PASS", "verdict_WARN",
            "verdict_FAIL", "multiplex_designed_products",
            "multiplex_cross_products", "multiplex_cross_exact",
            "multiplex_cross_mismatched", "cross_combinations",
            "mpx01_tm_range_c", "mpx02f_designed_load", "mpx02f_cross_load",
            "mpx02f_ratio", "worst_3p_cross_dimer_kcal_mol"]
    tsv = os.path.join(args.out, "profile_comparison.tsv")
    with open(tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for n in names:
            c = d[n]
            w.writerow([n, c["pairs"], c["binding_sites"], c["offtarget_total"],
                        c["offtarget_exact"], c["pairs_with_offtarget"],
                        c["verdict"].get("PASS", 0), c["verdict"].get("WARN", 0),
                        c["verdict"].get("FAIL", 0), c["designed"], c["cross"],
                        c["cross_exact"], c["cross_mismatched"],
                        c["cross_combinations"], c["mpx01"], c["designed_load"],
                        c["cross_load"],
                        round(c["load_ratio"], 3) if c["load_ratio"] is not None else None,
                        c["worst_dimer"]])

    for n in names:
        c = d[n]
        ratio = f"{c['load_ratio']:.2f}" if c["load_ratio"] is not None else "olculemedi"
        print(f"{n:12s} cift disi urun {c['offtarget_total']:>5} "
              f"({c['pairs_with_offtarget']}/{c['pairs']} cift) | "
              f"capraz urun {c['cross']:>4} (tam eslesen {c['cross_exact']}) | "
              f"cfDNA yuk orani {ratio}")
    print("esikler kalibre degildir (REP-03)")
    print("ciktilar:")
    for p in (png, tsv):
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
