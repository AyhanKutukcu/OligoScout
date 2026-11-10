#!/usr/bin/env python3
"""Üretim dimer modelinde tek bazlık bulge etkisini karşılaştır."""
from __future__ import annotations

import argparse
import csv
from itertools import combinations_with_replacement
from pathlib import Path
import sys

REF = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REF))
from oligoscout_thermo import align, profiles
from thermo_report import read_panel


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def energy(result):
    return result.dG37_kcal if result.found else None


def decision(result, warn, fail):
    value = energy(result)
    return "FAIL" if value is not None and value < fail else (
        "WARN" if value is not None and value < warn else "INFO")


def symmetric(a, b, mode, enabled, forward):
    reverse = align.best_dimer(b, a, mode, allow_bulge=enabled)
    return forward.found == reverse.found and (
        not forward.found or abs(forward.dG37_kcal - reverse.dG37_kcal) <= 1e-9)


def compare(name, path):
    labels = {}
    for label, seq in read_panel(str(path))[0]:
        labels.setdefault(align.validate(seq), set()).add(label)
    if not labels:
        raise ValueError(f"Boş panel: {path}")
    qc = profiles.get(name).qc
    for mode in align.MODES:
        warn, fail = ((qc.dimer_dg_warn, qc.dimer_dg_fail) if mode == "any"
                      else (qc.dimer3p_dg_warn, qc.dimer3p_dg_fail))
        for a, b in combinations_with_replacement(sorted(labels), 2):
            old = align.best_dimer(a, b, mode, allow_bulge=False)
            new = align.best_dimer(a, b, mode)
            start, end = energy(old), energy(new)
            delta = end - start if start is not None and end is not None else None
            old_sym = symmetric(a, b, mode, False, old)
            new_sym = symmetric(a, b, mode, True, new)
            check = {("any", True): "DIM-01", ("3prime", True): "DIM-02",
                     ("any", False): "DIM-03", ("3prime", False): "DIM-04"}[mode, a == b]
            yield dict(profile=name, mode=mode, primer_a="/".join(sorted(labels[a])),
                       primer_b="/".join(sorted(labels[b])), self_dimer=a == b, check=check,
                       dG37_off=start, dG37_on=end, delta_delta_G37=delta, bulge_only=start is None and end is not None,
                       severity_off=decision(old, warn, fail), severity_on=decision(new, warn, fail),
                       dim08_off="INFO" if old_sym else "FAIL", dim08_on="INFO" if new_sym else "FAIL")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ctdna", type=Path, default=REF / "out/design_ctdna/design_panel.tsv")
    parser.add_argument("--qpcr-sybr", type=Path, default=REF / "out/design_qpcr_sybr/design_panel.tsv")
    parser.add_argument("--out", type=Path, default=REF / "out/bulge")
    args = parser.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    args.out.mkdir(parents=True, exist_ok=True)
    all_rows, summaries, changes, dim06 = [], [], [], []
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for i, (profile, path) in enumerate((("ctdna", args.ctdna), ("qpcr_sybr", args.qpcr_sybr))):
        rows = list(compare(profile, path))
        all_rows.extend(rows)
        for j, mode in enumerate(align.MODES):
            group = [row for row in rows if row["mode"] == mode]
            changed = set()
            for k, row in enumerate(group):
                for check, before, after in ((row["check"], row["severity_off"], row["severity_on"]),
                                              ("DIM-08", row["dim08_off"], row["dim08_on"])):
                    if before != after:
                        changed.add(k)
                        changes.append({**{key: row[key] for key in ("profile", "mode", "primer_a", "primer_b",
                                                                     "dG37_off", "dG37_on")},
                                        "check": check, "before": before, "after": after})
            deltas = [row["delta_delta_G37"] for row in group if row["delta_delta_G37"] is not None]
            summaries.append(dict(profile=profile, mode=mode, combinations=len(group),
                                  lower_dg37=sum(value < -1e-9 for value in deltas),
                                  bulge_only=sum(row["bulge_only"] for row in group),
                                  changed_dim_decisions=len(changed)))
            ax = axes[i, j]
            ax.hist(deltas, bins=18, color="#32768a")
            ax.set(title=f"{profile} · {mode}", xlabel="ΔΔG₃₇ (kcal/mol)",
                   ylabel="Kombinasyon sayısı (log)", yscale="log")
            ax.text(0.02, 0.95, f"İki modelde puanlanan: {len(deltas)}/{len(group)}",
                    transform=ax.transAxes, va="top", fontsize=8)
            if mode == "3prime":
                for key in ("dG37_off", "dG37_on"):
                    worst = min((row for row in group if not row["self_dimer"] and row[key] is not None),
                                key=lambda row: (row[key], row["primer_a"], row["primer_b"]), default=None)
                    if worst:
                        dim06.append(dict(profile=profile, model=key, primer_a=worst["primer_a"],
                                          primer_b=worst["primer_b"], dG37=worst[key], severity="INFO"))
    write_tsv(args.out / "bulge_karsilastirma.tsv", summaries, summaries[0].keys())
    write_tsv(args.out / "bulge_detay.tsv", all_rows, all_rows[0].keys())
    write_tsv(args.out / "dim_karar_degisiklikleri.tsv", changes,
              ("profile", "mode", "primer_a", "primer_b", "dG37_off", "dG37_on", "check", "before", "after"))
    write_tsv(args.out / "dim06_en_guclu.tsv", dim06,
              ("profile", "model", "primer_a", "primer_b", "dG37", "severity"))
    fig.savefig(args.out / "bulge_etkisi.png", dpi=160)
    plt.close(fig)
    print(f"{len(all_rows)} profil/mod/kombinasyon satırı, {len(changes)} DIM karar değişikliği.")
    print("DIM-06 rapor kararı INFO. DIM-05 ve DIM-07 için mevcut Python panelinde karar yok.")


if __name__ == "__main__":
    main()
