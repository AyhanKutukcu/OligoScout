#!/usr/bin/env python3
"""TH-02 ve TH-04 için eşik değiştirmeyen karar desteği."""
from __future__ import annotations

import argparse
from pathlib import Path
from statistics import correlation
import sys

REF = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REF))
from oligoscout_thermo import concentration, engine, salt
from thermo_report import read_panel
from bulge_compare import write_tsv


def salt_status(spread, warn, fail):
    return "FAIL" if spread > fail else "WARN" if spread > warn else "INFO"


def terminal_status(energy, band, fail):
    return "FAIL" if energy < fail else "INFO" if band[0] <= energy <= band[1] else "WARN"


def measure(label, path, cond, profile):
    for primer_id, seq in read_panel(str(path))[0]:
        model = concentration.build(seq, cond)
        tm = salt.spread(model.tm_1m_K, model.f_gc, model.duplex.n_bp, cond)
        a, b = tm["owczarzy_2008"], tm["von_ahsen_2001"]
        spread = abs(a - b)
        dg = engine._terminal_dG37(seq, profile.terminal_window_nt)
        yield dict(panel=label, primer_id=primer_id, length_nt=len(seq), gc_pct=round(model.f_gc*100, 4),
                   owczarzy_2008_tm_c=round(a, 6), von_ahsen_2001_tm_c=round(b, 6), spread_c=round(spread, 6),
                   dG37_last5_kcal_mol=round(dg, 6),
                   th02=salt_status(spread, profile.salt_spread_warn_c, profile.salt_spread_fail_c),
                   th02_hypothetical_3c=salt_status(spread, 3.0, profile.salt_spread_fail_c),
                   th04=terminal_status(dg, profile.dg3p_warn_range, profile.dg3p_fail),
                   th04_hypothetical_band=terminal_status(dg, (-10.0, -5.0), profile.dg3p_fail))


def pearson(rows, key):
    x = [r[key] for r in rows]
    y = [r["spread_c"] for r in rows]
    return f"{correlation(x, y):.3f}" if len(set(x)) > 1 and len(set(y)) > 1 else "tanımsız"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chr22", type=Path, default=REF.parent / "OligoScout/data/benchmarks/grch38_chr22_panel_32.tsv")
    parser.add_argument("--latest", type=Path, required=True, help="son ctDNA/chr3 panelinin selected_panel.tsv dosyası")
    parser.add_argument("--out", type=Path, default=REF / "out/kalibrasyon")
    args = parser.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cond, profile = salt.Conditions(), engine.Profile()
    groups = {label: list(measure(label, path, cond, profile)) for label, path in (
        ("Depo chr22", args.chr22), ("Son ctDNA/chr3", args.latest))}
    rows = [row for group in groups.values() for row in group]
    if not all(groups.values()):
        raise ValueError("Her iki panel de en az bir primer içermeli.")
    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "th02_th04_kanit.tsv", rows, rows[0].keys())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.7), constrained_layout=True)
    for (label, group), color in zip(groups.items(), ("#126c8c", "#cc633b")):
        ax1.scatter([r["gc_pct"] for r in group], [r["spread_c"] for r in group],
                    s=[25 + 6*(r["length_nt"]-min(r["length_nt"] for r in rows)) for r in group],
                    label=label, color=color, alpha=0.75)
    ax1.axhline(profile.salt_spread_warn_c, color="#555555", linestyle="--", label="Mevcut WARN eşiği")
    ax1.set(title="TH-02 · Tuz modeli yayılımı", xlabel="GC (%) · işaret boyu primer uzunluğu",
            ylabel="|Tm farkı| (°C)")
    ax1.legend(fontsize=8)
    ax2.axvspan(*profile.dg3p_warn_range, color="#97bca5", alpha=0.45, label="Mevcut INFO bandı")
    ax2.hist([[r["dG37_last5_kcal_mol"] for r in group] for group in groups.values()], bins=14,
             label=list(groups), color=["#126c8c", "#cc633b"], histtype="step", linewidth=2)
    ax2.set(title="TH-04 · 3′ uç stabilitesi", xlabel="Son 5 bazın ΔG₃₇ (kcal/mol)", ylabel="Primer sayısı")
    ax2.legend(fontsize=8)
    fig.savefig(args.out / "esik_kanit.png", dpi=170)
    plt.close(fig)

    lines = ["# TH-02 ve TH-04 karar desteği", "", "REP-03 · KALİBRE DEĞİL. Üretim eşikleri korunuyor.", "",
             "Koşullar 60 °C, 200 nM primer, 50 mM monovalan tuz, 1,5 mM Mg ve 0,6 mM dNTP.",
             "GC ve uzunluk ilişkileri betimleyici Pearson katsayılarıdır. Sabit uzunlukta katsayı tanımsızdır.",
             "", "| Panel | n | TH-02 WARN | TH-04 WARN | Yayılım–GC r | Yayılım–uzunluk r |",
             "|---|---:|---:|---:|---:|---:|"]
    for label, group in groups.items():
        lines.append(f"| {label} | {len(group)} | {sum(r['th02']=='WARN' for r in group)} | "
                     f"{sum(r['th04']=='WARN' for r in group)} | {pearson(group, 'gc_pct')} | {pearson(group, 'length_nt')} |")
    salt_warn = sum(r["th02"] == "WARN" for r in rows)
    salt_remove = sum(r["th02"] == "WARN" and r["th02_hypothetical_3c"] == "INFO" for r in rows)
    terminal_remove = sum(r["th04"] == "WARN" and r["th04_hypothetical_band"] == "INFO" for r in rows)
    lines.extend(["", f"1. Yalnız Owczarzy 2008'i kanonik olarak raporlamak **{salt_warn}** TH-02 WARN kaydını "
                  "kaldırır. Tek modelle modeller arası yayılım tanımsız olur. TH-04 etkisi 0.",
                  f"2. İki modeli koruyup TH-02 WARN eşiğini 2,0 °C'den 3,0 °C'ye çıkarmak **{salt_remove}** "
                  "primeri WARN'dan çıkarır. TH-04 INFO bandını ayrıca −9…−6'dan −10…−5 kcal/mol'e genişletmek "
                  f"**{terminal_remove}** primeri TH-04 WARN'dan çıkarır. FAIL eşikleri sabit tutuldu.",
                  "", "Sayılar primer × panel satırlarıdır. Karar, deneysel veriyle kullanıcıya aittir."])
    (args.out / "esik_karar_notu.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{len(rows)} primer, TH-02 senaryoları {salt_warn}/{salt_remove}, TH-04 bant senaryosu {terminal_remove}.")


if __name__ == "__main__":
    main()
