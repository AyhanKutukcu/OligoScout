#!/usr/bin/env python3
"""Multipleks capraz reaktivite raporu (MPX-01, MPX-02, MPX-03).

Girdi : oligoscout_panel'in <prefix>.bindings.tsv dosyasi + panel TSV'si
Cikti : multiplex_summary.tsv   MPX-01/03 ozetleri
        cross_products.tsv      capraz urunler (tasarlanmamis kombinasyonlar)
        dimer_matrix.tsv        n(n+1)/2 dimer kombinasyonu

Yeni genom taramasi yapilmaz: capraz urunler mevcut baglanma kayitlarindan
kurulur, dolayisiyla aramanin eksiksizlik garantisi korunur.
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import multiplex as mpx  # noqa: E402
from oligoscout_thermo import profiles  # noqa: E402
from oligoscout_thermo import salt  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="oligoscout_panel oneki")
    ap.add_argument("--out", required=True)
    ap.add_argument("--panel", help="panel TSV'si; tasarlanmis kombinasyonlar "
                                    "ve beklenen urun boylari buradan okunur")
    ap.add_argument("--min-amplicon", type=int, default=50)
    ap.add_argument("--max-amplicon", type=int, default=3000)
    ap.add_argument("--max-per-combo", type=int, default=200,
                    help="kombinasyon basina kaydedilecek en fazla urun")
    ap.add_argument("--tm-range-warn", type=float, default=3.0)
    ap.add_argument("--tm-range-fail", type=float, default=5.0)
    ap.add_argument("--ce-gap", type=int, default=10)
    ap.add_argument("--fragment-bp", type=int, default=0,
                    help="cfDNA fragman uzunlugu (orn. 166). Verilirse her "
                         "urun CT-02 butun-fragman olasiligiyla agirliklanir; "
                         "parcalanmis sablonda uzun urunler zaten olusamaz.")
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    cond = salt.Conditions(monovalent_M=args.na / 1e3,
                           magnesium_M=args.mg / 1e3,
                           dntp_M=args.dntp / 1e3,
                           primer_M=args.primer_nm / 1e9)

    sites, owners = mpx.load_bindings(args.prefix + ".bindings.tsv")
    primers: dict[str, str] = {s.primer_id: s.primer for s in sites}

    # Tasarlanmis kombinasyonlar panel dosyasindan, dizi uzerinden kurulur.
    designed: set[frozenset[str]] = set()
    expected_lengths: list[int] = []
    if args.panel:
        with open(args.panel, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                low = {k.lower(): v for k, v in row.items()}
                f = low.get("primer1") or low.get("forward_primer")
                rv = low.get("primer2") or low.get("reverse_primer")
                if f and rv:
                    designed.add(frozenset({f, rv}))
                length = low.get("expected_amplicon_length")
                if length:
                    expected_lengths.append(int(length))
    else:
        # panel verilmediyse bindings'teki pair_id eslesmelerinden tureyelim
        roles: dict[str, dict[str, str]] = collections.defaultdict(dict)
        for seq, labels in owners.items():
            for label in labels:
                pid, role = label.rsplit("_", 1)
                roles[pid][role] = seq
        for pid, rr in roles.items():
            if "F" in rr and "R" in rr:
                designed.add(frozenset({rr["F"], rr["R"]}))

    # Sayim EKSIKSIZ; --max-per-combo yalnizca dosyaya yazilan satirlari kirpar.
    n_total = n_designed = n_cross = n_cross_exact = 0
    load_designed = load_cross = load_cross_exact = 0.0
    combo_counts: collections.Counter = collections.Counter()
    combo_exact: collections.Counter = collections.Counter()
    kept: list = []
    per_combo_kept: collections.Counter = collections.Counter()

    for prod in mpx.iter_products(sites, designed, args.min_amplicon,
                                  args.max_amplicon):
        n_total += 1
        w = (profiles.intact_fragment_probability(prod.length, args.fragment_bp)
             if args.fragment_bp else 0.0)
        if prod.designed_pair:
            n_designed += 1
            load_designed += w
            continue
        n_cross += 1
        load_cross += w
        combo = frozenset({prod.left_id, prod.right_id})
        combo_counts[combo] += 1
        if prod.mismatches == 0:
            n_cross_exact += 1
            load_cross_exact += w
            combo_exact[combo] += 1
        if per_combo_kept[combo] < args.max_per_combo:
            per_combo_kept[combo] += 1
            kept.append(prod)

    def label(seq: str) -> str:
        return "/".join(sorted(owners.get(seq, {seq[:10]})))

    # --- ciktilar -----------------------------------------------------------
    cp_path = os.path.join(args.out, "cross_products.tsv")
    with open(cp_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["left_primer_id", "right_primer_id", "chrom", "start",
                    "end", "length_bp", "total_mismatches", "kind"])
        for prod in sorted(kept, key=lambda x: (x.mismatches, -x.length))[:5000]:
            w.writerow([label(prod.left_id), label(prod.right_id), prod.chrom,
                        prod.start, prod.end, prod.length, prod.mismatches,
                        "cross"])

    cc_path = os.path.join(args.out, "cross_combinations.tsv")
    with open(cc_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["primer_a", "primer_b", "n_products", "n_exact_match"])
        for combo, n in combo_counts.most_common():
            ab = sorted(combo)
            # tek elemanli kombinasyon: primer KENDISIYLE urun veriyor
            # (iki iplikte de baglanip tek primerle amplikon)
            a, b = (ab[0], ab[0]) if len(ab) == 1 else (ab[0], ab[1])
            w.writerow([label(a), label(b), n, combo_exact[combo]])

    dm_path = os.path.join(args.out, "dimer_matrix.tsv")
    dm = mpx.dimer_matrix(primers, cond)
    for row in dm:
        row["primer_a"] = label(row["primer_a"])
        row["primer_b"] = label(row["primer_b"])
    with open(dm_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dm[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(dm)

    tmh = mpx.tm_homogeneity(primers, cond)
    if expected_lengths:
        sep = mpx.size_separation(expected_lengths, args.ce_gap)
        sep_note = "panelde beyan edilen hedef urun boylari"
    else:
        sep = {"n_products": 0, "min_gap_bp": None, "conflicts": 0}
        sep_note = "not_evaluated: --panel ile beklenen boylar verilmedi"

    tm_sev = ("FAIL" if tmh["tm_range_c"] > args.tm_range_fail
              else "WARN" if tmh["tm_range_c"] > args.tm_range_warn else "INFO")
    mpx02_sev = "FAIL" if n_cross else "INFO"
    mpx03_sev = "WARN" if (sep["conflicts"] or 0) else "INFO"

    worst_dimer = min(
        (r for r in dm if r["dG37_3prime"] is not None and not r["self"]),
        key=lambda r: (r["dG37_3prime"], r["primer_a"], r["primer_b"]),
        default=None)   # esitlik bozucu: rapor deterministik olmali

    sm_path = os.path.join(args.out, "multiplex_summary.tsv")
    with open(sm_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["check", "value", "threshold", "severity", "note"])
        w.writerow(["MPX-01", f'{tmh["tm_range_c"]} C',
                    f'warn >{args.tm_range_warn}, fail >{args.tm_range_fail}',
                    tm_sev,
                    f'{label(tmh["tm_min_primer"])} {tmh["tm_min_c"]} C .. '
                    f'{label(tmh["tm_max_primer"])} {tmh["tm_max_c"]} C'])
        w.writerow(["MPX-02", f"{n_cross} capraz urun",
                    "fail: herhangi bir capraz urun", mpx02_sev,
                    f"toplam urun {n_total}, tasarlanmis {n_designed}, "
                    f"tam eslesen capraz {n_cross_exact}, "
                    f"capraz kombinasyon {len(combo_counts)}"])
        w.writerow(["MPX-03", f'min bosluk {sep["min_gap_bp"]} bp',
                    f"warn <{args.ce_gap} bp", mpx03_sev,
                    f'{sep["conflicts"]} catisma, {sep["n_products"]} urun; {sep_note}'])
        if args.fragment_bp:
            ratio = (load_cross_exact / load_designed) if load_designed else None
            w.writerow(["MPX-02f", f"{load_cross_exact:.2f} / "
                                   f"{load_designed:.2f}",
                        f"cfDNA {args.fragment_bp} bp ile agirlikli",
                        "FAIL" if (ratio or 0) >= 0.5 else
                        "WARN" if (ratio or 0) >= 0.1 else "INFO",
                        f"tam eslesen capraz yuk / tasarlanmis yuk = "
                        f"{ratio:.2f}" if ratio is not None else "tasarlanmis yuk 0"])
        if worst_dimer:
            w.writerow(["DIM-06",
                        f'{worst_dimer["dG37_3prime"]} kcal/mol',
                        "en kotu 3'-uclu capraz dimer", "INFO",
                        f'{worst_dimer["primer_a"]} x {worst_dimer["primer_b"]}, '
                        f'{worst_dimer["mismatches_3prime"]} ic uyumsuzluk'])

    print(f"primer sayisi        : {len(primers)}")
    print(f"baglanma yeri        : {len(sites)}")
    print(f"toplam urun          : {n_total}")
    print(f"  tasarlanmis cift   : {n_designed}")
    print(f"  CAPRAZ (MPX-02)    : {n_cross}  -> {mpx02_sev}")
    print(f"    tam eslesen      : {n_cross_exact}")
    n_self = sum(1 for c in combo_counts if len(c) == 1)
    print(f"    kombinasyon      : {len(combo_counts)} / "
          f"{len(primers)*(len(primers)+1)//2} olasi")
    print(f"    tek primer       : {n_self} primer kendisiyle urun veriyor")
    if args.fragment_bp:
        print(f"MPX-02f cfDNA yuku   : tasarlanmis {load_designed:.2f}, "
              f"tam eslesen capraz {load_cross_exact:.2f} "
              f"({args.fragment_bp} bp fragman)")
    print(f"MPX-01 Tm araligi    : {tmh['tm_range_c']} C  -> {tm_sev}")
    print(f"MPX-03 min boy farki : {sep['min_gap_bp']} bp  -> {mpx03_sev}")
    if worst_dimer:
        print(f"en kotu 3'-uclu dimer: {worst_dimer['dG37_3prime']} kcal/mol "
              f"({worst_dimer['primer_a']} x {worst_dimer['primer_b']})")
    print("ciktilar:")
    for p in (sm_path, cp_path, cc_path, dm_path):
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
