#!/usr/bin/env python3
"""Diferansiyel test icin dizi listesi uretir (id, sequence, partner).

Kaynaklar (varsa sirayla eklenir):
  * depo chr22 kiyaslama paneli (32 cift -> 64 primer, her biri esiyle)
  * verilen ek panel TSV'leri (--panel, tekrarlanabilir)
  * thermo_report.DEMO_PANEL kenar durumlari (GC/AT zengin, self-komplementer,
    kisa dizi) -- esi yoktur, yalniz tek-oligo sutunlarini zorlar

Cikti: id<TAB>sequence<TAB>partner (bos partner = yalniz tek-oligo sutunlari)
Liste deterministiktir: ayni girdi ayni dosyayi uretir (REP-01).
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from oligoscout_thermo import nn  # noqa: E402
from thermo_report import DEMO_PANEL, read_panel  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(HERE)
CHR22 = os.path.join(REPO_ROOT, "OligoScout", "data", "benchmarks",
                     "grch38_chr22_panel_32.tsv")


def rows_from_panel(path: str, label: str) -> list[tuple[str, str, str]]:
    _, pairs = read_panel(path)
    out: list[tuple[str, str, str]] = []
    for pair_id, fwd, rev in pairs:
        out.append((f"{label}_{pair_id}_F", fwd, rev))
        out.append((f"{label}_{pair_id}_R", rev, fwd))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--panel", action="append", default=[],
                    help="ek panel TSV'si (tekrarlanabilir)")
    ap.add_argument("--chr22", default=CHR22,
                    help="depo chr22 kiyaslama paneli")
    ap.add_argument("--no-demo", action="store_true",
                    help="kenar durumu dizilerini ekleme")
    args = ap.parse_args()

    rows: list[tuple[str, str, str]] = []
    if args.chr22 and os.path.isfile(args.chr22):
        rows += rows_from_panel(args.chr22, "CHR22")
    for index, path in enumerate(args.panel, start=1):
        if not os.path.isfile(path):
            raise SystemExit(f"panel bulunamadi: {path}")
        rows += rows_from_panel(path, f"P{index}")
    if not args.no_demo:
        rows += [(f"EDGE_{pid}", seq, "") for pid, seq in DEMO_PANEL]

    seen: set[str] = set()
    clean: list[tuple[str, str, str]] = []
    skipped: list[str] = []
    for pid, seq, partner in rows:
        if pid in seen:
            continue
        try:
            seq = nn.validate(seq)
            partner = nn.validate(partner) if partner else ""
        except ValueError as exc:
            skipped.append(f"{pid}: {exc}")
            continue
        seen.add(pid)
        clean.append((pid, seq, partner))

    if not clean:
        raise SystemExit("dizi bulunamadi")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("id\tsequence\tpartner\n")
        for pid, seq, partner in clean:
            fh.write(f"{pid}\t{seq}\t{partner}\n")

    with_partner = sum(1 for _, _, p in clean if p)
    print(f"dizi            : {len(clean)}")
    print(f"esi olan        : {with_partner}")
    if skipped:
        print(f"atlanan         : {len(skipped)}")
        for item in skipped[:5]:
            print(f"  {item}")
    print(f"cikti           : {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
