#!/usr/bin/env python3
"""SPC-08 icin birincil 24 kromozom DISINDAKI GRCh38.p14 dizilerini ayirir.

Girdi : downloads/GCF_000001405.40_GRCh38.p14_genomic.fna.gz (RefSeq, tam derleme)
Cikti : OligoScout/data/annotations/grch38_extra_contigs.fa      (dizi)
        OligoScout/data/annotations/grch38_extra_contigs.tsv     (ad, rol, uzunluk, aciklama)

Roller: mitochondrion (NC_012920.1), alt_locus, fix_patch, novel_patch,
unlocalized, unplaced. Birincil NC_000001..NC_000024 atlanir (onlar PPFM
indeksinde zaten taranir).

    python3 cli/build_extra_contigs.py [--gcf YOL] [--out-dir YOL]
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(HERE)
DEFAULT_GCF = os.path.join(ROOT, "downloads", "GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
DEFAULT_OUT = os.path.join(ROOT, "OligoScout", "data", "annotations")
PRIMARY = re.compile(r"^NC_0000(0[1-9]|1[0-9]|2[0-4])\.")


def role(header: str) -> str:
    h = header.lower()
    if "mitochondrion" in h:
        return "mitochondrion"
    if "alternate locus" in h:
        return "alt_locus"
    if "patch" in h:
        return "novel_patch" if "novel" in h else "fix_patch"
    if "unlocalized" in h:
        return "unlocalized"
    if "unplaced" in h:
        return "unplaced"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gcf", default=DEFAULT_GCF)
    ap.add_argument("--out-dir", default=DEFAULT_OUT)
    a = ap.parse_args()
    if not os.path.isfile(a.gcf):
        print(f"RefSeq derlemesi yok: {a.gcf}", file=sys.stderr)
        return 1
    os.makedirs(a.out_dir, exist_ok=True)
    fa = os.path.join(a.out_dir, "grch38_extra_contigs.fa")
    tsv = os.path.join(a.out_dir, "grch38_extra_contigs.tsv")
    rows, keep, name, length, desc = [], False, None, 0, ""
    with gzip.open(a.gcf, "rt") as src, open(fa, "w", newline="\n") as out:
        def flush():
            if keep and name:
                rows.append((name, role(desc), length, desc))
        for line in src:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                desc = line[1:].strip()[len(name):].strip()
                keep = not PRIMARY.match(name)
                length = 0
                if keep:
                    out.write(f">{name}\n")
                continue
            if keep:
                out.write(line)
                length += len(line.strip())
        flush()
    with open(tsv, "w", newline="\n") as fh:
        fh.write("contig\trole\tlength\tdescription\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    by = {}
    for _, r, n, _ in rows:
        by.setdefault(r, [0, 0])
        by[r][0] += 1
        by[r][1] += n
    for r, (c, n) in sorted(by.items()):
        print(f"{r:14s} {c:4d} dizi {n:>12,d} bp")
    print(f"yazildi: {fa}\n         {tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
