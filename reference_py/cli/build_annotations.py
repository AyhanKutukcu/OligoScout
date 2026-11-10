#!/usr/bin/env python3
"""Ham anotasyonlari kompakt aralik tablolarina indirger (SPC-06, SPC-07).

NEDEN. GENCODE GTF (~59 MB) ve UCSC genomicSuperDups her panel_qc kosusunda
ayristirmak icin cok buyuktur. Bu betik bir kez calisir ve yalnizca gereken
araliklari cikarir:

  SPC-06  psodogen araliklari : GTF'de feature=gene ve gene_type'inda
          'pseudogene' gecen kayitlar -> pseudogene_intervals.tsv
  SPC-07  segmental duplikasyon : genomicSuperDups -> segdup_intervals.tsv

Cikti sema: chrom<TAB>start<TAB>end<TAB>name (0-tabanli, bitis haric).

Kullanim (depo kokunden):
    python3 reference_py/cli/build_annotations.py
    python3 reference_py/cli/build_annotations.py --annotations <klasor>

Kaynak dosyalar yoksa ilgili cikti atlanir (cikis kodu 0).
"""

from __future__ import annotations

import argparse
import glob
import gzip
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(HERE)
DEFAULT_DIR = os.path.join(REPO_ROOT, "OligoScout", "data", "annotations")


def _open(path: str):
    return (gzip.open(path, "rt", encoding="utf-8", errors="replace")
            if path.endswith(".gz") else open(path, encoding="utf-8"))


def find_one(folder: str, *patterns: str) -> str | None:
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(folder, pat)))
        if hits:
            return hits[0]
    return None


def build_pseudogenes(gtf: str, out: str) -> int:
    """GTF feature=gene + gene_type icinde 'pseudogene' -> (chrom,start,end,name)."""
    n = 0
    with _open(gtf) as fh, open(out, "w", encoding="utf-8", newline="\n") as w:
        w.write("chrom\tstart\tend\tname\n")
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = cols[8]
            if "pseudogene" not in attrs:
                continue
            gtype = _attr(attrs, "gene_type")
            if "pseudogene" not in (gtype or ""):
                continue
            name = _attr(attrs, "gene_name") or _attr(attrs, "gene_id") or "?"
            start = int(cols[3]) - 1          # GTF 1-tabanli -> 0-tabanli
            end = int(cols[4])
            w.write(f"{cols[0]}\t{start}\t{end}\t{name}\n")
            n += 1
    return n


def _attr(attrs: str, key: str) -> str | None:
    marker = key + ' "'
    i = attrs.find(marker)
    if i < 0:
        return None
    i += len(marker)
    j = attrs.find('"', i)
    return attrs[i:j] if j > 0 else None


def build_transcript_exons(gtf: str, out: str) -> int:
    """GTF feature=exon -> (chrom,start,end,gene_name,transcript_id) (RT-03).

    Izoform kapsami icin ekson araliklari. 0-tabanli, bitis haric. Buyuktur
    (~1.3M satir) ama panel_qc/isoform_qc gende filtreleyerek okur."""
    n = 0
    with _open(gtf) as fh, open(out, "w", encoding="utf-8", newline="\n") as w:
        w.write("chrom\tstart\tend\tgene_name\ttranscript_id\n")
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            attrs = cols[8]
            gene = _attr(attrs, "gene_name") or _attr(attrs, "gene_id") or "?"
            tid = _attr(attrs, "transcript_id") or "?"
            start = int(cols[3]) - 1          # GTF 1-tabanli -> 0-tabanli
            end = int(cols[4])
            w.write(f"{cols[0]}\t{start}\t{end}\t{gene}\t{tid}\n")
            n += 1
    return n


def build_segdups(superdups: str, out: str) -> int:
    """genomicSuperDups (bin,chrom,start,end,...) -> (chrom,start,end,name)."""
    n = 0
    with _open(superdups) as fh, open(out, "w", encoding="utf-8", newline="\n") as w:
        w.write("chrom\tstart\tend\tname\n")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            # UCSC bin formati: 0=bin 1=chrom 2=chromStart 3=chromEnd
            try:
                start, end = int(cols[2]), int(cols[3])
            except ValueError:
                continue
            w.write(f"{cols[1]}\t{start}\t{end}\tsegdup\n")
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations", default=DEFAULT_DIR,
                    help="ham + kompakt anotasyon klasoru")
    args = ap.parse_args()
    folder = args.annotations
    if not os.path.isdir(folder):
        print(f"anotasyon klasoru yok: {folder}", file=sys.stderr)
        return 0

    gtf = find_one(folder, "*.annotation.gtf.gz", "*.gtf.gz", "*.gtf")
    sd = find_one(folder, "genomicSuperDups*.txt.gz", "genomicSuperDups*.txt",
                  "*SuperDups*")

    if gtf:
        out = os.path.join(folder, "pseudogene_intervals.tsv")
        n = build_pseudogenes(gtf, out)
        print(f"SPC-06 psodogen: {n} aralik -> {out}")
        out = os.path.join(folder, "transcript_exons.tsv")
        n = build_transcript_exons(gtf, out)
        print(f"RT-03 ekson    : {n} ekson -> {out}")
    else:
        print("GTF bulunamadi; SPC-06 / RT-03 atlandi")
    if sd:
        out = os.path.join(folder, "segdup_intervals.tsv")
        n = build_segdups(sd, out)
        print(f"SPC-07 segdup  : {n} aralik -> {out}")
    else:
        print("genomicSuperDups bulunamadi; SPC-07 atlandi")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
