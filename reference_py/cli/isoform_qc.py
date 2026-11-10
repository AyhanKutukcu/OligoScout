#!/usr/bin/env python3
"""RT-03 izoform kapsami: primer ciftlerinin hedefledigi transkript izoformlari.

Girdi : --pairs TSV (pair_id, gene, chrom, fwd_start, fwd_end, rev_start, rev_end;
        koordinatlar GENOM, 0-tabanli yari-acik)
        --exons transcript_exons.tsv (cli/build_annotations.py ciktisi; GENCODE)
Cikti : isoform_qc.tsv            cift basina RT-03 (hedeflenen/hedeflenmeyen)
        isoform_qc_records.jsonl  REP-02 kanit kaydi

Transkriptom FM-index GEREKMEZ; GENCODE ekson anotasyonu yeterlidir. Ekson
dosyasi yalniz panelde gecen genler icin bellege okunur (tam dosya taranir ama
yalniz ilgili genler tutulur)."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import isoform, params  # noqa: E402


def read_tsv(path: str) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_exons_for_genes(path: str, genes: set[str]) -> isoform.TranscriptIndex:
    """transcript_exons.tsv'yi tarar, yalniz `genes` icin ekson tutar."""
    idx = isoform.TranscriptIndex()
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["gene_name"]
            if gene not in genes:
                continue
            tid = row["transcript_id"]
            g = idx.by_gene.setdefault(gene, {})
            chrom, exons = g.get(tid, (row["chrom"], []))
            exons.append((int(row["start"]), int(row["end"])))
            g[tid] = (chrom, exons)
    return idx


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True, help="cift koordinat TSV'si")
    ap.add_argument("--exons", required=True, help="transcript_exons.tsv")
    ap.add_argument("--out", default="out/isoform")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    pairs = read_tsv(args.pairs)
    genes = {r["gene"] for r in pairs if r.get("gene")}
    index = load_exons_for_genes(args.exons, genes)
    digest = params.registry_digest()

    qc_path = os.path.join(args.out, "isoform_qc.tsv")
    rec_path = os.path.join(args.out, "isoform_qc_records.jsonl")
    import collections
    verdicts = collections.Counter()
    with open(qc_path, "w", newline="", encoding="utf-8") as fq, \
         open(rec_path, "w", encoding="utf-8") as fj:
        qw = csv.writer(fq, delimiter="\t", lineterminator="\n")
        qw.writerow(["pair_id", "gene", "n_isoforms", "n_targeted",
                     "targeted", "untargeted", "severity"])
        for r in pairs:
            rec = isoform.check_rt03(
                r["gene"], r["chrom"],
                (int(r["fwd_start"]), int(r["fwd_end"])),
                (int(r["rev_start"]), int(r["rev_end"])), index)
            v = rec["value"]
            verdicts[rec["severity"]] += 1
            qw.writerow([r.get("pair_id", ""), v["gene"],
                         v.get("n_isoforms", 0), v.get("n_targeted", 0),
                         ";".join(v.get("targeted", [])),
                         ";".join(v.get("untargeted", [])), rec["severity"]])
            fj.write(json.dumps({"pair_id": r.get("pair_id", ""), "check": rec,
                                 "param_set_sha256": digest},
                                ensure_ascii=False, sort_keys=True) + "\n")

    print(f"cift sayisi      : {len(pairs)}")
    print(f"gen sayisi       : {len(genes)}  (ekson indeksinde: {index.genes()})")
    print(f"RT-03 kararlari  : " +
          ", ".join(f"{k}={v}" for k, v in sorted(verdicts.items())))
    print(f"ciktilar         : {qc_path}")
    print(f"                   {rec_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
