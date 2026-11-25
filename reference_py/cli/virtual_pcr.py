#!/usr/bin/env python3
"""Sanal PCR (Kalendar et al. 2024 virtualPCR kapsami) — kullanici kaliplari ve GRCh38.

Ornekler:
  # plazmit (dairesel), dejenere primer, 5' kuyruk, amplikon dizileri
  python3 cli/virtual_pcr.py --template pUC19.fa --circular all \\
      --primer M13F=gtaaaacgacggccagtGTAAAACGACGGCCAGT --primer RT2=CARATGGAYGTNAARAC \\
      --mode all --extract --out out/vpcr

  # oligo tablosu (id, sequence, role, pair, tail5, tail3), problar ve baglantili arama
  python3 cli/virtual_pcr.py --template genomes.fa --oligos oligos.tsv --link L,R,0,5000

  # insan genomu (24 kromozom, FM-index; dejenere primer dogal olarak aranir)
  python3 cli/virtual_pcr.py --genome grch38 --primer RTp=CARATGGAYGTNAARAC \\
      --primer RTm=CATRTCRTCNACRTA --mode all --min-product 200 --max-product 500

Varsayilan kural virtualPCR ile ayni: toplam <= 3 uyumsuzluk, son 7 bazda <= 1;
urun 50-5000 bp. Ta = en dusuk primer-yer Tm'i + ln(urun uzunlugu).
Ciktilar: sites.tsv, products.tsv, links.tsv, summary.json.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "reference_py"))

from oligoscout_thermo import salt, vpcr  # noqa: E402

DATA = os.path.join(ROOT, "OligoScout", "data")
PANEL_BINARY = os.environ.get("OLIGOSCOUT_PANEL_BINARY",
                              os.path.join(ROOT, "OligoScout", "build", "release", "oligoscout_panel"))
PPFM = os.path.join(DATA, "indexes", "ppfm")
GENOME_FA = os.path.join(DATA, "genomes", "GRCh38.p14.fa")
GENOME_MANIFEST = os.path.join(DATA, "benchmarks", "reference", "GRCh38.p14.primary24.manifest.tsv")


def read_oligos(args) -> list[vpcr.Oligo]:
    out = []
    if args.oligos:
        with open(args.oligos, encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
                r = {k.lower(): (v or "") for k, v in r.items()}
                out.append(vpcr.parse_oligo(r.get("id") or r.get("oligo_id"), r["sequence"],
                                            r.get("role") or "primer", r.get("pair", ""),
                                            r.get("tail5", ""), r.get("tail3", "")))
    for spec, role in [(p, "primer") for p in args.primer] + [(p, "probe") for p in args.probe]:
        oid, _, seq = spec.partition("=")
        if not seq:
            raise SystemExit(f"--{role} ID=DIZI bekleniyor: {spec}")
        out.append(vpcr.parse_oligo(oid, seq, role))
    ids = [o.id for o in out]
    if len(set(ids)) != len(ids):
        raise SystemExit("oligo kimlikleri tekil olmali")
    return out


def genome_sites(oligos, rule, args):
    """GRCh38 24 kromozom: oligoscout_panel (FM-index), her oligo bir kez."""
    if not (os.path.isfile(PANEL_BINARY) and os.access(PANEL_BINARY, os.X_OK)):
        raise SystemExit(f"oligoscout_panel yok: {PANEL_BINARY}")
    if rule.max_mismatch > 3:
        raise SystemExit("genom aramasinda en fazla 3 uyumsuzluk")
    with tempfile.TemporaryDirectory(prefix="oligoscout-vpcr-g-") as tmp:
        panel = os.path.join(tmp, "panel.tsv")
        with open(panel, "w") as fh:
            fh.write("pair_id\tprimer1\tprimer2\n")
            for o in oligos:
                if len(o.core) > 64:
                    raise SystemExit(f"{o.id}: genom aramasinda en fazla 64 nt")
                fh.write(f"{o.id}\t{o.core}\t{o.core}\n")
        prefix = os.path.join(tmp, "g")
        subprocess.run([PANEL_BINARY, os.path.join(PPFM, "manifest.sa8.tsv"), PPFM, panel, prefix,
                        "--mode", "sensitive", "--max-mismatch", str(rule.max_mismatch),
                        "--terminal-window", str(rule.terminal_window),
                        "--terminal-max", str(rule.terminal_max),
                        "--min-amplicon", "1", "--max-amplicon", "1", "--single-primer", "off",
                        "--cache", "1"], check=True, capture_output=True, text=True)
        sites = []
        with open(prefix + ".bindings.tsv", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["role"] != "F":
                    continue
                mask = r["mask"]
                L = len(mask)
                term = sum(1 for i in range(max(0, L - rule.terminal_window), L) if mask[i] == "X")
                sites.append(vpcr.Site(r["pair_id"], r["chrom"], int(r["start"]), int(r["end_exclusive"]),
                                       r["strand"], int(r["mismatches"]), term, mask,
                                       r["site_plus_strand"]))
    lengths = {}
    with open(GENOME_MANIFEST, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lengths[r["chromosome"]] = (int(r["length"]), False)
    return sites, lengths


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--template", action="append", default=[], help="FASTA (tekrarlanabilir)")
    ap.add_argument("--genome", choices=("grch38",), help="insan genomu, FM-index")
    ap.add_argument("--oligos", help="TSV: id, sequence[, role, pair, tail5, tail3]")
    ap.add_argument("--primer", action="append", default=[], help="ID=DIZI (kucuk harf = kuyruk)")
    ap.add_argument("--probe", action="append", default=[], help="ID=DIZI (yalniz baglanma yerleri)")
    ap.add_argument("--link", action="append", default=[], help="A,B,min_gap,max_gap")
    ap.add_argument("--mode", choices=("pairs", "all"), help="pairs: ayni cift; all: tum primerler")
    ap.add_argument("--circular", default="", help="all ya da kayit adlari (virgulle)")
    ap.add_argument("--max-mismatch", type=int, default=3)
    ap.add_argument("--terminal-window", type=int, default=7)
    ap.add_argument("--terminal-max", type=int, default=1)
    ap.add_argument("--template-n", choices=("reject", "mismatch", "match"), default="reject")
    ap.add_argument("--g-tolerant", type=int, default=None, metavar="K",
                    help="en fazla K adet G.T, G.G, G.A uyumsuzlugu sinira sayilmaz (raporlanir)")
    ap.add_argument("--genome-engine", choices=("auto", "fm", "scan"), default="auto",
                    help="fm: 24 kromozom FM-index; scan: GRCh38 FASTA kaba kuvvet taramasi")
    ap.add_argument("--min-product", type=int, default=50)
    ap.add_argument("--max-product", type=int, default=5000)
    ap.add_argument("--engine", choices=("auto", "cpp", "python"), default="auto")
    ap.add_argument("--extract", action="store_true", help="urun dizilerini yaz")
    ap.add_argument("--no-thermo", action="store_true", help="Tm/Ta hesaplama")
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    ap.add_argument("--ta", type=float, default=60.0, help="kosul kaydi icin; Ta urun basina hesaplanir")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--out", default="out/vpcr")
    args = ap.parse_args()
    if not args.template and not args.genome:
        ap.error("--template ya da --genome gerekli")

    oligos = read_oligos(args)
    if not oligos:
        ap.error("en az bir oligo gerekli (--oligos / --primer / --probe)")
    omap = {o.id: o for o in oligos}
    mode = args.mode or ("pairs" if any(o.pair for o in oligos) else "all")
    rule = vpcr.ScanRule(args.max_mismatch, args.terminal_window, args.terminal_max, args.template_n,
                         args.g_tolerant is not None, args.g_tolerant if args.g_tolerant is not None else 2)
    circular = "all" if args.circular == "all" else [c for c in args.circular.split(",") if c]
    cond = salt.Conditions(monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
                           dntp_M=args.dntp / 1e3, primer_M=args.primer_nm / 1e9, anneal_c=args.ta)

    sites, records, rec_seq, engine_used = [], {}, {}, ""
    if args.genome:
        need_scan = (rule.g_tolerant or rule.max_mismatch > 3 or rule.template_n != "reject"
                     or any(len(o.core) > 64 for o in oligos))
        g_engine = args.genome_engine if args.genome_engine != "auto" else ("scan" if need_scan else "fm")
        if g_engine == "fm":
            if need_scan:
                ap.error("bu kural FM-index ile aranamaz; --genome-engine scan kullanin")
            sites, records = genome_sites(oligos, rule, args)
            engine_used = "oligoscout_panel (FM-index, sensitive rule, IUPAC)"
        else:
            s, rl = vpcr.scan_cpp(GENOME_FA, oligos, rule, (), threads=args.threads, primary_only=True)
            sites += s
            records.update(rl)
            engine_used = "oligoscout_vpcr_scan on GRCh38.p14 primary 24 (exhaustive)"
    for fa in args.template:
        use_cpp = args.engine == "cpp" or (args.engine == "auto" and vpcr.default_scan_binary())
        recs = vpcr.read_fasta(fa, circular)
        for r in recs:
            rec_seq[r.name] = r
        if use_cpp:
            s, rl = vpcr.scan_cpp(fa, oligos, rule, circular, threads=args.threads)
            engine_used = engine_used or "oligoscout_vpcr_scan (exhaustive, bit-parallel)"
        else:
            s = vpcr.scan_python(recs, oligos, rule)
            rl = {r.name: (r.length, r.circular) for r in recs}
            engine_used = engine_used or "vpcr.scan_python (exhaustive)"
        sites += s
        records.update(rl)

    products = vpcr.assemble(sites, omap, records, args.min_product, args.max_product, mode)
    probe_sites = [s for s in sites if omap[s.oligo].role == "probe"]
    genome_ref = None
    if args.genome and args.extract:
        sys.path.insert(0, os.path.join(ROOT, "web"))
        from local_reference import Reference
        genome_ref = Reference(GENOME_FA, GENOME_MANIFEST)

    tm_cache = {}

    def tm_of(site):
        key = (site.oligo, site.record, site.start, site.strand)
        if key not in tm_cache:
            tm_cache[key] = vpcr.site_tm(omap[site.oligo], site, cond)
        return tm_cache[key]

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "sites.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["oligo_id", "role", "record", "start", "end_exclusive", "strand", "mismatches",
                    "terminal_mismatches", "mask", "site_plus_strand", "wraps", "tm_c", "three_prime_paired"])
        for s in sites:
            t = {} if args.no_thermo else tm_of(s)
            w.writerow([s.oligo, omap[s.oligo].role, s.record, s.start, s.end, s.strand, s.mismatches,
                        s.terminal_mismatches, s.mask, s.site, int(s.wraps), t.get("tm_c", ""),
                        t.get("three_prime_paired", "")])
    cols = ["product_id", "record", "topology", "orientation", "left_oligo", "right_oligo", "start",
            "end_exclusive", "wraps", "template_length", "amplicon_length", "left_mismatches",
            "right_mismatches", "tm_left_c", "tm_right_c", "ta_c", "probes_inside"]
    if args.extract:
        cols += ["template_sequence", "amplicon_sequence"]
    with open(os.path.join(args.out, "products.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for i, p in enumerate(products, 1):
            tl = None if args.no_thermo else tm_of(p.left).get("tm_c")
            tr = None if args.no_thermo else tm_of(p.right).get("tm_c")
            row = [f"P{i}", p.record, "circular" if p.circular else "linear", p.orientation,
                   p.left.oligo, p.right.oligo, p.start, p.end, int(p.wraps), p.template_length,
                   p.length, p.left.mismatches, p.right.mismatches, tl if tl is not None else "",
                   tr if tr is not None else "", vpcr.annealing_temperature(tl, tr, p.length) or "",
                   ";".join(vpcr.probes_inside(p, probe_sites))]
            if args.extract:
                if p.record in rec_seq:
                    seqs = vpcr.amplicon_sequence(p, rec_seq[p.record], omap)
                elif genome_ref is not None and p.length <= 5000:
                    tmpl = genome_ref.read(p.record, p.start, p.end).upper()
                    rec = vpcr.Record(p.record, "N" * p.start + tmpl)
                    seqs = vpcr.amplicon_sequence(p, rec, omap)
                else:
                    seqs = {"template_sequence": "", "amplicon_sequence": ""}
                row += [seqs["template_sequence"], seqs["amplicon_sequence"]]
            w.writerow(row)
    links = []
    for spec in args.link:
        a, b, lo, hi = spec.split(",")
        links += vpcr.find_links(sites, a, b, int(lo), int(hi))
    with open(os.path.join(args.out, "links.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["record", "a", "b", "strand", "a_start", "b_start", "gap", "insert_start", "insert_end"])
        for k in links:
            w.writerow([k.record, k.a.oligo, k.b.oligo, k.strand, k.a.start, k.b.start, k.gap,
                        k.insert_start, k.insert_end])
    by_rec = collections.Counter(p.record for p in products)
    summary = {
        "engine": engine_used, "mode": mode,
        "rule": {"max_mismatch": rule.max_mismatch, "terminal_window": rule.terminal_window,
                 "terminal_max": rule.terminal_max, "template_n": rule.template_n},
        "product_range": [args.min_product, args.max_product],
        "records": {k: {"length": v[0], "topology": "circular" if v[1] else "linear"} for k, v in records.items()},
        "oligos": {o.id: {"core": o.core, "tail5": o.tail5, "tail3": o.tail3, "role": o.role,
                          "pair": o.pair, "degeneracy": o.degeneracy,
                          "sites": sum(1 for s in sites if s.oligo == o.id)} for o in oligos},
        "products": len(products), "products_by_record": dict(by_rec),
        "products_wrapping_origin": sum(p.wraps for p in products),
        "links": len(links),
        "ta_formula": "min(Tm_left, Tm_right) + ln(amplicon length) (Kalendar et al. 2024)",
        "conditions": cond.as_dict(),
    }
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    print(f"{len(sites)} yer, {len(products)} urun, {len(links)} baglanti -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
