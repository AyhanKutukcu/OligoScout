#!/usr/bin/env python3
"""Anchor blind-spot analysis from an exhaustive (unanchored) <=3-mismatch scan.

Inputs
  --hits      exhaustive_scan output (primer_id = primer sequence)
  --panel     pairs: pair_id, primer1, primer2 (design_panel.tsv or the Test #81 panel)
  --oligoscout-bindings  optional panel.bindings.tsv from oligoscout_panel (same panel)

Site classes (mask is 5'->3'; the last characters are the 3' end)
  anchored      no mismatch in the 3' 12 nt            -> OligoScout contract (enumerated)
  pb_relevant   <=1 mismatch in the 3' 5 nt            -> Primer-BLAST default rule would still
                                                          treat the site as a potential target
  blind_spot    not anchored but pb_relevant           -> what an anchored search cannot see but
                                                          Primer-BLAST's rule would count

Products (per pair): plus-strand site upstream of a minus-strand site, same chromosome, 50-3000 bp.
  hetero         F->  <-R  and  R->  <-F   (what oligoscout_panel assembles)
  single_primer  F->  <-F  and  R->  <-R   (NOT assembled by oligoscout_panel; Primer-BLAST checks them)
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
from pathlib import Path

ANCHOR, LAST5, LO, HI = 12, 5, 50, 3000


def load_sites(path):
    by_primer = collections.defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = r["mask"]
            site = (r["chrom"], int(r["start"]), int(r["end_exclusive"]), r["strand"], int(r["mismatches"]), m)
            by_primer[r["primer_id"]].append(site)
    return by_primer


def klass(mask):
    anchored = "X" not in mask[-ANCHOR:]
    pb = mask[-LAST5:].count("X") <= 1
    return anchored, pb


def products(sites):
    """sites: list of (chrom, start, end, strand, mm, mask, primer). Returns list of product tuples."""
    out = []
    by_chrom = collections.defaultdict(lambda: ([], []))
    for s in sites:
        by_chrom[s[0]][0 if s[3] == "+" else 1].append(s)
    for chrom, (plus, minus) in by_chrom.items():
        minus.sort(key=lambda s: s[2])
        ends = [s[2] for s in minus]
        import bisect
        for p in plus:
            i = bisect.bisect_left(ends, p[1] + LO)
            while i < len(minus) and minus[i][2] - p[1] <= HI:
                m = minus[i]
                if m[1] >= p[2]:          # same geometry as oligoscout_panel: no overlap
                    out.append((chrom, p[1], m[2], m[2] - p[1], p, m))
                i += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--oligoscout-bindings")
    ap.add_argument("--oligoscout-products")
    ap.add_argument("--label", default="panel")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    sites = load_sites(a.hits)
    pairs = list(csv.DictReader(open(a.panel), delimiter="\t"))

    # ---- site level
    n = collections.Counter()
    mm_pos = collections.Counter()           # 3'-distance of mismatches in blind-spot sites
    for seq, ss in sites.items():
        for s in ss:
            anc, pb = klass(s[5])
            n["all"] += 1; n["anchored"] += anc; n["pb_relevant"] += pb
            n["blind_spot"] += (not anc) and pb
            n["anchored_not_pb"] += anc and not pb
            if not anc and pb:
                L = len(s[5])
                for i, c in enumerate(s[5]):
                    if c == "X" and i >= L - ANCHOR:
                        mm_pos[L - i] += 1
            if not anc and pb and s[5][-1] == "X":
                n["blind_spot_terminal_3p_mismatch"] += 1

    # ---- verification against OligoScout's own enumeration
    verify = None
    if a.oligoscout_bindings:
        os_set = set()
        for r in csv.DictReader(open(a.oligoscout_bindings), delimiter="\t"):
            os_set.add((r["primer"], r["chrom"], int(r["start"]), r["strand"], int(r["mismatches"])))
        ex_set = {(seq, s[0], s[1], s[3], s[4]) for seq, ss in sites.items() for s in ss if klass(s[5])[0]}
        verify = {"oligoscout_unique_sites": len(os_set), "exhaustive_anchored_sites": len(ex_set),
                  "intersection": len(os_set & ex_set), "only_oligoscout": len(os_set - ex_set),
                  "only_exhaustive": len(ex_set - os_set)}
        verify["exact_match"] = verify["only_oligoscout"] == 0 and verify["only_exhaustive"] == 0

    # ---- product level
    rows = []
    tot = collections.Counter()
    for p in pairs:
        ss = []
        for k in ("primer1", "primer2"):
            for s in sites.get(p[k], []):
                ss.append(s + (k,))
        prods = products(ss)
        c = collections.Counter()
        for pr in prods:
            a1, b1 = klass(pr[4][5]); a2, b2 = klass(pr[5][5])
            single = pr[4][6] == pr[5][6] or p["primer1"] == p["primer2"]
            maxmm = max(pr[4][4], pr[5][4])
            c["all_le3mm"] += 1
            c["single_primer_le3mm"] += single
            c["oligoscout"] += a1 and a2 and not single            # oligoscout_panel contract
            c["anchored_single_primer"] += a1 and a2 and single     # anchored but not assembled
            c["primer_blast_rule"] += b1 and b2
            pbx = (b1 and b2) and not (a1 and a2)   # oligoscout_panel >= 0.4 assembles F-F/R-R too
            c["pb_not_oligoscout"] += pbx
            c["pb_not_oligoscout_max1mm_per_site"] += pbx and maxmm <= 1
            c["pb_not_oligoscout_max2mm_per_site"] += pbx and maxmm <= 2
            c["oligoscout_not_pb"] += (a1 and a2) and not (b1 and b2)
        keys = ("all_le3mm", "single_primer_le3mm", "oligoscout", "anchored_single_primer", "primer_blast_rule",
                "pb_not_oligoscout", "pb_not_oligoscout_max1mm_per_site", "pb_not_oligoscout_max2mm_per_site",
                "oligoscout_not_pb")
        row = {"pair_id": p["pair_id"], **{k: c[k] for k in keys}}
        rows.append(row); tot.update(c)
        tot["pairs_with_pb_extra"] += c["pb_not_oligoscout"] > 0
        tot["pairs_with_pb_extra_max2mm"] += c["pb_not_oligoscout_max2mm_per_site"] > 0
    with open(out / "per_pair.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

    prod_verify = None
    if a.oligoscout_products:
        prows = list(csv.DictReader(open(a.oligoscout_products), delimiter="\t"))
        has_single = bool(prows) and "orientation" in prows[0]
        cnt = collections.Counter(r["pair_id"] for r in prows)
        want = (lambda r: r["oligoscout"] + r["anchored_single_primer"]) if has_single else (lambda r: r["oligoscout"])
        mism = [(r["pair_id"], want(r), cnt.get(r["pair_id"], 0)) for r in rows if want(r) != cnt.get(r["pair_id"], 0)]
        prod_verify = {"oligoscout_products": sum(cnt.values()),
                       "exhaustive_anchored_products": tot["oligoscout"] + (tot["anchored_single_primer"] if has_single else 0),
                       "includes_single_primer": has_single,
                       "single_primer_in_oligoscout": sum(1 for r in prows if r.get("orientation") in ("F-F", "R-R")),
                       "pairs_mismatching": mism[:20], "exact_match": not mism}

    summary = {
        "label": a.label, "primers": len(sites), "pairs": len(pairs),
        "contract": {"max_mismatches": 3, "anchor_nt": ANCHOR, "pb_rule": "<=1 mismatch in 3' 5 nt",
                     "product_bp": [LO, HI]},
        "sites": dict(n),
        "blind_spot_fraction_of_pb_relevant": round(n["blind_spot"] / n["pb_relevant"], 4) if n["pb_relevant"] else 0,
        "blind_spot_mismatch_distance_from_3p": dict(sorted(mm_pos.items())),
        "products": {k: tot[k] for k in ("all_le3mm", "single_primer_le3mm", "oligoscout", "anchored_single_primer",
                                         "primer_blast_rule", "pb_not_oligoscout", "pb_not_oligoscout_max1mm_per_site",
                                         "pb_not_oligoscout_max2mm_per_site", "oligoscout_not_pb",
                                         "pairs_with_pb_extra", "pairs_with_pb_extra_max2mm")},
        "verify_sites_vs_oligoscout": verify,
        "verify_products_vs_oligoscout": prod_verify,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
