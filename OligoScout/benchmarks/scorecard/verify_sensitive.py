#!/usr/bin/env python3
"""oligoscout_panel --mode sensitive  vs  brute-force exhaustive scan.

Every site with <=3 mismatches and <=1 mismatch in the 3'-terminal 5 nt (the
Primer-BLAST default rule) must be enumerated, with the same mismatch mask, and
the per-pair product counts must equal the brute-force Primer-BLAST-rule count.

    python3 verify_sensitive.py --hits hits.tsv --bindings panel.bindings.tsv \
        --products panel.products.tsv --per-pair results/<x>/per_pair.tsv \
        --contract panel.contract.tsv --label <x> --out results/<x>/sensitive_verify.json
"""
import argparse
import collections
import csv
import json

import blindspot_analysis as B


def main():
    ap = argparse.ArgumentParser()
    for k in ("hits", "bindings", "products", "per_pair", "contract", "label", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args()
    sites = B.load_sites(a.hits)
    ex = {(seq, s[0], s[1], s[3], s[4], s[5]) for seq, ss in sites.items() for s in ss if B.klass(s[5])[1]}
    got = {(r["primer"], r["chrom"], int(r["start"]), r["strand"], int(r["mismatches"]), r["mask"])
           for r in csv.DictReader(open(a.bindings), delimiter="\t")}
    want = {r["pair_id"]: int(r["primer_blast_rule"]) for r in csv.DictReader(open(a.per_pair), delimiter="\t")}
    cnt = collections.Counter(r["pair_id"] for r in csv.DictReader(open(a.products), delimiter="\t"))
    contract = {r["field"]: r["value"] for r in csv.DictReader(open(a.contract), delimiter="\t")}
    out = {"label": a.label, "search_mode": contract.get("search_mode"),
           "search_seconds": float(contract.get("search_seconds", "nan")),
           "sites_bruteforce": len(ex), "sites_oligoscout": len(got),
           "only_bruteforce": len(ex - got), "only_oligoscout": len(got - ex),
           "products_bruteforce": sum(want.values()), "products_oligoscout": sum(cnt.values()),
           "pairs_mismatching": [(k, v, cnt.get(k, 0)) for k, v in want.items() if v != cnt.get(k, 0)][:20]}
    out["exact_match"] = out["only_bruteforce"] == out["only_oligoscout"] == 0 and not out["pairs_mismatching"]
    json.dump(out, open(a.out, "w"), indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
