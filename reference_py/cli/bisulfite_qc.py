#!/usr/bin/env python3
"""BIS-04 iplik ozgullugu: CT/GA donusmus genom taramasindan cift kararlari.

Girdi : <ct_prefix>.products.tsv ve <ga_prefix>.products.tsv
        (her biri `run_bisulfite_panel.sh` -> oligoscout_panel ciktisi)
Cikti : bisulfite_qc.tsv            cift basina BIS-04 karari
        bisulfite_qc_records.jsonl  REP-02 kanit kaydi

Konvansiyon: bir cift, DONUSTURULMUS genomun yalniz amaclanan donusum ipliginde
(CT ya da GA) bir urun vermeli; diger donusumde urun iplik-caprazdir (risk).
Amaclanan donusum --intended ile ya da panel 'conversion' sutunuyla verilir.

NOT: Bu tarama DONUSTURULMUS genoma karsidir; donusmemis (orijinal) genoma karsi
ozgulluk ayri bir panel aramasidir (standart SPC ailesi).
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import bisulfite, params  # noqa: E402


def read_tsv(path: str) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def count_by_pair(products: list[dict]) -> tuple[dict, dict]:
    total = collections.Counter()
    intended = collections.Counter()
    for p in products:
        pid = p["pair_id"]
        total[pid] += 1
        if int(p.get("total_mismatches", 0)) == 0:
            intended[pid] += 1
    return total, intended


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ct-prefix", required=True,
                    help="CT donusturulmus taramanin cikti oneki (örn. out/bis/CT/panel)")
    ap.add_argument("--ga-prefix", required=True,
                    help="GA donusturulmus taramanin cikti oneki")
    ap.add_argument("--out", default="out/bisulfite")
    ap.add_argument("--intended", default="CT", choices=["CT", "GA"],
                    help="amaclanan donusum ipligi (panel 'conversion' sutunu yoksa)")
    ap.add_argument("--panel", help="opsiyonel: pair_id<TAB>conversion (CT/GA) TSV'si")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    ct = read_tsv(args.ct_prefix + ".products.tsv")
    ga = read_tsv(args.ga_prefix + ".products.tsv")
    ct_total, ct_int = count_by_pair(ct)
    ga_total, ga_int = count_by_pair(ga)

    intended_of: dict[str, str] = {}
    if args.panel:
        for row in read_tsv(args.panel):
            low = {k.lower(): v for k, v in row.items()}
            pid = low.get("pair_id")
            conv = (low.get("conversion") or "").upper()
            if pid and conv in ("CT", "GA"):
                intended_of[pid] = conv

    pair_ids = sorted(set(ct_total) | set(ga_total) | set(ct_int) | set(ga_int),
                      key=lambda x: (len(x), x))
    digest = params.registry_digest()

    qc_path = os.path.join(args.out, "bisulfite_qc.tsv")
    rec_path = os.path.join(args.out, "bisulfite_qc_records.jsonl")
    verdicts = collections.Counter()
    with open(qc_path, "w", newline="", encoding="utf-8") as fq, \
         open(rec_path, "w", encoding="utf-8") as fj:
        qw = csv.writer(fq, delimiter="\t", lineterminator="\n")
        qw.writerow(["pair_id", "intended_conversion", "ct_products",
                     "ga_products", "ct_exact", "ga_exact",
                     "bis04_severity", "note"])
        for pid in pair_ids:
            conv = intended_of.get(pid, args.intended)
            rec = bisulfite.check_bis04(
                ct_int.get(pid, 0), ct_total.get(pid, 0),
                ga_int.get(pid, 0), ga_total.get(pid, 0), conv)
            verdicts[rec["severity"]] += 1
            qw.writerow([pid, conv, ct_total.get(pid, 0), ga_total.get(pid, 0),
                         ct_int.get(pid, 0), ga_int.get(pid, 0),
                         rec["severity"], rec["note"]])
            fj.write(json.dumps({"pair_id": pid, "check": rec,
                                 "param_set_sha256": digest},
                                ensure_ascii=False, sort_keys=True) + "\n")

    print(f"cift sayisi      : {len(pair_ids)}")
    print(f"BIS-04 kararlari : " +
          ", ".join(f"{k}={v}" for k, v in sorted(verdicts.items())))
    print(f"parametre digest : {digest[:16]}...")
    print(f"ciktilar         : {qc_path}")
    print(f"                   {rec_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
