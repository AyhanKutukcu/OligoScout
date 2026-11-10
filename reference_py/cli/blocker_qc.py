#!/usr/bin/env python3
"""CT-05 (WT blocker / ddPCR prob ayrimi) ve ilgili dimer kontrolleri.

Allel-ayrim tahlili reaktiflerini alir (mutant prob, WT prob; istege bagli
ileri/geri primer ve blocker) ve rapor uretir:

  CT-05   mutant vs WT prob Tm farki (kanal ayrimi)
  DIM-05  prob x primer dimerleri (verilirse)
  blocker x prob/primer dimerleri (verilirse)

Diziler MODIFIYE notasyonda olabilir: +baz = LNA, /uc/ = 5'/3' isareti.
LNA stabilizasyonu MODELLENMEZ; Tm'ler DNA iskeletinin alt sinirlaridir.

Kullanim:
    python3 cli/blocker_qc.py --mutant CCC+CT+C+T+C+TGAC \\
        --wt C+CT+C+T+ATG+A+CC \\
        --forward AGTCCTTTCTGTGTCTTCCCA --reverse TCCAGCATCTCCAGCAGC \\
        --blocker AGCATCTCCAGCAGCAGG+T+CA+T+AGAGGG/3Phos/ \\
        --out out/blocker

Cikti: ct05_qc.tsv, ct05_qc_records.jsonl
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import oligomod, params, salt, seqpair  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutant", required=True, help="mutant prob (modifiye)")
    ap.add_argument("--wt", required=True, help="WT prob (modifiye)")
    ap.add_argument("--forward", default="", help="ileri primer")
    ap.add_argument("--reverse", default="", help="geri primer")
    ap.add_argument("--blocker", default="", help="WT blocker (modifiye)")
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "out", "blocker"))
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--salt-model", default="owczarzy_2008",
                    choices=list(salt.MODELS))
    args = ap.parse_args()

    cond = salt.Conditions(monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
                           dntp_M=args.dntp / 1e3, primer_M=args.primer_nm / 1e9,
                           anneal_c=args.ta)
    os.makedirs(args.out, exist_ok=True)

    checks: list[dict] = []
    checks.append({**seqpair.check_ct05(args.mutant, args.wt, cond=cond,
                                        salt_model=args.salt_model),
                   "role": "probe_pair"})

    mut = oligomod.parse(args.mutant)
    wt = oligomod.parse(args.wt)
    fwd = oligomod.parse(args.forward).dna if args.forward else ""
    rev = oligomod.parse(args.reverse).dna if args.reverse else ""
    blk = oligomod.parse(args.blocker).dna if args.blocker else ""

    # DIM-05: her prob x primer (varsa). DNA iskeleti uzerinde.
    if fwd and rev:
        for label, probe in (("mutant", mut.dna), ("wt", wt.dna)):
            checks.append({**seqpair.check_primer_probe_dimers(fwd, rev, probe),
                           "role": f"DIM-05_{label}"})
        if blk:
            checks.append({**seqpair.check_primer_probe_dimers(fwd, rev, blk),
                           "role": "DIM-05_blocker"})

    digest = params.registry_digest()
    records = [{**c, "conditions": cond.as_dict(), "param_set_sha256": digest,
                "engine_version": "oligoscout-thermo-ref/blocker"}
               for c in checks]

    tsv_path = os.path.join(args.out, "ct05_qc.tsv")
    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["check", "role", "severity", "value", "threshold", "note"])
        for c in checks:
            w.writerow([c["id"], c.get("role", ""), c["severity"],
                        json.dumps(c["value"], ensure_ascii=False,
                                   separators=(",", ":")),
                        json.dumps(c["threshold"], ensure_ascii=False,
                                   separators=(",", ":")),
                        c.get("note", "")])
    rec_path = os.path.join(args.out, "ct05_qc_records.jsonl")
    with open(rec_path, "w", encoding="utf-8", newline="\n") as fh:
        for r in records:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")

    for c in checks:
        print(f"{c['id']:7s} {c.get('role',''):16s} {c['severity']:5s} "
              + json.dumps(c["value"], ensure_ascii=False))
    print("LNA stabilizasyonu modellenMEDI; Tm'ler DNA iskeleti alt siniri")
    print("esikler kalibre degildir (REP-03)")
    print("ciktilar:")
    for p in (tsv_path, rec_path):
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
