#!/usr/bin/env python3
"""OligoScout primer/pair QC on PrimerBank pairs with an experimental outcome.

Input : results/primerbank/pairs.json — PrimerBank mouse pairs collected from
        pga.mgh.harvard.edu (search by gene symbol, species Mouse) with the
        validation page's sequencing/BLAST result:
        {gene, acc, id, amp, f, ftm, r, rtm, validated, identity, match_pos, match_len}
Label : 'works' = product sequenced and its BLAST top hit is the expected
        transcript (match position 1, identity >= 90 %).
Output: results/primerbank/{per_pair.tsv, summary.json}

What it measures: how often OligoScout's sequence-intrinsic checks (TH, DIM,
SS, SEQ/CMP, PAIR-01) raise WARN/FAIL on primers that worked in qPCR. It does
NOT measure specificity (mouse genome not indexed) and has no failure class
large enough for a discrimination statistic.
"""
from __future__ import annotations

import collections
import csv
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "reference_py"))
from oligoscout_thermo import engine, profiles, seqpair  # noqa: E402

RES = HERE / "results" / "primerbank"


def main():
    pairs = json.loads((RES / "pairs.json").read_text(encoding="utf-8"))
    prof = profiles.get("rt_qpcr")
    eng = engine.ThermoEngine(profile=prof.qc)
    rows, sev_by_check = [], collections.defaultdict(collections.Counter)
    tm_diff = []
    for p in pairs:
        works = bool(p.get("validated")) and p.get("match_pos") == 1 and (p.get("identity") or 0) >= 90
        outcome = ("sequence_confirmed" if works else
                   "blast_missing" if p.get("match_pos") is None else "expected_not_top_hit")
        checks = []
        tms = []
        for role, s, pb_tm in (("F", p["f"], p["ftm"]), ("R", p["r"], p["rtm"])):
            ev = eng.evaluate(role, s)
            tms.append(ev["thermo"].get("tm_c") if isinstance(ev.get("thermo"), dict) else None)
            checks += [(c["id"], c["severity"]) for c in ev["checks"]]
            checks += [(c["id"], c["severity"]) for c in seqpair.check_primer(s, seqpair.SeqPairConstraints())]
        pe = eng.evaluate_pair(p["id"], p["f"], p["r"])
        checks += [(c["id"], c["severity"]) for c in pe["checks"]]
        worst = max((sv for _, sv in checks), key={"INFO": 0, "PASS": 0, "WARN": 1, "FAIL": 2}.get, default="INFO")
        for cid, sv in checks:
            sev_by_check[cid][sv] += 1
        for t, pb in zip(tms, (p["ftm"], p["rtm"])):
            if t is not None:
                tm_diff.append(t - pb)
        rows.append({"primerbank_id": p["id"], "gene": p["gene"], "validated": int(bool(p.get("validated"))),
                     "works": int(works), "outcome": outcome, "identity": p.get("identity", ""), "amplicon": p["amp"],
                     "f": p["f"], "r": p["r"], "pb_tm_f": p["ftm"], "pb_tm_r": p["rtm"],
                     "os_tm_f": round(tms[0], 2) if tms[0] else "", "os_tm_r": round(tms[1], 2) if tms[1] else "",
                     "worst": worst,
                     "fail_checks": ",".join(sorted({c for c, s in checks if s == "FAIL"})),
                     "warn_checks": ",".join(sorted({c for c, s in checks if s == "WARN"}))})
    with (RES / "per_pair.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    ok = [r for r in rows if r["works"]]
    worst = collections.Counter(r["worst"] for r in ok)
    fail_ids = collections.Counter(c for r in ok for c in r["fail_checks"].split(",") if c)
    warn_ids = collections.Counter(c for r in ok for c in r["warn_checks"].split(",") if c)
    n = len(ok)
    by_outcome = {}
    for oc in ("sequence_confirmed", "expected_not_top_hit", "blast_missing"):
        grp = [r for r in rows if r["outcome"] == oc]
        if grp:
            c = collections.Counter(r["worst"] for r in grp)
            by_outcome[oc] = {"pairs": len(grp), "worst": dict(c),
                              "fail_rate": round(c.get("FAIL", 0) / len(grp), 3),
                              "median_amplicon": statistics.median(r["amplicon"] for r in grp)}
    summ = {
        "source": "PrimerBank (pga.mgh.harvard.edu), mouse, searched by gene symbol; validation pages read 2026-11-22",
        "pairs_collected": len(rows), "genes": len({r["gene"] for r in rows}),
        "with_validation_page": sum(r["validated"] for r in rows), "works_sequence_confirmed": n,
        "profile": "rt_qpcr (default conditions)",
        "by_outcome": by_outcome,
        "worst_verdict_on_working_pairs": dict(worst),
        "fail_rate_on_working_pairs": round(worst.get("FAIL", 0) / n, 3) if n else None,
        "pairs_with_fail_by_check": dict(fail_ids.most_common()),
        "pairs_with_warn_by_check": dict(warn_ids.most_common(12)),
        "tm_oligoscout_minus_primerbank": {"n": len(tm_diff), "mean": round(statistics.mean(tm_diff), 2),
                                           "sd": round(statistics.pstdev(tm_diff), 2)} if tm_diff else None,
        "limits": ["'blast_missing' (no sequencing/BLAST result) is a weak failure proxy: PrimerBank notes short "
                   "products (<75 bp) often fail sequencing; too few for sensitivity/AUC",
                   "specificity not evaluated: OligoScout indexes GRCh38 only",
                   "gene list is a convenience sample of well-studied genes"],
    }
    # TH-03 depends on the annealing temperature: how many working pairs fail it at other Ta
    from oligoscout_thermo import salt
    ta_scan = {}
    for ta in (60.0, 58.0, 56.0):
        e2 = engine.ThermoEngine(cond=salt.Conditions(anneal_c=ta), profile=prof.qc)
        bad = 0
        for p in pairs:
            if not next(r for r in rows if r["primerbank_id"] == p["id"])["works"]:
                continue
            sev = [c["severity"] for s_ in (p["f"], p["r"]) for c in e2.evaluate("x", s_)["checks"] if c["id"] == "TH-03"]
            bad += "FAIL" in sev
        ta_scan[f"{ta:g}"] = {"working_pairs_failing_TH03": bad, "of": n}
    summ["th03_by_annealing_c"] = ta_scan
    summ["reading"] = ("PrimerBank recommends 60 C annealing. At 60 C OligoScout's TH-03 (bound fraction < 0.80) fails "
                       "pairs that were sequence-confirmed, because OligoScout's Tm under 50 mM Na+/1.5 mM Mg2+ sits "
                       "near 60-64 C for these primers; TH-02 (salt-model spread > 2 C) warns on every primer. "
                       "Both thresholds are uncalibrated (REP-03) and this sample argues for relaxing them.")
    (RES / "summary.json").write_text(json.dumps(summ, indent=1), encoding="utf-8")
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
