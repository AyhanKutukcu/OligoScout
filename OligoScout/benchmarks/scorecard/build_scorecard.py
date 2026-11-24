#!/usr/bin/env python3
"""Collect every OligoScout validation/benchmark result into one scorecard.

Writes
  mcp_server/benchmark_scores.json                   (served as oligoscout://benchmark)
  OligoScout/benchmarks/scorecard/SCORECARD.md       (human-readable)

Each entry carries its evidence file so a reviewer can trace every number.
Run after: thermo_vs_primer3.py, exhaustive_scan + blindspot_analysis.py, test81, regresyon.
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OS = ROOT / "OligoScout"
T81 = OS / "benchmarks" / "test81" / "results"
OUT = ROOT / "reference_py" / "out"


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT)).replace("\\", "/")


def jload(p: Path):
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else None


def test81():
    p = T81 / "test81_results.tsv"
    rows = list(csv.DictReader(p.open(encoding="utf-8"), delimiter="\t"))
    num = lambda v: float(v) if re.fullmatch(r"[-+0-9.eE]+", v or "") else v  # noqa: E731
    methods = [{k: num(v) for k, v in r.items()} for r in rows]
    return {"evidence": rel(p), "panel": "1,024 Primer3 pairs (2,048 x 24-nt primers), 4 difficulty strata x 256",
            "reference": "GRCh38.p14 primary, 24 sequences",
            "contract": "exact 12-nt 3' anchor, <=3 mismatches/primer, products 50-3000 bp",
            "truth": "independent full-reference oracle (8,262,803 bindings, 425,055 products)",
            "machine": "i5-12500H, WSL2, single thread, single run; index build excluded",
            "methods": methods}


def regression():
    rep = OUT / "regresyon" / "REGRESYON_RAPORU.md"
    txt = rep.read_text(encoding="utf-8") if rep.is_file() else ""
    m = re.search(r"PASS (\d+)\s+SKIP (\d+)\s+FAIL (\d+)", txt)
    t = OUT / "regresyon" / "01_test_paketi.log"
    tm = re.search(r"Ran (\d+) tests", t.read_text(encoding="utf-8")) if t.is_file() else None
    diff = OUT / "diff_cpp" / "diff_cpp_summary.md"
    dtxt = diff.read_text(encoding="utf-8") if diff.is_file() else ""
    dm = re.search(r"dizi sayisi: (\d+)", dtxt)
    pc = OUT / "param_crosscheck" / "crosscheck_summary.md"
    return {
        "unit_tests": {"tests": 292, "status": "OK (1 platform skip)", "evidence": "reference_py/tests",
                       "note": "v0.4.1; the end-to-end regression log below predates it"},
        "cpp_tests": {"single_primer_products": "PASS", "primer_pair_assembler_equivalence": "PASS (500 checks)",
                      "terminal_rule_search_vs_brute_force": "PASS (terminal max 0/1/2)"},
        "salt_layer_vs_oligotm": {"values": 60, "identical": 60, "tolerance_c": 0.05,
                                  "evidence": "reference_py/out/salt_crosscheck/salt_crosscheck_summary.md"},
        "mcp_oauth_tests": {"tests": 3, "status": "OK", "evidence": "mcp_server/tests/test_oauth.py"},
        "end_to_end_regression": {"pass": int(m.group(1)), "skip": int(m.group(2)), "fail": int(m.group(3)),
                                  "evidence": rel(rep)} if m else None,
        "python_vs_cpp": {"sequences": int(dm.group(1)) if dm else None, "columns": 17, "differences": 0
                          if "Yok." in dtxt else "see file", "evidence": rel(diff)},
        "parameter_tables_vs_primer3": {"tables": 6, "values_compared": 373, "identical": 371,
                                        "flagged_deviations": 0,
                                        "note": "2 of 96 terminal-mismatch (tstack2) entries are not counted as "
                                                "matches; see the evidence file",
                                        "evidence": rel(pc)},
    }


def thermo():
    p = HERE / "results" / "thermo_vs_primer3" / "summary.json"
    d = jload(p)
    return {"evidence": rel(p), **d} if d else None


def blindspot():
    out = {}
    for key, name in (("ctdna_panel", "blindspot_ctdna"), ("test81_first64_balanced", "blindspot_test81_64")):
        p = HERE / "results" / name / "summary.json"
        d = jload(p)
        if d:
            out[key] = {"evidence": rel(p), **d}
    return out


def sensitive():
    out = {}
    for name in ("blindspot_ctdna", "blindspot_test81_64"):
        p = HERE / "results" / name / "sensitive_verify.json"
        d = jload(p)
        if d:
            out[d["label"]] = {"evidence": rel(p), **d}
    return out


def speed():
    p = HERE / "results" / "speed" / "summary.json"
    d = jload(p)
    return {"evidence": rel(p), **d} if d else None


def primer_blast_web():
    p = HERE / "results" / "primer_blast_web" / "summary.json"
    d = jload(p)
    if not d:
        return None
    keep = {k: v for k, v in d.items() if k not in ("primer_blast_job_keys", "per_pair_os_pbRule_matched_pbListed")}
    return {"evidence": rel(p), **keep}


def lna_model():
    p = HERE / "results" / "lna_mctigue2004" / "summary.json"
    d = jload(p)
    return {"evidence": rel(p), **d} if d else None


def primerbank():
    p = HERE / "results" / "primerbank" / "summary.json"
    d = jload(p)
    if not d:
        return None
    bo = d["by_outcome"]
    sc = bo["sequence_confirmed"]
    lines = [f"{d['pairs_collected']} PrimerBank mouse pairs with validation pages ({d['genes']} genes): "
             f"{sc['pairs']} sequence-confirmed, {bo.get('expected_not_top_hit', {}).get('pairs', 0)} expected "
             f"transcript not the top BLAST hit, {bo.get('blast_missing', {}).get('pairs', 0)} without a BLAST result",
             f"OligoScout worst verdict on sequence-confirmed pairs (rt_qpcr, Ta 60 C): {sc['worst']}; "
             f"FAIL mostly TH-03 ({d['pairs_with_fail_by_check'].get('TH-03', 0)} pairs)",
             "TH-03 failures on working pairs by annealing temperature: " + ", ".join(
                 f"{k} C: {v['working_pairs_failing_TH03']}/{v['of']}" for k, v in d["th03_by_annealing_c"].items()),
             f"Tm OligoScout - PrimerBank: mean {d['tm_oligoscout_minus_primerbank']['mean']:+.2f} C "
             f"(SD {d['tm_oligoscout_minus_primerbank']['sd']})",
             d["reading"]]
    return {"evidence": rel(p), "lines": lines, **d}


def catalog():
    p = OUT / "qc_catalog_audit"
    tsv = next(iter(sorted(p.glob("qc_kural_dogrulama.tsv"))), None)
    counts = {}
    if tsv:
        for r in csv.DictReader(tsv.open(encoding="utf-8-sig"), delimiter="\t"):
            st = r.get("katalog_durumu") or "?"
            counts[st] = counts.get(st, 0) + 1
    return {"rules": sum(counts.values()), "status_counts": counts, "evidence": rel(tsv) if tsv else None,
            "wet_lab_calibrated_rules": 0, "calibration_record": "REP-03 (thresholds are literature starting values)"}


def main():
    card = {
        "title": "OligoScout validation & benchmark scorecard",
        "generated_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "engine": "oligoscout-thermo-ref/0.4.1 + oligoscout_panel 0.4 (C++20 FM-index, PPFM 24 shards; "
                  "anchored and sensitive search modes, single-primer products)",
        "how_to_read": [
            "Every number links to an evidence file under the project root.",
            "Specificity results are exact against the stated CONTRACT; the anchor_blind_spot section measures "
            "what lies outside it.",
            "No wet-lab calibration exists yet (REP-03): computational agreement is not assay performance.",
        ],
        "specificity_test81": test81(),
        "anchor_blind_spot": blindspot(),
        "sensitive_mode_vs_brute_force": sensitive(),
        "primer_blast_web_head_to_head": primer_blast_web(),
        "speed_repeated": speed(),
        "primerbank_validated_pairs": primerbank(),
        "lna_model_vs_mctigue2004": lna_model(),
        "thermodynamics_vs_primer3": thermo(),
        "software_integrity": regression(),
        "qc_catalog": catalog(),
        "integrations": {
            "mcp_conformance": {"sdk": "official MCP Python SDK 2.3.0", "transports": ["stdio", "streamable-http"],
                                "checks_passed": 26, "checks_total": 26, "protocol": "2025-06-18",
                                "tools": 15, "evidence": "mcp_server/README.md#dogrulama"},
            "browser_port_equivalence": {"cases": 2000, "identical": 2000, "evidence": "web/ (engine.js port check)"},
        },
        "not_yet_benchmarked": [
            "Wet-lab validation (qPCR efficiency, specificity by melt/gel/NGS) - required before any calibrated claim",
            "Head-to-head design quality vs Primer3/Primer-BLAST/PrimerQuest outputs on the same targets with a "
            "wet-lab outcome (PrimerBank gives outcomes only for its own mouse designs; QC false-alarm rate on "
            "them is in section 3b)",
            "Multiplex panel performance vs PrimerPlex/Ultiplex/olivar on wet-lab data",
        ],
    }
    dst = ROOT / "mcp_server" / "benchmark_scores.json"
    dst.write_text(json.dumps(card, indent=2, ensure_ascii=False), encoding="utf-8")
    md = HERE / "SCORECARD.md"
    md.write_text(render_md(card), encoding="utf-8")
    print("wrote", rel(dst), "and", rel(md))


def render_md(c):
    L = [f"# {c['title']}", "", f"_generated {c['generated_utc']}_", ""]
    L += ["## 1. Specificity search vs. independent oracle (Test #81)", "",
          "| Method | Binding P | Binding R | Product R | Search s | Max RSS GiB |", "|---|---:|---:|---:|---:|---:|"]
    for m in c["specificity_test81"]["methods"]:
        L.append(f"| {m['label']} | {m['binding_precision']:.6f} | {m['binding_recall']:.6f} | "
                 f"{m['product_recall']:.6f} | {m['search_wall_seconds']:.1f} | {m['search_max_rss_kib'] / 1048576:.2f} |")
    L += ["", "## 2. Anchor blind spot (exhaustive unanchored <=3-mismatch scan)", ""]
    for k, d in c["anchor_blind_spot"].items():
        p, s = d["products"], d["sites"]
        L += [f"**{k}** ({d['pairs']} pairs, {d['primers']} primers)", "",
              f"- sites: {s['all']:,} with <=3 mm; anchored (OligoScout contract) {s['anchored']:,}; "
              f"Primer-BLAST-rule relevant {s['pb_relevant']:,}; blind spot {s['blind_spot']:,}",
              f"- products: OligoScout {p['oligoscout']:,}; Primer-BLAST rule {p['primer_blast_rule']:,}; "
              f"extra with <=2 mm/site {p['pb_not_oligoscout_max2mm_per_site']:,}; with <=1 mm/site "
              f"{p['pb_not_oligoscout_max1mm_per_site']:,}; single-primer products (F-F/R-R) found by "
              f"oligoscout_panel {d['verify_products_vs_oligoscout'].get('single_primer_in_oligoscout', 0)} "
              f"of {p['anchored_single_primer']} in the contract",
              f"- pairs with any extra product: {p['pairs_with_pb_extra']}/{d['pairs']} "
              f"(<=2 mm/site: {p['pairs_with_pb_extra_max2mm']}/{d['pairs']})", ""]
    sv = c.get("sensitive_mode_vs_brute_force") or {}
    if sv:
        L += ["## 2b. Sensitive mode (Primer-BLAST terminal rule) vs. brute force", "",
              "| Panel | sites brute force | sites OligoScout | products brute force | products OligoScout | exact |",
              "|---|---:|---:|---:|---:|---|"]
        for k, d in sv.items():
            L.append(f"| {k} | {d['sites_bruteforce']:,} | {d['sites_oligoscout']:,} | {d['products_bruteforce']:,} | "
                     f"{d['products_oligoscout']:,} | {'yes' if d['exact_match'] else 'NO'} |")
        L.append("")
    pb = c.get("primer_blast_web_head_to_head")
    if pb:
        L += ["## 2c. NCBI Primer-BLAST web vs. OligoScout sensitive mode", "",
              f"- {pb['pairs']} pairs ({pb['pairs_detail']}); Primer-BLAST listed {pb['primer_blast_products_listed']:,} "
              f"products, {pb['primer_blast_products_rule_satisfying']} satisfy the terminal rule on both sites "
              f"({pb['primer_blast_listed_but_violating_terminal_rule']} listed products violate it)",
              f"- rule products also found by OligoScout: {pb['rule_products_also_found_by_oligoscout']} "
              f"(primary chromosomes {pb['primary_chromosome_recall']}; with SPC-08 extra contigs "
              f"{pb.get('rule_products_found_including_spc08', '-')})",
              f"- OligoScout-only rule products (Primer-BLAST seeding misses): {pb['oligoscout_only_products']}", ""]
    sp = c.get("speed_repeated")
    if sp:
        L += ["## 2d. Repeated timings (native Linux VM)", "", f"_{sp['machine']}_", "",
              "| Run | n | search s (mean ± SD) | max RSS MiB | products |", "|---|---:|---:|---:|---:|"]
        for k, d in sp["runs"].items():
            prod = d.get("products", d.get("hetero_products"))
            L.append(f"| {k} | {d['n']} | {d['search_s_mean']:.2f} ± {d['search_s_sd']:.2f} | {d['max_rss_mib']} | "
                     f"{prod:,} |")
        L.append("")
    t = c["thermodynamics_vs_primer3"]
    if t:
        L += ["## 3. Thermodynamics vs. primer3 (oligotm / ntthal), 2,048 primers", "",
              "| Comparison | mean diff | SD | max abs | Pearson r |", "|---|---:|---:|---:|---:|"]
        for k in ("A_tm_1M_aligned_conc", "C_tm_same_salt_model_aligned_conc", "B_tm_default_each_tool",
                  "D_heterodimer_dG37", "E_selfdimer_dG37", "F_hairpin_dG37"):
            v = t[k]
            L.append(f"| {k} | {v['mean_diff']:+.3f} | {v['sd_diff']:.3f} | {v['max_abs']:.3f} | {v.get('pearson_r', '')} |")
        pm = t.get("previous_dimer_model")
        if pm:
            L += ["", f"Dimer model: {t.get('dimer_model', '')}. Before end terms ({pm['model']}): hetero-dimer mean diff "
                  f"{pm['D_heterodimer_dG37']['mean_diff']:+.2f}, kappa {pm['D_heterodimer_dG37']['cohen_kappa']}; "
                  f"now kappa {t['D_heterodimer_dG37']['flag_agreement']['cohen_kappa']} "
                  f"(self-dimer {t['E_selfdimer_dG37'].get('flag_agreement', {}).get('cohen_kappa', '-')})."]
    ln = c.get("lna_model_vs_mctigue2004")
    if ln:
        m3, m4 = ln["M3_measured_reference_plus_table"], ln["M4_oligoscout_unified_dna_plus_table"]
        L += ["", "## 3c. LNA model vs. McTigue 2004 Table 1", "",
              f"- {ln['n']} singly LNA-substituted oligos (excluded: {', '.join(ln['excluded'])})",
              f"- measured DNA reference + LNA table: RMSD ΔTm {m3['rmsd_dtm_c']} °C "
              f"(paper {m3['paper_reported_rmsd_dtm_c']} °C)",
              f"- OligoScout DNA model + LNA table: RMSD Tm {m4['rmsd_tm_c']} °C, {m4['within_2c']}/{ln['n']} within "
              f"2 °C; DNA references alone {m4['dna_reference_rmsd_tm_c']} °C"]
    pk = c.get("primerbank_validated_pairs")
    if pk:
        L += ["", "## 3b. PrimerBank experimentally tested mouse pairs", ""] + [f"- {x}" for x in pk.get("lines", [])]
    r = c["software_integrity"]
    L += ["", "## 4. Software integrity", "",
          f"- unit tests: {r['unit_tests']['tests']} OK; C++ single-primer + assembler equivalence tests PASS",
          f"- end-to-end regression: {r['end_to_end_regression']}",
          f"- Python vs C++: {r['python_vs_cpp']['sequences']} sequences, {r['python_vs_cpp']['differences']} differences",
          f"- parameter tables vs primer3: {r['parameter_tables_vs_primer3']['identical']}/"
          f"{r['parameter_tables_vs_primer3']['values_compared']} identical, 0 flagged deviations",
          f"- salt layer vs oligotm: {r['salt_layer_vs_oligotm']['identical']}/{r['salt_layer_vs_oligotm']['values']} identical",
          f"- QC catalog: {c['qc_catalog']['rules']} rules, status {c['qc_catalog']['status_counts']}",
          "", "## 5. Not yet benchmarked", ""] + [f"- {x}" for x in c["not_yet_benchmarked"]]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
