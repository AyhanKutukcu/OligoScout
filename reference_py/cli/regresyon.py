#!/usr/bin/env python3
"""Tekrarlanabilir uctan uca regresyon: test paketi + butun CLI zinciri.

NEDEN. Depodaki CLI'ler tek tek calisiyor ama hepsini ayni girdilerle
siraya dizen, ciktilarinin semasini denetleyen ve tek rapor ureten
commit edilmis bir kosum yoktu. Bir degisiklikten sonra "her sey hala
calisiyor mu" sorusunun tek komutluk cevabi bu betiktir.

Her adim uc sonuctan birini alir:
  PASS   komut 0 dondu ve bekledigi dosyalari bekledigi basliklarla yazdi
  SKIP   girdisi yok (ornegin C++ tarama ciktisi ya da son panel dosyasi)
  FAIL   komut sifirdan farkli dondu ya da cikti/sema eksik

Bir adimin FAIL olmasi digerlerini durdurmaz. Cikis kodu: herhangi bir
FAIL varsa 1.

Kullanim:
    python3 cli/regresyon.py [--out out/regresyon] [--timeout 900]
                             [--latest <selected_panel.tsv>] [--hizli]

`--hizli` uzun suren tasarim adimlarini atlar (yalnizca denetim adimlari).
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # reference_py
REPO_ROOT = os.path.dirname(HERE)
PY = sys.executable


def rel(*parts: str) -> str:
    return os.path.join(HERE, *parts)


def header_of(path: str) -> list[str]:
    with open(path, encoding="utf-8-sig", newline="") as fh:
        line = fh.readline().rstrip("\r\n")
    return line.split("\t")


def build_steps(out_root: str, args) -> list[dict]:
    target = os.path.join(REPO_ROOT, "hedef_chr22.fa")
    design_out = os.path.join(out_root, "design_ctdna")
    panel_ctdna = rel("out", "panel_design_ctdna")
    panel_qpcr = rel("out", "panel_design_qpcr_sybr")
    prefix_ctdna = os.path.join(panel_ctdna, "panel")
    panel_qc_out = os.path.join(out_root, "panel_qc")
    dump = rel("out", "diff_cpp", "cpp_dump.tsv")
    seq_list = os.path.join(out_root, "diff_cpp", "diff_cpp_sequences.tsv")

    steps: list[dict] = [
        {"name": "01_test_paketi", "argv": [PY, "-m", "unittest", "discover",
                                            "-s", "tests"]},
        {"name": "02_param_transcribe", "argv": [PY, "cli/param_transcribe.py", "--help"],
         "note": "duman testi: alt komut girdisi olmadan calismaz"},
        {"name": "03_param_crosscheck",
         "argv": [PY, "cli/param_crosscheck.py", "--out",
                  os.path.join(out_root, "param_crosscheck")],
         "produces": [(os.path.join(out_root, "param_crosscheck",
                                    "crosscheck.tsv"),
                       ["table_id", "oracle", "key", "quantity"])]},
        {"name": "04_salt_crosscheck",
         "argv": [PY, "cli/salt_crosscheck.py", "--out",
                  os.path.join(out_root, "salt_crosscheck")],
         "ok_codes": (0, 2),
         "note": "cikis 2 = oracle calismadi (WSL disi); FAIL sayilmaz"},
        {"name": "05_terminal_audit",
         "argv": [PY, "cli/terminal_table_audit.py", "--out",
                  os.path.join(out_root, "terminal_audit")],
         "produces": [(os.path.join(out_root, "terminal_audit",
                                    "terminal_audit.tsv"),
                       ["key", "primer_3p", "allele_3p"])]},
        {"name": "06_thermo_report",
         "argv": [PY, "cli/thermo_report.py", "--out",
                  os.path.join(out_root, "thermo")],
         "produces": [(os.path.join(out_root, "thermo", "thermo_report.tsv"),
                       ["primer_id", "sequence", "length_nt"]),
                      (os.path.join(out_root, "thermo", "param_registry.tsv"),
                       ["table_id"])]},
        {"name": "07_bulge_compare",
         "argv": [PY, "cli/bulge_compare.py", "--out",
                  os.path.join(out_root, "bulge")]},
    ]

    if not args.hizli:
        steps += [
            {"name": "08_design_primers",
             "argv": [PY, "cli/design_primers.py", "--target", target,
                      "--profile", "ctdna", "--max-pairs", "50",
                      "--out", design_out],
             "requires": [target],
             "produces": [(os.path.join(design_out, "design_pairs.tsv"),
                           ["left_primer", "right_primer", "product_start",
                            "product_end"]),
                          (os.path.join(design_out, "design_panel.tsv"),
                           ["pair_id", "primer1", "primer2"])]},
            {"name": "08b_sequence_qc",
         "argv": [PY, "cli/sequence_qc.py", "--design", design_out,
                  "--target", target, "--profile", "ctdna", "--gene", "TP53",
                  "--contigs", os.path.join(REPO_ROOT, "OligoScout", "data",
                                            "annotations", "chrM.fa"),
                  "--out", os.path.join(out_root, "seqpair")],
         "requires": [os.path.join(design_out, "design_pairs.tsv"), target],
         "produces": [(os.path.join(out_root, "seqpair", "seqpair_qc.tsv"),
                       ["pair_id", "role", "check", "severity"]),
                      (os.path.join(out_root, "seqpair", "seqpair_panel.tsv"),
                       ["pair_id", "verdict"])]},
        {"name": "09_probe_qc",
             "argv": [PY, "cli/probe_qc.py", "--design", design_out,
                      "--target", target, "--profile", "ddpcr",
                      "--out", os.path.join(out_root, "probe")],
             "requires": [os.path.join(design_out, "design_pairs.tsv"), target],
             "produces": [(os.path.join(out_root, "probe", "probe_qc.tsv"),
                           ["pair_id", "check", "severity"])]},
            {"name": "10_amplicon_melt",
             "argv": [PY, "cli/amplicon_melt.py", "--design", design_out,
                      "--target", target, "--out",
                      os.path.join(out_root, "melt")],
             "requires": [os.path.join(design_out, "design_pairs.tsv")],
             "produces": [(os.path.join(out_root, "melt", "amplicon_melt.tsv"),
                           ["id", "length_bp", "tm_c", "qpcr02"])]},
        ]

    steps += [
        {"name": "10b_build_annotations",
         "argv": [PY, "cli/build_annotations.py"],
         "requires": [os.path.join(REPO_ROOT, "OligoScout", "data", "annotations")],
         "note": "ham GTF/superDups yoksa adim atlanir (SKIP)"},
        {"name": "11_panel_qc",
         "argv": [PY, "cli/panel_qc.py", "--prefix", prefix_ctdna,
                  "--out", panel_qc_out],
         "requires": [prefix_ctdna + ".bindings.tsv",
                      prefix_ctdna + ".products.tsv"],
         "produces": [(os.path.join(panel_qc_out, "panel_qc.tsv"),
                       ["pair_id", "tm_f_c", "sites_f"])]},
        {"name": "12_multiplex_qc",
         "argv": [PY, "cli/multiplex_qc.py", "--prefix", prefix_ctdna,
                  "--out", os.path.join(out_root, "multiplex")],
         "requires": [prefix_ctdna + ".products.tsv"]},
        {"name": "13_variant_qc",
         "argv": [PY, "cli/variant_qc.py", "--panel", panel_ctdna,
                  "--design", rel("out", "design_ctdna"), "--offline",
                  "--out", os.path.join(out_root, "variant")],
         "requires": [panel_ctdna,
                      rel("out", "design_ctdna", "design_panel.tsv")],
         "produces": [(os.path.join(out_root, "variant", "variant_qc.tsv"),
                       ["pair_id", "role", "check", "severity"])],
         "note": "--offline: yalniz yerel varyant onbellegi"},
        {"name": "14_rerank",
         "argv": [PY, "cli/rerank.py", "--design", rel("out", "design_ctdna"),
                  "--qc", panel_qc_out, "--out",
                  os.path.join(out_root, "rerank")],
         "requires": [os.path.join(panel_qc_out, "panel_qc.tsv"),
                      rel("out", "design_ctdna", "design_pairs.tsv")]},
        {"name": "15_compare_profiles",
         "argv": [PY, "cli/compare_profiles.py",
                  "--panel", f"ctdna={panel_ctdna}",
                  "--panel", f"qpcr_sybr={panel_qpcr}",
                  "--out", os.path.join(out_root, "compare")],
         "requires": [panel_ctdna, panel_qpcr],
         "produces": [(os.path.join(out_root, "compare",
                                    "profile_comparison.tsv"), [])]},
        {"name": "16_threshold_evidence",
         "argv": [PY, "cli/threshold_evidence.py", "--latest",
                  args.latest or "", "--out",
                  os.path.join(out_root, "kalibrasyon")],
         "requires": [args.latest] if args.latest else ["__yok__"],
         "note": "--latest verilmezse atlanir"},
        {"name": "16b_cpp_tablo_guncel",
         "argv": [PY, os.path.join(REPO_ROOT, "OligoScout", "tools",
                                   "gen_thermo_tables.py"), "--check"],
         "requires": [os.path.join(REPO_ROOT, "OligoScout", "tools",
                                   "gen_thermo_tables.py")],
         "note": "C++ tablo basligi params.py ile tutuyor mu"},
        {"name": "17_diff_cpp_sequences",
         "argv": [PY, "cli/diff_cpp_sequences.py", "--out", seq_list],
         "produces": [(seq_list, ["id", "sequence", "partner"])]},
        {"name": "18_diff_cpp",
         "argv": [PY, "cli/diff_cpp.py", "--sequences", seq_list,
                  "--dump", dump, "--out", os.path.join(out_root, "diff_cpp")],
         "requires": [dump, seq_list],
         "produces": [(os.path.join(out_root, "diff_cpp", "diff_cpp.tsv"),
                       ["id", "column", "sinif"])],
         "newer_than": [os.path.join(REPO_ROOT, "OligoScout", "src",
                                     "native_thermo.cpp"),
                        os.path.join(REPO_ROOT, "OligoScout", "include",
                                     "primerpair", "generated",
                                     "thermo_tables.inc")],
         "newer_key": dump,
         "note": "dump yoksa ya da C++ kaynagindan eskiyse atlanir; "
                 "uretmek icin RUN_DIFF_CPP.cmd"},
    ]
    return steps


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=rel("out", "regresyon"))
    ap.add_argument("--timeout", type=float, default=900.0)
    ap.add_argument("--latest", help="son panelin selected_panel.tsv dosyasi")
    ap.add_argument("--hizli", action="store_true",
                    help="uzun tasarim adimlarini atla")
    args = ap.parse_args()

    out_root = os.path.abspath(args.out)
    os.makedirs(out_root, exist_ok=True)
    steps = build_steps(out_root, args)

    results = []
    for step in steps:
        name = step["name"]
        missing = [p for p in step.get("requires", [])
                   if not (p and os.path.exists(p))]
        if not missing and step.get("newer_than"):
            key = step["newer_key"]
            newest = max((os.path.getmtime(x) for x in step["newer_than"]
                          if os.path.exists(x)), default=0.0)
            if os.path.getmtime(key) < newest:
                results.append({"name": name, "status": "SKIP", "code": "",
                                "sure_s": "",
                                "neden": "C++ dokumu kaynaktan eski; "
                                         "once RUN_DIFF_CPP.cmd"})
                print(f"SKIP {name}: dokum bayat")
                continue
        if missing:
            results.append({"name": name, "status": "SKIP", "code": "",
                            "sure_s": "", "neden": "girdi yok: " +
                            ", ".join(os.path.basename(m) or m
                                      for m in missing[:2])})
            print(f"SKIP {name}: girdi yok")
            continue

        log_path = os.path.join(out_root, f"{name}.log")
        started = time.time()
        try:
            done = subprocess.run(step["argv"], cwd=HERE, capture_output=True,
                                  text=True, timeout=args.timeout)
            code, stdout, stderr = done.returncode, done.stdout, done.stderr
        except subprocess.TimeoutExpired:
            code, stdout, stderr = 124, "", f"zaman asimi ({args.timeout}s)"
        elapsed = time.time() - started
        with open(log_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(f"$ {' '.join(step['argv'])}\n\n{stdout}\n--- stderr ---\n{stderr}\n")

        ok_codes = step.get("ok_codes", (0,))
        reason = ""
        status = "PASS" if code in ok_codes else "FAIL"
        if status == "PASS":
            for path, columns in step.get("produces", []):
                if not os.path.exists(path):
                    status, reason = "FAIL", f"cikti yok: {os.path.basename(path)}"
                    break
                if columns:
                    head = header_of(path)
                    eksik = [c for c in columns if c not in head]
                    if eksik:
                        status = "FAIL"
                        reason = (f"{os.path.basename(path)} basliginda eksik: "
                                  f"{', '.join(eksik)}")
                        break
        else:
            reason = (stderr.strip().splitlines() or ["cikis kodu " + str(code)])[-1][:160]

        results.append({"name": name, "status": status, "code": str(code),
                        "sure_s": f"{elapsed:.1f}", "neden": reason})
        print(f"{status} {name}  ({elapsed:.1f}s, cikis {code})"
              + (f"  {reason}" if reason else ""))

    n_fail = sum(1 for r in results if r["status"] == "FAIL")
    n_skip = sum(1 for r in results if r["status"] == "SKIP")
    n_pass = sum(1 for r in results if r["status"] == "PASS")

    tsv_path = os.path.join(out_root, "regresyon.tsv")
    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(results[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)

    md_path = os.path.join(out_root, "REGRESYON_RAPORU.md")
    with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# OligoScout uctan uca regresyon\n\n")
        fh.write(f"- python: {sys.version.split()[0]}\n")
        fh.write(f"- adim: {len(results)}  |  PASS {n_pass}  SKIP {n_skip}  "
                 f"FAIL {n_fail}\n\n")
        fh.write("| Adim | Sonuc | Cikis | Sure (s) | Not |\n|---|---|---:|---:|---|\n")
        for r in results:
            fh.write(f"| {r['name']} | {r['status']} | {r['code']} | "
                     f"{r['sure_s']} | {r['neden']} |\n")
        fh.write("\nPASS = komut 0 dondu ve beklenen dosyalari beklenen "
                 "basliklarla yazdi. SKIP = girdisi yok. Her adimin tam "
                 "ciktisi `<adim>.log` dosyasindadir.\n")
        if n_fail:
            fh.write("\n## Basarisiz adimlar\n\n")
            for r in results:
                if r["status"] == "FAIL":
                    fh.write(f"### {r['name']}\n\n```text\n")
                    log = os.path.join(out_root, f"{r['name']}.log")
                    if os.path.exists(log):
                        with open(log, encoding="utf-8") as lf:
                            fh.write("".join(lf.readlines()[-40:]))
                    fh.write("```\n\n")

    print(f"\nPASS {n_pass}  SKIP {n_skip}  FAIL {n_fail}")
    print(f"rapor: {md_path}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
