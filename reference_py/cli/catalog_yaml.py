#!/usr/bin/env python3
"""QC katalogu: YAML tek dogruluk kaynagi <-> uretilen TSV'ler.

    python3 cli/catalog_yaml.py migrate   # bir kerelik: mevcut TSV/JS'den catalog/qc_catalog.yaml uret
    python3 cli/catalog_yaml.py export    # YAML -> out/durum/qc_katalog_kapsami.tsv, web/src/qc_rules.tsv
    python3 cli/catalog_yaml.py check     # YAML gecerli mi, uretilen dosyalar guncel mi (cikis 1 = fark)

YAML alanlari (kural basina):
  id, family, status (etiketli_uygulandi | etiketli_kismi | ...), spec_title, where,
  evidence_level (regression | unit_test | archive | code_reference | none),
  evidence (regression/archive/test/code atiflari), calibration (REP-03 durumu),
  name {tr,en}, check {tr,en}, why {tr,en}, scope_note (varsa).
PyYAML yalniz bu arac icin gerekir; motor calisirken YAML okunmaz.
"""
from __future__ import annotations

import csv
import io
import re
import sys
from pathlib import Path

import yaml

PY = Path(__file__).resolve().parents[1]
ROOT = PY.parent
YAML_PATH = PY / "catalog" / "qc_catalog.yaml"
KAPSAM = PY / "out" / "durum" / "qc_katalog_kapsami.tsv"
AUDIT = PY / "out" / "qc_catalog_audit" / "qc_kural_dogrulama.tsv"
RULES = ROOT / "web" / "src" / "qc_rules.tsv"
SCOPE_JS = ROOT / "web" / "src" / "catalog.js"
CALIBRATION = {"status": "uncalibrated",
               "record": "REP-03",
               "note": "esikler literatur baslangic degerleri; islak laboratuvar kalibrasyonu yok"}
STATUSES = {"etiketli_uygulandi", "etiketli_kismi", "etiketli_kismi_tablo_pending", "wsl_script_hazir", "bekliyor"}
LEVELS = ("regression", "unit_test", "archive", "code_reference", "none")


def read_tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def split(v):
    return [x.strip() for x in (v or "").split(";") if x.strip()]


def scope_notes():
    text = SCOPE_JS.read_text(encoding="utf-8")
    block = text[text.index("const QC_SCOPE_NOTES"):text.index("};")]
    return {m.group(1): m.group(2) for m in re.finditer(r"'([A-Z0-9]+-[0-9]+)':\s*'((?:[^'\\]|\\.)*)'", block)}


def migrate():
    kapsam = {r["kural"]: r for r in read_tsv(KAPSAM)}
    audit = {r["kural"]: r for r in read_tsv(AUDIT)}
    rules = {r["id"]: r for r in read_tsv(RULES)}
    notes = scope_notes()
    out = []
    for rid in kapsam:
        k, a, t = kapsam[rid], audit.get(rid, {}), rules.get(rid, {})
        ev = {"regression": split(a.get("yeni_regresyon_kaydi")), "archive": split(a.get("arsiv_kaydi")),
              "unit_test": split(a.get("testte_kimlik_veya_sembol")), "code": split(a.get("kodda_kimlik_veya_sembol"))}
        level = ("regression" if ev["regression"] else "unit_test" if ev["unit_test"] else
                 "archive" if ev["archive"] else "code_reference" if ev["code"] else "none")
        item = {"id": rid, "family": k["aile"], "status": k["durum"], "spec_title": k["spec_basligi"],
                "where": k["nerede"], "evidence_level": level, "evidence": ev, "calibration": dict(CALIBRATION),
                "name": {"tr": t.get("name_tr", a.get("sade_ad", "")), "en": t.get("name_en", "")},
                "check": {"tr": t.get("check_tr", ""), "en": t.get("check_en", "")},
                "why": {"tr": t.get("why_tr", ""), "en": t.get("why_en", "")}}
        if rid in notes:
            item["scope_note"] = notes[rid].replace("\\'", "'")
        out.append(item)
    doc = {"schema": "oligoscout-qc-catalog/1", "rules": out}
    YAML_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(YAML_PATH, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# OligoScout QC katalogu — TEK DOGRULUK KAYNAGI.\n"
                 "# Uretilen dosyalar: python3 reference_py/cli/catalog_yaml.py export\n")
        yaml.safe_dump(doc, fh, allow_unicode=True, sort_keys=False, width=110)
    print(f"yazildi: {YAML_PATH} ({len(out)} kural)")


def load():
    doc = yaml.safe_load(YAML_PATH.read_text(encoding="utf-8"))
    rules = doc["rules"]
    errors = []
    ids = [r["id"] for r in rules]
    if len(ids) != len(set(ids)):
        errors.append("yinelenen kural kimligi")
    for r in rules:
        for f in ("id", "family", "status", "evidence_level", "name", "calibration"):
            if f not in r:
                errors.append(f"{r.get('id')}: '{f}' eksik")
        if r.get("status") not in STATUSES:
            errors.append(f"{r['id']}: bilinmeyen durum {r.get('status')}")
        if r.get("evidence_level") not in LEVELS:
            errors.append(f"{r['id']}: bilinmeyen kanit duzeyi {r.get('evidence_level')}")
    if errors:
        raise SystemExit("YAML gecersiz:\n  " + "\n  ".join(errors))
    return rules


def render(rules):
    k = io.StringIO()
    w = csv.writer(k, delimiter="\t", lineterminator="\n")
    w.writerow(["kural", "aile", "spec_basligi", "durum", "nerede"])
    for r in rules:
        w.writerow([r["id"], r["family"], r["spec_title"], r["status"], r["where"]])
    q = io.StringIO()
    w = csv.writer(q, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    w.writerow(["id", "name_tr", "check_tr", "why_tr", "name_en", "check_en", "why_en"])
    for r in sorted(rules, key=lambda x: x["id"]):
        w.writerow([r["id"], r["name"]["tr"], r["check"]["tr"], r["why"]["tr"],
                    r["name"]["en"], r["check"]["en"], r["why"]["en"]])
    return {KAPSAM: k.getvalue(), RULES: q.getvalue()}


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "check"
    if cmd == "migrate":
        migrate()
        return 0
    rules = load()
    files = render(rules)
    if cmd == "export":
        for path, text in files.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            bom = path == KAPSAM
            path.write_text(("﻿" if bom and False else "") + text, encoding="utf-8")
            print(f"yazildi: {path}")
        return 0
    stale = []
    for path, text in files.items():
        cur = path.read_text(encoding="utf-8-sig").replace("\r\n", "\n") if path.exists() else ""
        if cur != text:
            stale.append(str(path))
    levels = {}
    for r in rules:
        levels[r["evidence_level"]] = levels.get(r["evidence_level"], 0) + 1
    print(f"{len(rules)} kural; kanit duzeyi {levels}; durum "
          f"{ {s: sum(r['status'] == s for r in rules) for s in sorted({r['status'] for r in rules})} }")
    if stale:
        print("GUNCEL DEGIL: " + ", ".join(stale) + " -> catalog_yaml.py export")
        return 1
    print("uretilen dosyalar YAML ile tutarli")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
