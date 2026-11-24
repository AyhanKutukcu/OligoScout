"""OligoScout tool surface for MCP / REST.

Design rules (from the evidence-record philosophy of the engine):
* Every number in a response comes from an engine record; summaries are
  generated deterministically from those records, never by the language model.
* Anything that could not be evaluated is returned explicitly in
  ``not_evaluated`` with a reason -- never silently skipped.
* Every response carries provenance: parameter-set SHA-256, engine version,
  reference build and the REP-03 calibration notice.
"""
from __future__ import annotations

import collections
import csv
import gzip
import json
import os
import re
import tempfile
import threading
import time
import uuid
from dataclasses import asdict
from pathlib import Path

from protocol import Prompt, Registry, Resource, Tool, ToolInputError

HERE = Path(__file__).resolve().parent
ROOT = Path(os.getenv("OLIGOSCOUT_ROOT", HERE.parent)).resolve()
import sys  # noqa: E402

sys.path.insert(0, str(ROOT / "web"))
sys.path.insert(0, str(ROOT / "reference_py"))

import local_engine as LE  # noqa: E402  (the same backend used by oligoscout.org)
from local_reference import CHROMOSOMES, parse_region  # noqa: E402
from oligoscout_thermo import (align, engine, nn, oligomod, params, profiles,  # noqa: E402
                               salt)

if os.getenv("OLIGOSCOUT_PANEL_BINARY"):
    LE.BINARY = Path(os.environ["OLIGOSCOUT_PANEL_BINARY"])

VERSION = "0.4.0"
CALIBRATION = ("REP-03: most thresholds are literature starting values, NOT calibrated against "
               "wet-lab data; LNA contributions are not modelled (Tm is a lower bound). "
               "Research use only - confirm experimentally.")
PROFILE_IDS = sorted(profiles.PROFILES)
SALT_MODELS = list(salt.MODELS)
SEV_ORDER = {"FAIL": 3, "WARN": 2, "INFO": 1, "PASS": 1}

# ---------------------------------------------------------------------------
# engine wrapper
# ---------------------------------------------------------------------------
class Engine(LE.LocalEngine):
    """LocalEngine with two extensions: raw-sequence templates and optional genome screen."""

    CUSTOM = "custom"

    def __init__(self):
        super().__init__()
        self._custom = None
        self._scan_ok = None
        self.scan_enabled = True

    # raw FASTA support: the design routine reads its template through read_region
    def read_region(self, region):
        if region == "@custom" and self._custom is not None:
            return self.CUSTOM, 0, self._custom
        return super().read_region(region)

    def variant_window(self, chrom, start, end):
        if chrom == self.CUSTOM:
            return None, {}, "", "custom sequence has no genome coordinates"
        return super().variant_window(chrom, start, end)

    def scan_available(self):
        if self._scan_ok is None:
            ok = LE.BINARY.is_file() and os.access(LE.BINARY, os.X_OK) and self.index_status["primary"]["ready"]
            if ok:
                try:
                    import subprocess
                    out = subprocess.run([str(LE.BINARY)], capture_output=True, text=True, timeout=20)
                    ok = "Usage" in (out.stdout + out.stderr)
                except (OSError, subprocess.SubprocessError):
                    ok = False
            self._scan_ok = bool(ok)
        return self._scan_ok

    def screen_design(self, source, candidates, count, work, stage):
        if self.scan_enabled and self.scan_available() and source["template"]["chrom"] != self.CUSTOM:
            return super().screen_design(source, candidates, count, work, stage)
        reason = ("genome scan disabled by request" if not self.scan_enabled else
                  "custom sequence: no genome coordinates" if source["template"]["chrom"] == self.CUSTOM else
                  "C++ oligoscout_panel or PPFM index not runnable on this host")
        stats = {"scanned": 0, "passed": 0, "rejected": 0, "unscanned": len(candidates),
                 "requested": count, "skipped_reason": reason}
        return candidates[:count], stats, {"pairs": {}, "skipped": reason}


_ENGINE = None
_LOCK = threading.Lock()
_DESIGN_LOCK = threading.Lock()


def eng():
    global _ENGINE
    with _LOCK:
        if _ENGINE is None:
            _ENGINE = Engine()
        return _ENGINE


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def provenance(extra=None):
    e = eng()
    p = {"param_set_sha256": e.digest, "engine_version": engine.ENGINE_VERSION, "mcp_version": VERSION,
         "reference": "GRCh38.p14 primary (24 chromosomes)" if e.reference else None,
         "calibration": CALIBRATION}
    if extra:
        p.update(extra)
    return p


def worst(records):
    best = "INFO"
    for r in records:
        s = r.get("severity")
        if SEV_ORDER.get(s, 0) > SEV_ORDER.get(best, 0):
            best = s
    return best


def counts(records):
    c = collections.Counter(r.get("severity") for r in records if r.get("severity"))
    return {k: c.get(k, 0) for k in ("FAIL", "WARN", "INFO")}


def brief(value, limit=220):
    if isinstance(value, dict):
        value = {k: v for k, v in value.items() if k not in ("alignment", "p_unpaired_by_base", "curve")}
    s = json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return s if len(s) <= limit else s[:limit] + "…"


def flagged(records, top=12):
    rows = [r for r in records if r.get("severity") in ("WARN", "FAIL")]
    rows.sort(key=lambda r: -SEV_ORDER[r["severity"]])
    return [{"id": r.get("id"), "role": r.get("role"), "severity": r["severity"],
             "value": brief(r.get("value")), "threshold": r.get("threshold"), "note": r.get("note", "")}
            for r in rows[:top]]


def conditions(args):
    return LE.reaction(args.get("conditions") or {})


def clean_dna(seq, label, lo=8, hi=80, modified=False):
    if not isinstance(seq, str) or not seq.strip():
        raise ToolInputError(f"{label} sequence required")
    raw = re.sub(r"\s+", "", seq)
    parsed = oligomod.parse(raw) if modified else None
    dna = parsed.dna if parsed else raw.upper()
    if not re.fullmatch(r"[ACGT]+", dna):
        bad = sorted(set(re.sub(r"[ACGT]", "", dna)))
        raise ToolInputError(f"{label}: unsupported characters {bad} (SEQ-01; use A/C/G/T, LNA as +N)")
    if not lo <= len(dna) <= hi:
        raise ToolInputError(f"{label}: length must be {lo}-{hi} nt (got {len(dna)})")
    return dna, parsed


def parse_fasta(text, limit=60_000):
    if not isinstance(text, str) or not text.strip():
        raise ToolInputError("template sequence is empty")
    name, seq = "", []
    for line in text.splitlines():
        if line.startswith(">"):
            if name or seq:
                break
            name = line[1:].strip()
        else:
            seq.append(re.sub(r"[\s\d]", "", line))
    s = "".join(seq).upper().replace("U", "T")
    bad = set(re.sub(r"[ACGTN]", "", s))
    if bad:
        raise ToolInputError(f"template: unsupported characters {sorted(bad)}")
    if not 40 <= len(s) <= limit:
        raise ToolInputError(f"template must be 40-{limit} nt")
    return name or "custom", s


def template_source(args, profile="sanger"):
    """Resolve a GRCh38 region or a raw sequence into a LocalEngine 'source'."""
    e = eng()
    cond = conditions(args)
    model = args.get("salt_model", "owczarzy_2008")
    cons = asdict(LE.get_constraints({}, profiles.get(profile)))
    if args.get("region"):
        try:
            chrom, offset, seq = e.read_region(args["region"])
        except ValueError as exc:
            raise ToolInputError(str(exc))
        label = None
    elif args.get("template_sequence"):
        label, seq = parse_fasta(args["template_sequence"])
        chrom, offset = Engine.CUSTOM, 0
    else:
        return None
    return {"template": LE.template_data(chrom, offset, seq, "primary"), "conditions": cond,
            "profile": profile, "salt_model": model, "constraints": cons, "label": label}


def locate(source, fwd, rev, pid="Q1"):
    seq = source["template"]["sequence"]
    a = seq.find(fwd)
    rc = nn.reverse_complement(rev)
    b = seq.find(rc, a + 1 if a >= 0 else 0)
    if a < 0 or b < 0:
        missing = [n for n, ok in (("forward", a >= 0), ("reverse", b >= 0)) if not ok]
        raise ToolInputError(f"{' and '.join(missing)} primer not found (exact match) on the template "
                             f"strand(s); SEQ-05 cannot be satisfied. Check orientation (reverse primer is "
                             f"given 5'->3' on the minus strand) or the region.")
    return {"id": pid, "F": fwd, "R": rev, "start": a, "end": b + len(rc)}


def genomic(source, i):
    t = source["template"]
    return f"{t['chrom']}:{t['offset'] + i + 1:,}" if t["chrom"] != Engine.CUSTOM else f"pos {i + 1}"


# ---------------------------------------------------------------------------
# jobs (long genome scans) -- keeps chat clients responsive
# ---------------------------------------------------------------------------
class Jobs:
    def __init__(self):
        self.items, self.lock = collections.OrderedDict(), threading.Lock()

    def run(self, label, fn, wait):
        jid = uuid.uuid4().hex[:16]
        rec = {"job_id": jid, "tool": label, "status": "running", "stage": "started", "started": time.time()}
        with self.lock:
            self.items[jid] = rec
            while len(self.items) > 50:
                self.items.popitem(last=False)

        def stage(text):
            rec["stage"] = text

        def work():
            try:
                rec["result"] = fn(stage)
                rec["status"] = "done"
            except ValueError as exc:
                rec["status"], rec["error"] = "error", str(exc)
            except Exception as exc:  # pragma: no cover
                rec["status"], rec["error"] = "error", f"engine error: {exc}"
            rec["finished"] = time.time()

        th = threading.Thread(target=work, daemon=True)
        th.start()
        th.join(wait)
        return rec

    def get(self, jid):
        with self.lock:
            return self.items.get(jid)


JOBS = Jobs()


def job_response(rec, render):
    if rec["status"] == "done":
        return render(rec["result"])
    if rec["status"] == "error":
        return {"text": f"Analysis failed: {rec['error']}", "error": True,
                "structured": {"status": "error", "error": rec["error"], "job_id": rec["job_id"]}}
    el = round(time.time() - rec["started"])
    return {"text": (f"Still running ({rec['stage']}, {el} s). Call oligoscout_job_status with "
                     f"job_id='{rec['job_id']}' to collect the result."),
            "structured": {"status": "running", "job_id": rec["job_id"], "stage": rec["stage"],
                           "elapsed_s": el, "provenance": provenance()}}


# ---------------------------------------------------------------------------
# catalog
# ---------------------------------------------------------------------------
def _load_catalog():
    rows = {}
    kat = ROOT / "reference_py" / "out" / "durum" / "qc_katalog_kapsami.tsv"
    if kat.is_file():
        with kat.open(encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                parts = (r.get("spec_basligi") or "").split(" / ")
                rows[r["kural"]] = {"id": r["kural"], "family": r["aile"], "scope": parts[0] if parts else "",
                                    "title": parts[1] if len(parts) > 1 else r["kural"],
                                    "method": " / ".join(parts[2:]), "status": r["durum"],
                                    "code": r.get("nerede", ""), "calibration": "literature_default (uncalibrated)"}
    audit = ROOT / "reference_py" / "out" / "qc_catalog_audit" / "qc_kural_dogrulama.tsv"
    if audit.is_file():
        with audit.open(encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["kural"] in rows:
                    rows[r["kural"]]["plain_name_tr"] = r.get("sade_ad", "")
                    rows[r["kural"]]["tests"] = r.get("testte_kimlik_veya_sembol", "")
    cat_js = ROOT / "web" / "src" / "catalog.js"
    if cat_js.is_file():
        for rid, note in re.findall(r"'([A-Z]+-\d+)':\s*'([^']+)'", cat_js.read_text(encoding="utf-8")):
            if rid in rows:
                rows[rid]["scope_note"] = note
    for rid, note in {"REP-03": "Calibration status declaration (every record carries it)",
                      "VAR-08": "Amplicon-level variant summary (indels/length change inside the product)"}.items():
        if rid in rows and rows[rid]["title"] == rid:
            rows[rid]["title"] = note
    return rows


CATALOG = None


def catalog():
    global CATALOG
    if CATALOG is None:
        CATALOG = _load_catalog()
    return CATALOG


# ---------------------------------------------------------------------------
# shared schemas
# ---------------------------------------------------------------------------
COND_SCHEMA = {"type": "object", "description": "Reaction conditions (defaults: 50 mM Na+/K+, 1.5 mM Mg2+, "
               "0.6 mM dNTP, 200 nM primer, Ta 60 C).", "additionalProperties": False, "properties": {
                   "na": {"type": "number", "minimum": 0.01, "maximum": 1000, "description": "monovalent mM"},
                   "mg": {"type": "number", "minimum": 0, "maximum": 100, "description": "Mg2+ mM"},
                   "dntp": {"type": "number", "minimum": 0, "maximum": 100, "description": "total dNTP mM"},
                   "pnm": {"type": "number", "minimum": 0.01, "maximum": 100000, "description": "primer nM"},
                   "ta": {"type": "number", "minimum": 0, "maximum": 100, "description": "annealing C"}}}
SALT = {"type": "string", "enum": SALT_MODELS, "default": "owczarzy_2008"}
PROFILE = {"type": "string", "enum": PROFILE_IDS, "default": "qpcr_sybr",
           "description": "Application profile: sets product size, Tm/GC/length windows and active rule families."}
REGION = {"type": "string", "maxLength": 60, "description": "GRCh38 region, 1-based inclusive, e.g. "
          "'chr7:55191719-55191900' (max 5000 bp). Use oligoscout_resolve_target to turn a gene symbol into coordinates."}
TEMPLATE = {"type": "string", "maxLength": 70000, "description": "Alternative to region: FASTA or raw DNA "
            "(40-60000 nt). No genome coordinates -> no gnomAD layer; genome scan still possible for the primers."}
DNA = {"type": "string", "minLength": 8, "maxLength": 200}


# ---------------------------------------------------------------------------
# tool handlers
# ---------------------------------------------------------------------------
def t_status(args):
    e = eng()
    s = e.status()
    cat = catalog()
    impl = sum(r["status"] == "etiketli_uygulandi" for r in cat.values())
    data = {"reference_ready": s["reference"]["ready"], "genome_scan_available": e.scan_available(),
            "indexes": {k: v.get("ready") for k, v in s["indexes"].items()},
            "annotations": s["annotations"], "gnomad": {"api": s["variants"]["api"],
                                                        "cached_variants": s["variants"]["count"]},
            "profiles": PROFILE_IDS, "salt_models": SALT_MODELS,
            "qc_catalog": {"rules": len(cat), "implemented": impl,
                           "partial": [r for r, v in cat.items() if v["status"] != "etiketli_uygulandi"]},
            "provenance": provenance()}
    txt = (f"OligoScout engine {engine.ENGINE_VERSION} · params {e.digest[:12]}… · reference "
           f"{'ready' if data['reference_ready'] else 'MISSING'} · genome scan "
           f"{'available' if data['genome_scan_available'] else 'NOT available'} · QC catalog {impl}/{len(cat)} "
           f"implemented. {CALIBRATION}")
    return {"text": txt, "structured": data}


def t_profiles(args):
    out = [profiles.get(p).as_dict() for p in PROFILE_IDS]
    lines = [f"- {p['profile']}: {p['description']} · product {p['product_bp'][0]}-{p['product_bp'][1]} bp · "
             f"Tm {p['tm_window_c'][0]}-{p['tm_window_c'][1]} C · length {p['length_nt'][0]}-{p['length_nt'][1]} nt"
             for p in out]
    return {"text": "Application profiles:\n" + "\n".join(lines), "structured": {"profiles": out}}


def t_rule(args):
    cat = catalog()
    rid, q, fam = (args.get("rule_id") or "").upper().strip(), (args.get("query") or "").lower(), args.get("family")
    if rid:
        r = cat.get(rid)
        if not r:
            raise ToolInputError(f"unknown rule {rid}; rules look like TH-02, SPC-03, VAR-01")
        hits = [r]
    else:
        hits = [r for r in cat.values() if (not fam or r["family"].lower() == fam.lower()) and
                (not q or q in json.dumps(r, ensure_ascii=False).lower())][:25]
    lines = [f"{r['id']} [{r['family']}] {r['title']} — {r['method']} · status {r['status']} · "
             f"{r['calibration']}" + (f" · scope note: {r['scope_note']}" if r.get("scope_note") else "")
             for r in hits]
    return {"text": "\n".join(lines) or "No rule matched.",
            "structured": {"rules": hits, "total_rules": len(cat), "calibration": CALIBRATION}}


_EXON_CACHE = {}
_STRAND = {}


def gene_info(gene):
    """(strand, MANE Select transcript id) for a gene from the local GENCODE v47 GTF (cached)."""
    key = gene.upper()
    if key in _STRAND:
        return _STRAND[key]
    gtf = LE.ANNOTATIONS / "gencode.v47.annotation.gtf.gz"
    strand = mane = None
    if gtf.is_file():
        needle = f'gene_name "{gene}";'
        lines = []
        try:
            import shutil
            import subprocess
            if shutil.which("zcat") and shutil.which("grep"):
                out = subprocess.run(f"zcat '{gtf}' | grep -F -i '{needle}' | grep -m3 -E '\\tgene\\t|MANE_Select'",
                                     shell=True, capture_output=True, text=True, timeout=180).stdout
                lines = out.splitlines()
            else:
                with gzip.open(gtf, "rt", encoding="utf-8") as fh:
                    for line in fh:
                        if needle.lower() in line.lower():
                            lines.append(line)
                            if "MANE_Select" in line:
                                break
        except Exception:  # pragma: no cover
            lines = []
        for line in lines:
            cols = line.split("\t")
            if len(cols) > 8:
                strand = strand or cols[6]
                if "MANE_Select" in cols[8] and mane is None:
                    m = re.search(r'transcript_id "([^"]+)"', cols[8])
                    mane = m.group(1) if m else None
    _STRAND[key] = (strand, mane)
    return _STRAND[key]


def t_resolve(args):
    gene = (args.get("gene") or "").strip()
    if not re.fullmatch(r"[A-Za-z0-9_.\-]{1,40}", gene):
        raise ToolInputError("gene symbol required, e.g. EGFR")
    path = LE.ANNOTATIONS / "transcript_exons.tsv"
    if not path.is_file():
        return {"text": "GENCODE exon table not installed (run RUN_ANNOTATIONS.cmd).", "error": True,
                "structured": {"not_evaluated": {"resolve": "transcript_exons.tsv missing"}}}
    key = gene.upper()
    if key not in _EXON_CACHE:
        rows = []
        with path.open(encoding="utf-8") as fh:
            next(fh)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 5 and c[3].upper() == key:
                    rows.append((c[0], int(c[1]), int(c[2]), c[4]))
        _EXON_CACHE[key] = rows
    rows = _EXON_CACHE[key]
    if not rows:
        return {"text": f"Gene '{gene}' not found in GENCODE v47 (check the official HGNC symbol).",
                "structured": {"gene": gene, "transcripts": []}}
    tx = collections.defaultdict(list)
    for chrom, a, b, tid in rows:
        tx[tid].append((a, b))
    exon = args.get("exon")
    strand, mane = gene_info(gene)
    out = []
    for tid, ex in sorted(tx.items(), key=lambda kv: (kv[0] != mane, -len(kv[1]))):
        ex.sort()
        bio = ex[::-1] if strand == "-" else ex  # transcript (biological) exon order
        item = {"transcript_id": tid, "mane_select": tid == mane, "chrom": rows[0][0], "strand": strand, "exons": len(ex),
                "span": f"{rows[0][0]}:{ex[0][0] + 1}-{ex[-1][1]}",
                "exon_regions": [f"{rows[0][0]}:{a + 1}-{b}" for a, b in bio][:60]}
        if exon:
            if strand is None:
                item["exon_note"] = "strand unknown: exon numbering not reliable, ask the user"
            elif 1 <= exon <= len(bio):
                a, b = bio[exon - 1]
                item["exon_region"] = f"{rows[0][0]}:{a + 1}-{b}"
                item["exon_note"] = f"exon {exon} in transcript order ({strand} strand)"
            else:
                item["exon_note"] = f"transcript has only {len(bio)} exons"
        out.append(item)
    top = out[0]
    txt = (f"{gene}: {len(out)} GENCODE v47 transcripts on {top['chrom']} ({strand or '?'} strand); "
           f"{'MANE Select' if top.get('mane_select') else 'longest (no MANE Select found)'} {top['transcript_id']} "
           f"({top['exons']} exons, {top['span']})." +
           (f" Exon {exon} (transcript order): {top.get('exon_region', 'n/a')}." if exon else "") +
           " Transcripts differ in exon content; confirm the transcript with the user if it matters." +
           " Pass a region of <=5000 bp to the design tools.")
    return {"text": txt, "structured": {"gene": gene, "transcripts": out[:20], "source": "GENCODE v47 (local)"}}


def t_thermo(args):
    cond = conditions(args)
    model = args["salt_model"]
    eobj = engine.ThermoEngine(cond=LE.condition_model(cond), salt_model=model)
    results, lines = [], []
    for i, raw in enumerate(args["sequences"], 1):
        dna, parsed = clean_dna(raw, f"sequence {i}", modified=True)
        r = eobj.evaluate(f"S{i}", dna)
        lna = len(parsed.lna_positions) if parsed else 0
        th = r["thermo"]
        recs = [{**c, "role": f"S{i}"} for c in r["checks"]]
        spread = next((c["value"] for c in r["checks"] if c["id"] == "TH-02"), {})
        item = {"id": f"S{i}", "sequence": dna, "length_nt": len(dna), "gc_pct": round(100 * r["gc_fraction"], 2),
                "tm_c": th["tm_c"], "tm_1M_c": th["tm_1M_c"], "dH_kcal_mol": th["dH_kcal_mol"],
                "dS_cal_mol_K": th["dS_1M_cal_mol_K"], "dG37_kcal_mol": th["dG37_kcal_mol"],
                "salt_model": th["salt_model"], "salt_model_spread": spread,
                "worst_severity": r["worst_severity"], "flags": flagged(recs),
                "records": recs if args.get("detail") == "full" else None,
                "not_evaluated": r.get("not_evaluated", {})}
        if lna:
            item["lna_count"] = lna
            item["not_evaluated"] = {**item["not_evaluated"], "LNA": "LNA stabilisation not modelled; Tm is a DNA-backbone LOWER BOUND"}
        results.append(item)
        lines.append(f"S{i} {dna} ({len(dna)} nt, GC {item['gc_pct']}%): Tm {th['tm_c']:.2f} C "
                     f"[{model}], ΔG37 {th['dG37_kcal_mol']:.2f} kcal/mol, worst {r['worst_severity']}"
                     + (f", {lna} LNA -> lower bound" if lna else ""))
    return {"text": "\n".join(lines) + f"\nConditions {cond}. {CALIBRATION}",
            "structured": {"results": results, "conditions": cond, "provenance": provenance()}}


def t_dimer(args):
    a, _ = clean_dna(args["sequence_a"], "sequence_a", modified=True)
    b = clean_dna(args["sequence_b"], "sequence_b", modified=True)[0] if args.get("sequence_b") else a
    out = {}
    for mode in ("any", "3prime"):
        d = align.best_dimer(a, b, mode)
        out[mode] = d.as_dict() if d.found else {"found": False}
    txt = []
    for mode, d in out.items():
        if d.get("found") is False:
            txt.append(f"{mode}: no duplex of >=2 bp")
        else:
            al = "\n".join(d["alignment"])
            txt.append(f"{mode}: ΔG37 {d['dG37_kcal_mol']:.2f} kcal/mol, helix {d['helix_len']} bp, "
                       f"3'-paired {d['three_prime_paired']}, mismatches {d['mismatch_count']}, bulges {d['bulge_count']}\n{al}")
    thr = {"any": {"warn": -6.0, "fail": -9.0}, "3prime": {"warn": -5.0, "fail": -7.0}}
    sev = {m: ("FAIL" if d.get("dG37_kcal_mol", 0) < thr[m]["fail"] else "WARN" if d.get("dG37_kcal_mol", 0) < thr[m]["warn"]
               else "INFO") for m, d in out.items()}
    return {"text": f"{'Self' if a == b else 'Hetero'}-dimer (DIM-0{1 if a == b else 3}/0{2 if a == b else 4}):\n"
                    + "\n\n".join(txt) + f"\nSeverity: any {sev['any']}, 3' {sev['3prime']}. Model: WC stacks + "
                    "<=1 internal mismatch (Allawi-Peyret) or one 1-nt bulge; 37 C, 1 M NaCl reference.",
            "structured": {"a": a, "b": b, "self": a == b, "dimers": out, "severity": sev, "thresholds_kcal_mol": thr,
                           "provenance": provenance()}}


def _render_design(res, top):
    pairs = res.get("pairs", [])
    spec = res.get("specificity") or {}
    rows, lines = [], []
    src = {"template": res["template"]}
    for p in pairs[:top]:
        d = p.get("design", {})
        qc = res.get("qc", {}).get(p["id"], {})
        recs = qc.get("records", [])
        sp = (spec.get("pairs") or {}).get(p["id"], {})
        summ = sp.get("summary") or {}
        var = qc.get("variants") or {}
        var_recs = var.get("records", []) if var.get("evaluated") else []
        row = {"id": p["id"], "forward": p["F"], "reverse": p["R"],
               "forward_pos": genomic(src, p["start"]), "reverse_pos": genomic(src, p["end"] - 1),
               "product_bp": p["end"] - p["start"], "tm_f_c": round(d.get("left", {}).get("tm_c", float("nan")), 2),
               "tm_r_c": round(d.get("right", {}).get("tm_c", float("nan")), 2),
               "delta_tm_c": round(d.get("delta_tm", 0), 3), "hetero_dimer_dG37": d.get("hetero_dg"),
               "pareto_rank": d.get("pareto_rank"), "uncalibrated_score": round(d.get("score", 0), 4),
               "qc_worst": worst(recs + var_recs), "qc_counts": counts(recs + var_recs),
               "flags": flagged(recs + var_recs, 8),
               "specificity": ({"products_total": summ.get("products_total"), "offtarget_total": summ.get("offtarget_total"),
                                "exact_sites_f": summ.get("exact_f"), "exact_sites_r": summ.get("exact_r"),
                                "spc01": summ.get("spc01"), "spc03": summ.get("spc03")} if summ else None),
               "variants": ({"evaluated": True, "source": var.get("source"), "worst": worst(var_recs),
                             "n_in_amplicon": len(var.get("variants", []))} if var.get("evaluated")
                            else {"evaluated": False, "reason": var.get("reason")})}
        rows.append(row)
        lines.append(f"| {p['id']} | {p['F']} | {p['R']} | {row['product_bp']} | {row['tm_f_c']}/{row['tm_r_c']} | "
                     f"{'specific' if summ and summ.get('offtarget_total') == 0 and summ.get('products_total') == 1 else ('n/a' if not summ else 'CHECK')} | "
                     f"{row['variants'].get('worst', 'n/e')} | {row['qc_worst']} |")
    ne = {}
    if spec.get("skipped"):
        ne["SPC-01..09"] = spec["skipped"]
    if not res.get("variant_coverage"):
        ne["VAR-01..08"] = res.get("variant_reason") or "no gnomAD coverage for this template"
    head = (f"Template {res['template']['region']} · profile {res['profile']} · {res['counts']['left']} F / "
            f"{res['counts']['right']} R candidates · {res['counts']['pairs']} pairs passed hard filters · "
            f"genome screen: {res.get('screening', {}).get('passed', 0)} specific of "
            f"{res.get('screening', {}).get('scanned', 0)} scanned.")
    table = ("| id | forward 5'->3' | reverse 5'->3' | bp | Tm F/R | specificity | variants | QC |\n|---|---|---|---|---|---|---|---|\n"
             + "\n".join(lines)) if lines else "No pair satisfied all constraints."
    return {"text": head + "\n" + table + (f"\nNot evaluated: {ne}" if ne else "") + f"\n{CALIBRATION}",
            "structured": {"template": {k: res["template"][k] for k in ("region", "chrom", "offset", "fasta_sha256")},
                           "profile": res["profile"], "conditions": res["conditions"], "constraints": res["constraints"],
                           "counts": res["counts"], "rejections": res.get("rejections"), "screening": res.get("screening"),
                           "pairs": rows, "not_evaluated": ne, "provenance": provenance()}}


def t_design(args):
    e = eng()
    profile = args["profile"]
    data = {"profile": profile, "top": args["top"], "conditions": args.get("conditions") or {},
            "salt_model": args["salt_model"], "variant_filter": args["avoid_variants"],
            "variant_af": args["variant_af"], "variant_scope": args["variant_scope"], "constraints": {},
            "search_mode": args["search_mode"]}
    cons = data["constraints"]
    if args.get("product_size"):
        cons["product_min"], cons["product_max"] = args["product_size"]
    if args.get("tm"):
        cons["tm_min"], cons["tm_opt"], cons["tm_max"] = args["tm"]
    if args.get("length"):
        cons["length_min"], cons["length_max"] = args["length"]
    for k in ("forward_range", "reverse_range", "must_cover", "target_range"):
        if args.get(k):
            data[{"forward_range": "left_region", "reverse_range": "right_region"}.get(k, k)] = args[k]
    if args.get("region"):
        data["region"] = args["region"]
    elif args.get("template_sequence"):
        name, seq = parse_fasta(args["template_sequence"])
        data["region"] = "@custom"
    else:
        raise ToolInputError("give either region or template_sequence")

    def run(stage):
        with _DESIGN_LOCK:  # the engine object holds per-call state for custom templates/scan toggle
            if data["region"] == "@custom":
                e._custom = seq
            e.scan_enabled = args["specificity"]
            try:
                with tempfile.TemporaryDirectory(prefix="oligoscout-mcp-") as tmp:
                    return LE.clean_json(e.design(data, stage, Path(tmp)))
            finally:
                e._custom, e.scan_enabled = None, True

    rec = JOBS.run("design", run, args["wait_seconds"])
    return job_response(rec, lambda res: _render_design(res, args["top"]))


def _pair_inputs(args):
    f = clean_dna(args["forward"], "forward", 12, 60)[0]
    r = clean_dna(args["reverse"], "reverse", 12, 60)[0]
    return f, r


def t_check_pair(args):
    f, r = _pair_inputs(args)
    src = template_source(args, args["profile"])
    if src is None:
        # primer-only QC: thermodynamics, dimers, composition (no template rules)
        eobj = engine.ThermoEngine(cond=LE.condition_model(conditions(args)), salt_model=args["salt_model"])
        recs = []
        for role, s in (("F", f), ("R", r)):
            recs += [{**c, "role": role} for c in eobj.evaluate("Q_" + role, s)["checks"]]
        recs += [{**c, "role": "pair"} for c in eobj.evaluate_pair("Q", f, r)["checks"]]
        ne = {"SEQ-04/05, PAIR-02/03/04, AMP-*, SS-04, VAR-*": "no template given (pass region or template_sequence)",
              "SPC-*": "use oligoscout_check_specificity"}
        return {"text": f"Primer-only QC: worst {worst(recs)}; counts {counts(recs)}.\n" +
                        "\n".join(f"- {x['id']} {x['role']} {x['severity']}: {x['value']}" for x in flagged(recs)) +
                        f"\nNot evaluated: {ne}\n{CALIBRATION}",
                "structured": {"forward": f, "reverse": r, "worst_severity": worst(recs), "counts": counts(recs),
                               "flags": flagged(recs, 30), "records": recs if args["detail"] == "full" else None,
                               "not_evaluated": ne, "provenance": provenance()}}
    pair = locate(src, f, r)
    e = eng()
    out = LE.clean_json(e.pair_qc(src, pair))
    recs = out["records"]
    var = out.get("variants") or {}
    vrecs = var.get("records", []) if var.get("evaluated") else []
    allr = recs + vrecs
    ne = {} if var.get("evaluated") else {"VAR-01..08": var.get("reason", "custom template")}
    ne["SPC-01..09"] = "genome-wide specificity is a separate call: oligoscout_check_specificity"
    txt = (f"Pair on {src['template']['region']} · product {pair['end'] - pair['start']} bp "
           f"({genomic(src, pair['start'])} -> {genomic(src, pair['end'] - 1)}) · Tm F {out['thermo']['F']['tm_c']:.2f} / "
           f"R {out['thermo']['R']['tm_c']:.2f} C · worst {worst(allr)} · {counts(allr)}\n" +
           "\n".join(f"- {x['id']} {x['role']} {x['severity']}: {x['value']} (threshold {x['threshold']})"
                     for x in flagged(allr)) + f"\nNot evaluated: {ne}\n{CALIBRATION}")
    return {"text": txt, "structured": {
        "pair": {**pair, "product_bp": pair["end"] - pair["start"]}, "template": src["template"]["region"],
        "thermo": out["thermo"], "worst_severity": worst(allr), "counts": counts(allr), "flags": flagged(allr, 40),
        "records": allr if args["detail"] == "full" else None,
        "variants": {"evaluated": var.get("evaluated", False), "source": var.get("source"),
                     "variants": var.get("variants", [])[:60], "reason": var.get("reason")},
        "not_evaluated": ne, "provenance": provenance()}}


def t_specificity(args):
    e = eng()
    if not e.scan_available():
        ne = {"SPC-01..09": "C++ oligoscout_panel or the 24-shard PPFM index is not runnable on this host"}
        return {"text": "Genome-wide specificity is not available on this server. " + ne["SPC-01..09"],
                "error": True, "structured": {"status": "not_evaluated", "not_evaluated": ne, "provenance": provenance()}}
    selected, sources = [], {}
    base = template_source(args, args["profile"])
    for i, p in enumerate(args["pairs"], 1):
        f, r = clean_dna(p["forward"], f"pairs[{i}].forward", 12, 60)[0], clean_dna(p["reverse"], f"pairs[{i}].reverse", 12, 60)[0]
        src = template_source(p, args["profile"]) if p.get("region") else base
        pid = re.sub(r"[^A-Za-z0-9_]", "_", p.get("id") or f"P{i:02d}")[:24]
        if src is not None:
            item = locate(src, f, r, pid)
        else:
            src = {"template": LE.template_data("NA", 0, "N", "primary"), "conditions": conditions(args),
                   "profile": args["profile"], "salt_model": "owczarzy_2008",
                   "constraints": asdict(LE.get_constraints({}, profiles.get(args["profile"])))}
            item = {"id": pid, "F": f, "R": r, "start": 0, "end": 0}
        item["target_id"] = pid
        sources[pid] = src
        selected.append(item)
    if len({p["id"] for p in selected}) != len(selected):
        raise ToolInputError("pair ids must be unique")
    first = next(iter(sources.values()))
    source = {**first, "target_sources": sources, "conditions": conditions(args)}
    opts = {"anchor": args["anchor_nt"], "mismatches": args["max_mismatches"],
            "min_product": args["product_range"][0], "max_product": args["product_range"][1],
            "search_mode": args["search_mode"]}

    def run(stage):
        with tempfile.TemporaryDirectory(prefix="oligoscout-mcp-scan-") as tmp:
            return LE.clean_json(e.scan(source, selected, Path(tmp), stage, opts))

    def render(res):
        rows, lines = [], []
        for p in selected:
            pr = res["pairs"].get(p["id"], {})
            s = pr.get("summary") or {}
            off = pr.get("offtarget", [])
            if not p["end"]:
                n_exact = sum(1 for o in off if o.get("total_mismatches") == 0)
                verdict = (f"NO REGION GIVEN: {s.get('products_total')} product(s) genome-wide, {n_exact} perfect-match; "
                           "the intended one cannot be identified without a region")
            else:
                verdict = ("SPECIFIC" if s.get("products_total") == 1 and len(pr.get("intended", [])) == 1
                           and s.get("offtarget_total") == 0 else
                           "OFF-TARGET PRODUCTS" if s.get("offtarget_total") else
                           "NO INTENDED PRODUCT FOUND" if not pr.get("intended") else "SEE PRODUCTS")
            prim = {role: {"sites": pr.get("bindings", {}).get(role, {}).get("n"),
                           "histogram": pr.get("bindings", {}).get(role, {}).get("hist"),
                           "closest_sites": pr.get("bindings", {}).get(role, {}).get("top", [])[:8]}
                    for role in ("F", "R")}
            rows.append({"id": p["id"], "verdict": verdict, "summary": s, "intended": pr.get("intended"),
                         "offtarget_products": off[:20], "primers": prim,
                         "records": [x for x in pr.get("records", [])] if args["detail"] == "full" else flagged(pr.get("records", []) +
                                     [x for role in ("F", "R") for x in pr.get("primers", {}).get(role, [])], 15)})
            lines.append(f"- {p['id']}: {verdict} · products {s.get('products_total')} (off-target "
                         f"{s.get('offtarget_total')}, exact off-target {s.get('offtarget_exact')}, single-primer "
                         f"{s.get('offtarget_single_primer', 'n/a')}) · F sites "
                         f"{prim['F']['sites']} · R sites {prim['R']['sites']}" +
                         "".join(f"\n    off-target {o['chrom']}:{o['amplicon_start']:,} {o['amplicon_length']} bp, "
                                 f"{o['total_mismatches']} mm, {o.get('orientation', 'F-R')}"
                                 + (" (single-primer product)" if o.get("orientation") in ("F-F", "R-R") else "")
                                 for o in off[:5]))
        c = res.get("contracts", {}).get("primary", {})
        rule = (f"SENSITIVE terminal rule (<= {opts['mismatches']} mismatches/primer, <= 1 in the 3'-terminal 5 nt, "
                f"no exact anchor)" if opts["search_mode"] == "sensitive" else
                f"exact 3' anchor {opts['anchor']} nt, <= {opts['mismatches']} mismatches/primer")
        txt = (f"Genome-wide scan of GRCh38.p14 (24 shards) · contract: {rule}, "
               f"products {opts['min_product']}-{opts['max_product']} bp · "
               f"index manifest {str(c.get('manifest_sha256', ''))[:12]}…\n" + "\n".join(lines) +
               ("\nProducts: F-R, R-F and single-primer F-F / R-R." if c.get("single_primer_products") == "on" else
                  "\nProducts: F-R and R-F only. NOT EVALUATED: single-primer products (F-F, R-R) - this "
                  "oligoscout_panel build predates them; rebuild with setup/run_panel.sh.")
               + (" Sensitive mode equals the brute-force Primer-BLAST-rule site set exactly (benchmark).\n"
                  if opts["search_mode"] == "sensitive" else
                  " Blind spot by contract: sites whose mismatches fall inside the 3' anchor are not enumerated; "
                  "rerun with search_mode='sensitive' for Primer-BLAST-rule coverage.\n")
               + CALIBRATION)
        return {"text": txt, "structured": {"pairs": rows, "contract": {**opts, **c},
                                            "multiplex": res.get("multiplex"), "provenance": provenance()}}

    rec = JOBS.run("specificity", run, args["wait_seconds"])
    return job_response(rec, render)


def t_variants(args):
    f, r = _pair_inputs(args)
    src = template_source(args, "sanger")
    if src is None or src["template"]["chrom"] == Engine.CUSTOM:
        raise ToolInputError("variant check needs a GRCh38 region (gnomAD coordinates)")
    pair = locate(src, f, r)
    out = LE.clean_json(eng().pair_variants(src, pair))
    if not out.get("evaluated"):
        return {"text": f"Variants not evaluated: {out.get('reason')}", "structured": {
            "evaluated": False, "not_evaluated": {"VAR-01..08": out.get("reason")}, "provenance": provenance()}}
    recs = out["records"]
    t = src["template"]
    vs = []
    for v in out["variants"]:
        idx = v["pos"] - 1 - t["offset"]
        where = ("F" if pair["start"] <= idx < pair["start"] + len(f) else
                 "R" if pair["end"] - len(r) <= idx < pair["end"] else "amplicon")
        d3 = (pair["start"] + len(f) - 1 - idx) if where == "F" else (idx - (pair["end"] - len(r))) if where == "R" else None
        vs.append({**{k: v.get(k) for k in ("rsid", "pos", "ref", "alt", "af", "genome_af", "exome_af")},
                   "under": where, "distance_from_3prime_nt": d3})
    vs.sort(key=lambda x: -(x["af"] or 0))
    risky = [x for x in vs if x["under"] != "amplicon" and (x["af"] or 0) >= 1e-3]
    def line(x):
        d = x["distance_from_3prime_nt"]
        where = f"under {x['under']}" + ("" if d is None else f", {d} nt from the 3' end")
        return f"- {x['rsid'] or x['pos']} {x['ref']}>{x['alt']} AF {(x['af'] or 0):.2e} {where}"
    txt = (f"gnomAD v4 ({out['source']}) · {len(vs)} variants in the amplicon · worst {worst(recs)}\n" +
           "\n".join(line(x) for x in (risky or vs[:5])) + f"\n{CALIBRATION}")
    return {"text": txt, "structured": {"evaluated": True, "source": out["source"], "worst_severity": worst(recs),
                                        "records": recs, "variants": vs[:80], "provenance": provenance()}}


def t_probe(args):
    f, r = _pair_inputs(args)
    src = template_source(args, "ddpcr")
    if src is None:
        raise ToolInputError("probe design needs region or template_sequence")
    pair = locate(src, f, r)
    data = {"probe": args.get("probe") or ""}
    out = LE.clean_json(eng().single("probe", src, pair, data))
    if out.get("status") == "no_candidate":
        return {"text": f"No probe candidate passed PRB rules; rejections {out['rejections']}",
                "structured": {"status": "no_candidate", **out, "provenance": provenance()}}
    res = out["result"]
    recs = res.get("checks", [])
    txt = (f"Probe {res.get('sequence')} ({res.get('length_nt')} nt) · Tm {res.get('tm_c')} C · ΔTm vs primers "
           f"{res.get('delta_tm_c')} · worst {res.get('worst_severity')}\n" +
           "\n".join(f"- {x['id']} {x['severity']}: {x['value']}" for x in flagged(recs)) + f"\n{CALIBRATION}")
    return {"text": txt, "structured": {"probe": res, "rejections": out.get("rejections"), "provenance": provenance()}}


def t_blocker(args):
    data = {k: args.get(k, "") for k in ("mutant", "wt", "forward", "reverse", "blocker")}
    data["conditions"] = args.get("conditions") or {}
    with tempfile.TemporaryDirectory(prefix="oligoscout-mcp-ct05-") as tmp:
        out = LE.clean_json(eng().oligo("blocker", data, Path(tmp)))
    recs = out["records"]
    ct = next((x for x in recs if x["id"] == "CT-05"), {})
    v = ct.get("value", {})
    txt = (f"CT-05 allele discrimination: mutant probe Tm {v.get('mutant_tm_c_dna_backbone')} C, WT probe "
           f"{v.get('wt_tm_c_dna_backbone')} C, ΔTm {v.get('delta_tm_c')} C -> {ct.get('severity')}. LNA modelled: "
           f"{v.get('lna_modelled')} (Tm values are DNA-backbone LOWER BOUNDS).\n" +
           "\n".join(f"- {x['id']} {x.get('role', '')} {x['severity']}: {brief(x['value'])}" for x in recs if x["id"] != "CT-05")
           + f"\n{CALIBRATION}")
    return {"text": txt, "structured": {"records": recs, "provenance": provenance()}}


def t_melt(args):
    f, r = _pair_inputs(args)
    src = template_source(args, "qpcr_sybr")
    if src is None:
        raise ToolInputError("melt needs region or template_sequence")
    pair = locate(src, f, r)
    out = LE.clean_json(eng().single("melt", src, pair, {}))
    s, chk = out["summary"], out["check"]
    curve = [c for c in out["curve"] if 60 <= c[0] <= 95][::2]
    txt = (f"Amplicon {pair['end'] - pair['start']} bp melt (Poland-type model): Tm {s.get('tm_c')} C, peaks "
           f"{s.get('peaks_c')} -> QPCR-02 {chk.get('severity')} ({chk.get('note', '')}). Loop entropy tables not "
           f"applied; sigma uncalibrated.\n{CALIBRATION}")
    return {"text": txt, "structured": {"summary": s, "check": chk, "curve_T_helicity_dthetadT": curve,
                                        "provenance": provenance()}}


def t_multiplex(args):
    from multiplex_design import run as mrun
    data = {"targets": args["targets"], "profile": args["profile"], "candidates": args["candidates_per_target"],
            "pool": args.get("pool") or {}, "conditions": args.get("conditions") or {},
            "salt_model": args["salt_model"], "constraints": {}}
    if args.get("product_size"):
        data["constraints"]["product_min"], data["constraints"]["product_max"] = args["product_size"]
    e = eng()
    if not e.scan_available():
        return {"text": "Multiplex design requires the genome scan (C++ + PPFM index), not available here.",
                "error": True, "structured": {"not_evaluated": {"MPX-*": "genome scan unavailable"}}}

    def run(stage):
        with tempfile.TemporaryDirectory(prefix="oligoscout-mcp-mpx-") as tmp:
            return LE.clean_json(mrun(e, data, stage, Path(tmp)))

    def render(res):
        txt = json.dumps({k: res.get(k) for k in ("pool", "summary", "selected") if k in res}, ensure_ascii=False)[:3500]
        return {"text": "Multiplex pool result (engine output):\n" + txt + f"\n{CALIBRATION}",
                "structured": {"result": {k: v for k, v in res.items() if k not in ("templates",)},
                               "provenance": provenance()}}

    rec = JOBS.run("multiplex", run, args["wait_seconds"])
    return job_response(rec, render)


def t_job(args):
    rec = JOBS.get(args["job_id"])
    if rec is None:
        raise ToolInputError("unknown or expired job_id (jobs live only in this server process)")
    renders = {"design": lambda res: _render_design(res, 50)}
    if rec["status"] != "done" or rec["tool"] not in renders:
        if rec["status"] == "done":
            return {"text": "Done.", "structured": {"status": "done", "result": rec["result"]}}
        return job_response(rec, None)
    return job_response(rec, renders[rec["tool"]])


# ---------------------------------------------------------------------------
# registry
# ---------------------------------------------------------------------------
GUIDE = (HERE / "GUIDE.md")
INSTRUCTIONS = (
    "OligoScout: thermodynamic PCR primer/probe design and 91-rule QC on GRCh38. Workflow: "
    "(1) oligoscout_resolve_target to turn a gene symbol into GRCh38 coordinates (ask the user which "
    "transcript/exon if ambiguous - never guess); (2) oligoscout_design_primers with a region <=5000 bp and the "
    "right profile (qpcr_sybr, rt_qpcr, ddpcr, ctdna, sanger, ngs_amplicon, bisulfite, standard_pcr); "
    "(3) oligoscout_check_specificity / oligoscout_check_variants / oligoscout_design_probe as needed. "
    "Report numbers only from tool outputs; state every item in not_evaluated; always relay the REP-03 "
    "calibration notice: thresholds are uncalibrated literature defaults and results need wet-lab validation. "
    "Read resource oligoscout://guide for interpretation rules.")


def build_registry():
    reg = Registry("oligoscout", "OligoScout primer design & QC", VERSION, INSTRUCTIONS)
    obj = lambda props, req=(): {"type": "object", "additionalProperties": False, "properties": props, "required": list(req)}
    wait = {"type": "integer", "minimum": 0, "maximum": 600, "default": 50,
            "description": "Seconds to wait synchronously; if the job is still running a job_id is returned for oligoscout_job_status."}
    detail = {"type": "string", "enum": ["summary", "full"], "default": "summary",
              "description": "'full' returns every evidence record (large)."}
    pair_props = {"forward": {**DNA, "description": "forward primer 5'->3'"},
                  "reverse": {**DNA, "description": "reverse primer 5'->3' (binds the minus strand)"},
                  "region": REGION, "template_sequence": TEMPLATE, "conditions": COND_SCHEMA, "salt_model": SALT}
    span = lambda d: {"type": "array", "items": {"type": "integer", "minimum": 1}, "minItems": 2, "maxItems": 2, "description": d}

    reg.tool(Tool("oligoscout_status", "Engine status",
                  "Report engine/version, parameter digest, which data layers are installed (GRCh38, PPFM index, "
                  "GENCODE, pseudogene/segdup, gnomAD) and QC catalog coverage. Call first when unsure what is available.",
                  obj({}), t_status))
    reg.tool(Tool("oligoscout_list_profiles", "Application profiles",
                  "List the 8 application profiles with product size, Tm, length and GC windows.", obj({}), t_profiles))
    reg.tool(Tool("oligoscout_explain_rule", "Explain QC rules",
                  "Look up QC rules of the 91-rule catalog (19 families) by id (e.g. TH-02, SPC-03, VAR-01), free-text "
                  "query or family. Returns method, implementation status, code location and calibration level.",
                  obj({"rule_id": {"type": "string", "maxLength": 10}, "query": {"type": "string", "maxLength": 80},
                       "family": {"type": "string", "maxLength": 30}}), t_rule))
    reg.tool(Tool("oligoscout_resolve_target", "Resolve gene to GRCh38",
                  "Resolve a gene symbol to GENCODE v47 transcripts and GRCh38 coordinates (optionally one exon). "
                  "Deterministic lookup in the local annotation; does not guess.",
                  obj({"gene": {"type": "string", "maxLength": 40}, "exon": {"type": "integer", "minimum": 1, "maximum": 500}},
                      ["gene"]), t_resolve))
    reg.tool(Tool("oligoscout_thermo", "Oligo thermodynamics",
                  "Nearest-neighbour thermodynamics for 1-50 oligos: Tm under the chosen salt model (+ spread across "
                  "Owczarzy 2004/2008 and von Ahsen 2001), ΔH/ΔS/ΔG37, θ at Ta, 3'-end stability, self-dimer and "
                  "hairpin (TH-01..05, DIM-01/02, SS-01/03). LNA '+N' and end modifications accepted (LNA -> lower bound).",
                  obj({"sequences": {"type": "array", "minItems": 1, "maxItems": 50, "items": {"type": "string", "maxLength": 200}},
                       "conditions": COND_SCHEMA, "salt_model": SALT, "detail": detail}, ["sequences"]), t_thermo))
    reg.tool(Tool("oligoscout_dimer", "Dimer alignment",
                  "Most stable self- or hetero-dimer (any position and 3'-extendable) with an IDT-style alignment; "
                  "WC stacks + one internal mismatch or 1-nt bulge (DIM-01..04).",
                  obj({"sequence_a": DNA, "sequence_b": DNA}, ["sequence_a"]), t_dimer))
    reg.tool(Tool("oligoscout_design_primers", "Design primers",
                  "Primer-BLAST-style design on a GRCh38 region or a pasted sequence: candidate generation, hard "
                  "filters, Pareto ranking, then (if available) genome-wide specificity screening, gnomAD variant "
                  "layer and full pair QC for the top pairs. Long runs return a job_id.",
                  obj({"region": REGION, "template_sequence": TEMPLATE, "profile": PROFILE,
                       "product_size": span("[min,max] product bp"), "tm": {"type": "array", "items": {"type": "number"},
                       "minItems": 3, "maxItems": 3, "description": "[min, opt, max] Tm C"},
                       "length": span("[min,max] primer length nt"),
                       "forward_range": span("[from,to] 1-based on template where the forward primer must lie"),
                       "reverse_range": span("[from,to] 1-based where the reverse primer must lie"),
                       "must_cover": span("[from,to] 1-based region the product must include"),
                       "top": {"type": "integer", "minimum": 1, "maximum": 50, "default": 5},
                       "avoid_variants": {"type": "boolean", "default": True},
                       "variant_af": {"type": "number", "minimum": 0, "maximum": 1, "default": 0.001},
                       "variant_scope": {"type": "string", "enum": ["all", "3p"], "default": "all"},
                       "specificity": {"type": "boolean", "default": True,
                                       "description": "screen candidates genome-wide with the C++ FM-index"},
                       "search_mode": {"type": "string", "enum": ["anchored", "sensitive"], "default": "anchored",
                                       "description": "anchored: exact 3' anchor (fast default). sensitive: Primer-BLAST "
                                       "terminal rule - every site with <=max mismatches and <=1 in the 3'-terminal 5 nt"},
                       "conditions": COND_SCHEMA, "salt_model": SALT, "wait_seconds": wait}), t_design,
                  idempotent=True, open_world=True))
    reg.tool(Tool("oligoscout_check_pair", "QC a primer pair",
                  "Run the full pair QC (SEQ, CMP, TH, SS, DIM, PAIR, AMP, VAR when coordinates exist) on an existing "
                  "primer pair. With region/template the primers are located on it; without, primer-only rules run.",
                  obj({**pair_props, "profile": PROFILE, "detail": detail}, ["forward", "reverse"]), t_check_pair,
                  open_world=True))
    reg.tool(Tool("oligoscout_check_specificity", "Genome-wide specificity",
                  "Scan 1-50 primer pairs against all 24 GRCh38 chromosomes with the C++ FM-index (exact 3' anchor by "
                  "default, or search_mode='sensitive' for the Primer-BLAST terminal rule): binding sites, off-target products, SPC-01..09, multiplex cross-products. "
                  "Optional region per pair identifies the intended product.",
                  obj({"pairs": {"type": "array", "minItems": 1, "maxItems": 50, "items": obj(
                      {"id": {"type": "string", "maxLength": 24}, "forward": DNA, "reverse": DNA, "region": REGION},
                      ["forward", "reverse"])},
                       "region": REGION, "profile": PROFILE, "anchor_nt": {"type": "integer", "minimum": 8, "maximum": 24, "default": 12},
                       "max_mismatches": {"type": "integer", "minimum": 0, "maximum": 3, "default": 3},
                       "search_mode": {"type": "string", "enum": ["anchored", "sensitive"], "default": "anchored",
                                       "description": "anchored: exact 3' anchor (fast default). sensitive: Primer-BLAST "
                                       "terminal rule - every site with <=max mismatches and <=1 in the 3'-terminal 5 nt"},
                       "product_range": {**span("[min,max] product bp"), "default": [50, 3000]},
                       "conditions": COND_SCHEMA, "detail": detail, "wait_seconds": wait}, ["pairs"]),
                  t_specificity))
    reg.tool(Tool("oligoscout_check_variants", "gnomAD variants under primers",
                  "SNPCheck-style: gnomAD v4 variants under the primers and in the amplicon with allele frequency, "
                  "distance from the 3' end and VAR-01..08 verdicts. Requires a GRCh38 region.",
                  obj({**pair_props}, ["forward", "reverse", "region"]), t_variants, open_world=True))
    reg.tool(Tool("oligoscout_design_probe", "Hydrolysis probe",
                  "Select (or evaluate a given) TaqMan/ddPCR hydrolysis probe inside the amplicon: PRB-01..05, SS-05, "
                  "DIM-05.", obj({**pair_props, "probe": {"type": "string", "maxLength": 200}},
                                 ["forward", "reverse"]), t_probe))
    reg.tool(Tool("oligoscout_check_blocker", "Allele discrimination (CT-05)",
                  "Mutant vs wild-type probe/blocker discrimination (ΔTm) and probe/blocker x primer dimers. LNA as +N, "
                  "end mods as /56-FAM/ etc. LNA is not modelled: Tm values are lower bounds.",
                  obj({"mutant": {"type": "string", "maxLength": 200}, "wt": {"type": "string", "maxLength": 200},
                       "forward": {"type": "string", "maxLength": 200}, "reverse": {"type": "string", "maxLength": 200},
                       "blocker": {"type": "string", "maxLength": 200}, "conditions": COND_SCHEMA}, ["mutant", "wt"]),
                  t_blocker))
    reg.tool(Tool("oligoscout_amplicon_melt", "Amplicon melt curve",
                  "Predicted amplicon melting profile (-dθ/dT peaks) for SYBR qPCR single-peak check (QPCR-02/AMP-02).",
                  obj({**pair_props}, ["forward", "reverse"]), t_melt))
    reg.tool(Tool("oligoscout_design_multiplex", "Multiplex panel design",
                  "Design 2-8 targets together and choose a compatible pool (Tm spread, cross-dimers, cross-products, "
                  "size separation). Requires the genome scan.",
                  obj({"targets": {"type": "array", "minItems": 2, "maxItems": 8, "items": obj(
                      {"name": {"type": "string", "maxLength": 60}, "region": REGION}, ["name", "region"])},
                       "profile": PROFILE, "candidates_per_target": {"type": "integer", "minimum": 1, "maximum": 4, "default": 3},
                       "product_size": span("[min,max] product bp"),
                       "pool": {"type": "object", "properties": {"tm_spread": {"type": "number"}, "dimer_any": {"type": "number"},
                                                                 "dimer_3p": {"type": "number"}, "product_gap": {"type": "integer"}}},
                       "conditions": COND_SCHEMA, "salt_model": SALT, "wait_seconds": wait}, ["targets"]), t_multiplex,
                  open_world=True))
    reg.tool(Tool("oligoscout_job_status", "Job status",
                  "Collect the result of a long-running design/specificity/multiplex job by job_id.",
                  obj({"job_id": {"type": "string", "maxLength": 32}}, ["job_id"]), t_job))

    reg.resource(Resource("oligoscout://guide", "Interpretation guide",
                          "How to read OligoScout outputs, severities, contracts and calibration limits.",
                          "text/markdown", lambda _: GUIDE.read_text(encoding="utf-8") if GUIDE.is_file() else INSTRUCTIONS))

    def cat_reader(rid):
        c = catalog()
        if rid:
            r = c.get(rid.upper())
            return json.dumps(r, ensure_ascii=False, indent=1) if r else None
        return json.dumps(list(c.values()), ensure_ascii=False)
    reg.resource(Resource("oligoscout://catalog", "QC rule catalog", "All 91 QC rules (JSON).", "application/json", cat_reader))
    reg.resource(Resource("oligoscout://profiles", "Application profiles", "Profile windows (JSON).", "application/json",
                          lambda _: json.dumps(t_profiles({})["structured"], ensure_ascii=False)))
    bench = HERE / "benchmark_scores.json"
    reg.resource(Resource("oligoscout://benchmark", "Benchmark scorecard", "Validation and benchmark results (JSON).",
                          "application/json", lambda _: bench.read_text(encoding="utf-8") if bench.is_file() else "{}"))

    reg.prompt(Prompt("design_assay", "Plan and run a primer/probe assay design with OligoScout.",
                      [{"name": "target", "description": "gene, exon, variant or region", "required": True},
                       {"name": "application", "description": "qPCR, ddPCR, ctDNA, Sanger, NGS, bisulfite…", "required": False}],
                      lambda a: {"description": "OligoScout assay design", "messages": [{"role": "user", "content": {
                          "type": "text", "text":
                          f"Design an assay for {a.get('target', '?')} ({a.get('application', 'unspecified application')}) "
                          "with OligoScout. Steps: resolve the target to GRCh38 (ask me if transcript/exon is ambiguous), "
                          "choose the matching profile and confirm product size/conditions, run oligoscout_design_primers, "
                          "then summarise the top pairs with specificity, variant risk and every not_evaluated item, and "
                          "end with the REP-03 calibration notice. Do not invent numbers."}}]}))
    return reg
