"""Run the existing OligoScout engines against the local human reference data."""
from __future__ import annotations

import collections
import csv
from dataclasses import asdict, fields
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import urllib.error

from local_reference import CHROMOSOMES, Reference, inspect_shards, parse_region, read_tsv

ROOT = Path(__file__).resolve().parents[1]
PYROOT = ROOT / "reference_py"
DATA = ROOT / "OligoScout" / "data"
ANNOTATIONS = DATA / "annotations"
BINARY = ROOT / "OligoScout" / "build" / "release" / "oligoscout_panel"
CONTIG_BINARY = Path(os.environ.get("OLIGOSCOUT_CONTIG_SCAN_BINARY")
                     or ROOT / "OligoScout" / "build" / "release" / "oligoscout_contig_scan")
EXTRA_CONTIGS = ANNOTATIONS / "grch38_extra_contigs.fa"
EXTRA_CONTIGS_TSV = ANNOTATIONS / "grch38_extra_contigs.tsv"
CACHE = Path(os.getenv("OLIGOSCOUT_VARIANT_CACHE", PYROOT / "out" / "variant_cache.jsonl"))
PROJECT_PANEL = PYROOT / "out" / "design_ctdna" / "design_panel.tsv"
PROJECT_REGION = "chr22:20000001-20001000"
GNOMAD_DATASET = "gnomad_r4"
GNOMAD_CACHE = PYROOT / "out" / "gnomad_api"
GNOMAD_MAX_AGE = 30 * 24 * 60 * 60
GNOMAD_LEFT_PAD = 1000
sys.path.insert(0, str(PYROOT))

from cli.isoform_qc import load_exons_for_genes
from cli.variant_qc import load_cache, parse_variants, request_region
from oligoscout_thermo import (
    amplicon, bisulfite, design, engine, isoform, nn, oligomod, params,
    probe, profiles, salt, seqpair, variants,
)


def finite(value, label, lo, hi, *, integer=False):
    if type(value) not in (int, float) or not math.isfinite(value) or not lo <= value <= hi:
        raise ValueError(f"{label} {lo}–{hi} aralığında olmalı.")
    if integer and int(value) != value:
        raise ValueError(f"{label} tam sayı olmalı.")
    return int(value) if integer else float(value)


def reaction(data):
    if not isinstance(data, dict):
        raise ValueError("Reaksiyon koşulları geçersiz.")
    defaults = {"na": (50, .01, 1000), "mg": (1.5, 0, 100), "dntp": (.6, 0, 100),
                "pnm": (200, .01, 100000), "ta": (60, 0, 100)}
    return {k: finite(data.get(k, default), k, lo, hi) for k, (default, lo, hi) in defaults.items()}


def condition_model(data):
    return salt.Conditions(monovalent_M=data["na"] / 1e3, magnesium_M=data["mg"] / 1e3,
                           dntp_M=data["dntp"] / 1e3, primer_M=data["pnm"] / 1e9,
                           anneal_c=data["ta"])


def condition_args(data):
    return [option for key, name in (("na", "na"), ("mg", "mg"), ("dntp", "dntp"),
                                     ("pnm", "primer-nm"), ("ta", "ta"))
            for option in (f"--{name}", str(data[key]))]


def dna(value, label, lo=12, hi=80, *, modified=False, optional=False):
    if optional and not value:
        return ""
    if not isinstance(value, str) or not value.strip() or len(value) > hi * 10:
        raise ValueError(f"{label} dizisini girin.")
    cleaned = re.sub(r"\s+", "", value)
    seq = oligomod.parse(cleaned).dna if modified else nn.validate(cleaned)
    nn.validate(seq)
    if not lo <= len(seq) <= hi:
        raise ValueError(f"{label} {lo}–{hi} baz olmalı.")
    return cleaned if modified else seq


def optional_interval(data, key, length):
    span = data.get(key)
    if span is None:
        return None
    if not isinstance(span, list) or len(span) != 2:
        raise ValueError("Konum kısıtı iki koordinat içermeli.")
    a = finite(span[0], "Aralık başı", 1, length, integer=True)
    b = finite(span[1], "Aralık sonu", a, length, integer=True)
    return a - 1, b


def get_constraints(settings, profile):
    if not isinstance(settings, dict):
        raise ValueError("Tasarım ayarları geçersiz.")
    pc = profile.constraints or design.DesignConstraints()
    ranges = {
        "length_min": (12, 60), "length_max": (12, 60), "tm_min": (0, 99),
        "tm_max": (0, 99), "tm_opt": (0, 99), "gc_min": (0, 1), "gc_max": (0, 1),
        "max_homopolymer": (2, 20), "gc_clamp_min": (0, 5), "gc_clamp_max": (0, 5),
        "self_dimer_dg": (-50, 0), "self_dimer3p_dg": (-50, 0),
        "hetero_dimer_dg": (-50, 0), "hetero_dimer3p_dg": (-50, 0),
        "product_min": (20, 5000), "product_max": (20, 5000), "max_delta_tm": (0, 40),
        "theta_min": (0, 1),
    }
    values = asdict(pc)
    for key, (lo, hi) in ranges.items():
        if key in settings:
            values[key] = finite(settings[key], key, lo, hi, integer=type(values[key]) is int)
    for key in ("length", "tm", "gc", "product", "gc_clamp"):
        if values[key + "_min"] > values[key + "_max"]:
            raise ValueError(f"{key}: alt sınır üst sınırı aşıyor.")
    if not values["tm_min"] <= values["tm_opt"] <= values["tm_max"]:
        raise ValueError("Optimum Tm, Tm penceresi içinde olmalı.")
    allowed = {item.name for item in fields(design.DesignConstraints)}
    return design.DesignConstraints(**{k: v for k, v in values.items() if k in allowed})


def clean_json(obj):
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if isinstance(obj, dict):
        return {key: clean_json(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean_json(value) for value in obj]
    return obj


def jsonl(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def write_tsv(path, names, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(names)
        writer.writerows(rows)


def normalized_row(row):
    return {key: parse_number(value) for key, value in row.items()}


def parse_number(value):
    if value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return float(value)
        except (TypeError, ValueError):
            return value


def run_command(command, timeout=720):
    try:
        process = subprocess.run(command, cwd=PYROOT, capture_output=True,
                                 text=True, errors="replace", timeout=timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        raise ValueError("Motor ayrılan süre içinde tamamlanamadı.") from exc
    if process.returncode:
        log = (process.stderr or process.stdout).strip().splitlines()
        raise ValueError(f"Motor çalışması başarısız ({process.returncode}): {' '.join(log[-4:])[-800:]}")


def run_cli(name, *arguments):
    run_command([sys.executable, str(PYROOT / "cli" / name), *map(str, arguments)])


class LocalEngine:
    def __init__(self):
        self.reference = None
        self.reference_error = ""
        try:
            self.reference = Reference(DATA / "genomes" / "GRCh38.p14.fa",
                                       DATA / "benchmarks" / "reference" / "GRCh38.p14.primary24.manifest.tsv")
        except (OSError, ValueError, KeyError) as exc:
            self.reference_error = str(exc)
        self.indexes = {
            "primary": DATA / "indexes" / "ppfm",
            "CT": DATA / "indexes" / "bisulfite" / "CT",
            "GA": DATA / "indexes" / "bisulfite" / "GA",
        }
        self.index_status = {key: inspect_shards(path) for key, path in self.indexes.items()}
        self.cache_error = ""
        self.gnomad_api = os.getenv("OLIGOSCOUT_GNOMAD_API", "1") != "0"
        self.gnomad_retry_at = 0
        self.gnomad_error = ""
        self.gnomad_windows = collections.OrderedDict()
        try:
            markers, records = load_cache(CACHE)
            self.variants = list({
                (variants.normalize_chrom(v.chrom), v.pos, v.ref, v.alt): v for v in records
            }.values())
            self.coverage = []
            for key in markers:
                match = re.fullmatch(r"(?:chr)?([^:]+):(\d+)-(\d+)", key)
                if match:
                    chrom, start, end = match.groups()
                    self.coverage.append((variants.normalize_chrom(chrom), int(start), int(end)))
        except (OSError, ValueError, TypeError, KeyError) as exc:
            self.variants, self.coverage = [], []
            self.cache_error = str(exc)
        self.digest = params.registry_digest()

    def status(self):
        return {
            "reference": {"ready": self.reference is not None, "build": "GRCh38.p14",
                          "chromosomes": list(CHROMOSOMES), "error": self.reference_error},
            "cpp": BINARY.is_file() and os.access(BINARY, os.X_OK),
            "indexes": self.index_status, "project": PROJECT_PANEL.is_file(),
            "annotations": {
                "pseudogene": (ANNOTATIONS / "pseudogene_intervals.tsv").is_file(),
                "segdup": (ANNOTATIONS / "segdup_intervals.tsv").is_file(),
                "isoform": (ANNOTATIONS / "transcript_exons.tsv").is_file(),
                "chrM": (ANNOTATIONS / "chrM.fa").is_file(),
            },
            "variants": {
                "ready": self.gnomad_api or bool(self.coverage), "count": len(self.variants),
                "api": self.gnomad_api, "dataset": GNOMAD_DATASET,
                "regions": [{"chrom": f"chr{c}", "start": a + 1, "end": b}
                            for c, a, b in variants_union(self.coverage)],
                "error": self.cache_error,
            },
            "param_set_sha256": self.digest,
        }

    def read_region(self, region):
        if self.reference is None:
            raise ValueError("Yerel GRCh38 referansı hazır değil.")
        chrom, start, end = parse_region(region)
        return chrom, start, self.reference.read(chrom, start, end)

    def covered(self, chrom, start, end):
        chrom = variants.normalize_chrom(chrom)
        for c, a, b in variants_union(self.coverage):
            if c == chrom and a <= start and end <= b:
                return True
        return False

    def variant_window(self, chrom, start, end):
        chrom = variants.normalize_chrom(chrom)
        key = (chrom, start, end)
        now = time.time()
        cached = self.gnomad_windows.get(key)
        if cached is not None and now - cached[0] < GNOMAD_MAX_AGE:
            self.gnomad_windows.move_to_end(key)
            return cached[1]
        path = GNOMAD_CACHE / f"{GNOMAD_DATASET}-{chrom}-{start}-{end}.json"
        rows = None
        cached_at = now
        from_disk = False
        try:
            cached_at = path.stat().st_mtime
            if now - cached_at < GNOMAD_MAX_AGE:
                payload = json.loads(path.read_text(encoding="utf-8"))
                if payload.get("dataset") == GNOMAD_DATASET and isinstance(payload.get("rows"), list):
                    rows = payload["rows"]
                    from_disk = True
        except (OSError, ValueError, AttributeError):
            pass
        if rows is None and self.gnomad_api and now >= self.gnomad_retry_at:
            try:
                rows = request_region(chrom, max(1, start + 1 - GNOMAD_LEFT_PAD), end,
                                      GNOMAD_DATASET, timeout=12)
                cached_at = now
                self.gnomad_error = ""
                try:
                    GNOMAD_CACHE.mkdir(parents=True, exist_ok=True)
                    temp = path.with_suffix(".tmp")
                    temp.write_text(json.dumps({"dataset": GNOMAD_DATASET, "rows": rows},
                                               separators=(",", ":")), encoding="utf-8")
                    temp.replace(path)
                except OSError:
                    pass
            except (OSError, ValueError, RuntimeError, TypeError, KeyError) as exc:
                self.gnomad_retry_at = now + 30
                if isinstance(exc, urllib.error.HTTPError):
                    self.gnomad_error = f"gnomAD API HTTP {exc.code} yanıtı verdi."
                else:
                    self.gnomad_error = "gnomAD API sorgusu tamamlanamadı."
        if rows is not None:
            parsed = parse_variants(rows, chrom, highest_af=True)
            records = variants.overlapping(parsed, chrom, start, end)
            frequencies = {
                (int(r["pos"]), r["ref"], r["alt"]): {
                    k: (r.get(k) or {}).get("af") for k in ("genome", "exome")
                }
                for r in rows
            }
            source = "gnomAD v4 API önbelleği" if from_disk else "gnomAD v4 API"
            result = (records, frequencies, source, "")
            self.gnomad_windows[key] = (cached_at, result)
            self.gnomad_windows.move_to_end(key)
            if len(self.gnomad_windows) > 16:
                self.gnomad_windows.popitem(last=False)
            return result
        if self.covered(chrom, start, end):
            return variants.overlapping(self.variants, chrom, start, end), {}, "Yerel gnomAD önbelleği", ""
        return None, {}, "", self.gnomad_error or "Bu bölge için gnomAD kaynağına erişilemedi."

    def project(self):
        chrom, offset, seq = self.read_region(PROJECT_REGION)
        if not PROJECT_PANEL.is_file():
            raise ValueError("Kayıtlı chr22 proje paneli bulunamadı.")
        pairs = []
        for row in read_tsv(PROJECT_PANEL):
            f, r = row["primer1"], row["primer2"]
            start, end = int(row["left_start"]) - offset, int(row["right_start"]) - offset
            if seq[start:start + len(f)] != f or seq[end - len(r):end] != nn.reverse_complement(r):
                raise ValueError("Kayıtlı panel koordinatları yerel GRCh38 referansıyla uyuşmuyor.")
            pairs.append({"id": row["pair_id"], "F": f, "R": r, "start": start, "end": end})
        return {"template": template_data(chrom, offset, seq, "primary"), "pairs": pairs,
                "conditions": reaction({}), "profile": "ctdna", "salt_model": "owczarzy_2008",
                "constraints": asdict(profiles.get("ctdna").constraints)}

    def design(self, data, stage, work):
        chrom, offset, original = self.read_region(data.get("region"))
        profile = data.get("profile", "sanger")
        if not isinstance(profile, str) or profile not in profiles.PROFILES:
            raise ValueError("Geçerli bir uygulama profili seçin.")
        conversion = data.get("conversion", "primary") if profile == "bisulfite" else "primary"
        if conversion not in ("primary", "CT", "GA") or (profile == "bisulfite" and conversion == "primary"):
            raise ValueError("Bisülfit tasarımı için CT veya GA seçin.")
        seq = original.replace("C", "T") if conversion == "CT" else original
        if conversion == "GA":
            seq = original.replace("G", "A")
        submitted = data.get("reference_sequence")
        if submitted is not None and (not isinstance(submitted, str) or submitted.upper() != original):
            raise ValueError("FASTA dizisi girilen GRCh38.p14 aralığıyla tam eşleşmiyor.")
        params_data = get_constraints(data.get("constraints", {}), profiles.get(profile))
        conditions = reaction(data.get("conditions", {}))
        cond = condition_model(conditions)
        model = data.get("salt_model", "owczarzy_2008")
        if not isinstance(model, str) or model not in salt.MODELS:
            raise ValueError("Geçerli bir tuz modeli seçin.")
        length = len(seq)
        target_range = optional_interval(data, "target_range", length) or (0, length)
        ranges = {}
        for key in ("left_region", "right_region"):
            a, b = optional_interval(data, key, length) or target_range
            low, high = max(a, target_range[0]), min(b, target_range[1])
            if high - low < params_data.length_min:
                raise ValueError("Primer konumu seçilen şablon aralığında yeterli alan bırakmıyor.")
            ranges[key] = low, high
        must_cover = optional_interval(data, "must_cover", length)
        if must_cover and not (target_range[0] <= must_cover[0] < must_cover[1] <= target_range[1]):
            raise ValueError("Kapsanacak hedef, tasarım aralığı içinde olmalı.")
        n_top = finite(data.get("top", 20), "Çift sayısı", 1, 50, integer=True)
        threshold = finite(data.get("variant_af", .001), "AF", 0, 1)
        scope = data.get("variant_scope", "all")
        if scope not in ("all", "3p"):
            raise ValueError("Varyant filtre kapsamı geçersiz.")
        variant_source, variant_reason = "", ""
        region_variants = None
        if conversion == "primary":
            stage("gnomAD · bölge varyantları alınıyor")
            region_variants, _, variant_source, variant_reason = self.variant_window(chrom, offset, offset + length)
            if region_variants is None and self.covered(chrom, offset + target_range[0], offset + target_range[1]):
                region_variants = self.variants
                variant_source, variant_reason = "Yerel gnomAD önbelleği", ""
        has_variants = region_variants is not None
        use_variants = data.get("variant_filter") is True and has_variants
        filter_candidate = None
        if use_variants:
            def filter_candidate(start, end, orientation):
                lo, hi = offset + start, offset + end
                for v in variants.overlapping(region_variants, chrom, lo, hi):
                    if v.af < threshold:
                        continue
                    v_start, v_end = v.pos - 1, v.pos - 1 + max(1, len(v.ref))
                    near = v_end > hi - 5 if orientation == "F" else v_start < lo + 5
                    if near or scope == "all":
                        return ("VAR-02" if near else "VAR-03", "Varyant filtresi",
                                f"{v.rsid or v.pos} AF={v.af:g}")
                return None
        stage("Python · aday primerler ve dimer filtreleri")
        start_clock = time.monotonic()
        left, l_rej = design.enumerate_candidates(seq, "F", params_data, cond, ranges["left_region"],
                                                   salt_model=model, candidate_filter=filter_candidate)
        right, r_rej = design.enumerate_candidates(seq, "R", params_data, cond, ranges["right_region"],
                                                    salt_model=model, candidate_filter=filter_candidate)
        stage("Python · çift kurgulama ve Pareto sıralaması")
        pairs, p_rej = design.build_pairs(left, right, params_data, must_cover=must_cover)
        ranked = design.rank_pairs(pairs, params_data)
        rejected = collections.Counter(r.rule for entries in (l_rej, r_rej, p_rej) for r in entries)
        candidates = [
            {"id": f"WEB_{i:03d}", "F": p.left.sequence, "R": p.right.sequence,
             "start": p.product_start, "end": p.product_end, "design": asdict(p)}
            for i, p in enumerate(ranked, 1)
        ]
        result = {
            "template": template_data(chrom, offset, seq, conversion),
            "conditions": conditions, "profile": profile, "salt_model": model,
            "constraints": asdict(params_data),
            "counts": {"left": len(left), "right": len(right), "pairs": len(candidates)},
            "rejections": dict(rejected), "variant_filter": use_variants,
            "variant_coverage": has_variants, "variant_source": variant_source,
            "variant_reason": variant_reason,
            "candidate_cap_reached": len(ranked) >= params_data.max_pairs,
            "search_mode": search_mode(data),
        }
        result["pairs"], result["screening"], result["specificity"] = self.screen_design(
            result, candidates, n_top, work, stage
        )
        result["qc"] = {}
        if result["pairs"]:
            stage("Python · seçilen çiftlerin kalite kontrolü")
            for pair in result["pairs"]:
                layer = result["specificity"]["pairs"].get(pair["id"], {}).get("variants")
                result["qc"][pair["id"]] = self.pair_qc(result, pair, layer)
        result["elapsed_ms"] = round(1000 * (time.monotonic() - start_clock))
        return result

    def screen_design(self, source, candidates, count, work, stage):
        selected, batches = [], []
        checked = passed = 0
        batch_size = max(32, count * 4)
        t = source["template"]
        options = {
            "min_product": min(50, source["constraints"]["product_min"]),
            "max_product": max(3000, source["constraints"]["product_max"]),
        }
        while checked < len(candidates) and len(selected) < count:
            batch = candidates[checked:checked + batch_size]
            folder = work / f"batch-{checked}"
            folder.mkdir()
            progress = f"Aday {checked + 1}–{checked + len(batch)} / {len(candidates)}"
            panel, contracts, settings = self.genome_search(
                source, batch, folder, lambda name: stage(f"{name} · {progress}"), options
            )
            kept = specific_pairs(t, batch, folder)
            passed += len(kept)
            kept = kept[:count - len(selected)]
            selected.extend(kept)
            if kept:
                batches.append((folder, {p["id"] for p in kept}))
            checked += len(batch)
            batch_size = min(512, batch_size * 2)
        stats = {"scanned": checked, "rejected": checked - passed, "passed": passed,
                 "unscanned": len(candidates) - checked, "requested": count}
        if not selected:
            return selected, stats, None

        final = work / "selected"
        final.mkdir()
        panel = self.write_scan_panel(source, selected, final)
        for key in contracts:
            target = final / key
            target.mkdir()
            first = batches[0][0] / key
            counts = {}
            for suffix in ("bindings", "products"):
                path = target / f"panel.{suffix}.tsv"
                counts[f"{suffix[:-1]}_records"] = copy_selected_rows(path, key, suffix, batches)
            counts["panel_pairs"] = len(selected)
            # The final evidence tables contain only the selected subset of these runs.
            contract = {r["field"]: r["value"] for r in read_tsv(first / "panel.contract.tsv")}
            for field in ("search_seconds", "unique_primer_searches", "shard_loads", "shard_evictions"):
                contract.pop(field, None)
            contract.update(counts)
            write_tsv(target / "panel.contract.tsv", ["field", "value"], contract.items())
        contracts = self.scan_contracts(final, contracts)
        evidence = self.evaluate_scan(source, selected, final, stage, panel, contracts, settings)
        return selected, stats, evidence

    @staticmethod
    def pair_source(source, pair):
        return source.get("target_sources", {}).get(pair.get("target_id"), source)

    def pair_variants(self, source, pair):
        source = self.pair_source(source, pair)
        t = source["template"]
        a, b = t["offset"] + pair["start"], t["offset"] + pair["end"]
        if t["conversion"] != "primary":
            return {"evaluated": False, "reason": "Bu önbellek dönüşmemiş GRCh38 referansına ait."}
        pool, frequency, dataset, reason = self.variant_window(
            t["chrom"], t["offset"], t["offset"] + len(t["sequence"])
        )
        if pool is None and self.covered(t["chrom"], a, b):
            pool, dataset = self.variants, "Yerel gnomAD önbelleği"
        if pool is None:
            return {"evaluated": False, "reason": reason}
        context = {"source": dataset, "dataset": GNOMAD_DATASET,
                   "af_policy": "max_genome_exome" if frequency else "as_cached"}
        f, r = pair["F"], pair["R"]
        sites = (
            variants.Site(f"{pair['id']}_F", t["chrom"], a, a + len(f), "+", f, f),
            variants.Site(f"{pair['id']}_R", t["chrom"], b - len(r), b, "-", r, r),
        )
        records = []
        for role, site in zip(("F", "R"), sites):
            records.extend({**c.as_dict(), "role": role} for c in variants.evaluate_site(site, pool, conditions=context))
        records.extend({**c.as_dict(), "role": "amplicon"} for c in variants.evaluate_amplicon(
            pair["id"], t["chrom"], a, b, sites, pool, conditions=context
        ))
        hits = variants.overlapping(pool, t["chrom"], a, b)
        return {"evaluated": True, "source": dataset, "records": records,
                "variants": [{**asdict(v), "genome_af": frequency.get((v.pos, v.ref, v.alt), {}).get("genome"),
                              "exome_af": frequency.get((v.pos, v.ref, v.alt), {}).get("exome")}
                             for v in hits]}

    def pair_qc(self, source, pair, variant_layer=None):
        source = self.pair_source(source, pair)
        t = source["template"]
        cond = condition_model(source["conditions"])
        f, r, a, b = pair["F"], pair["R"], pair["start"], pair["end"]
        seq = t["sequence"]
        amp = seq[a:b]
        tm_window = source["constraints"]["tm_min"], source["constraints"]["tm_max"]
        eng = engine.ThermoEngine(cond=cond, profile=engine.Profile(
            name=source["profile"], tm_window_c=tm_window), salt_model=source["salt_model"])
        f_e, r_e = eng.evaluate(pair["id"] + "_F", f), eng.evaluate(pair["id"] + "_R", r)
        pair_e = eng.evaluate_pair(pair["id"], f, r)
        con = seqpair.SeqPairConstraints(
            product_bp=(source["constraints"]["product_min"], source["constraints"]["product_max"]))
        records = []
        for role, s, out in (("F", f, f_e), ("R", r, r_e)):
            records += [{**c, "role": role} for c in out["checks"]]
            records += [{**c, "role": role} for c in seqpair.check_primer(s, con)]
            records.append({**seqpair.check_seq01(s), "role": role})
        records += [{**c, "role": "pair"} for c in pair_e["checks"]]
        records += [{**c, "role": "pair"} for c in seqpair.check_pair(
            f, r, a, b, amp, cond=cond, salt_model=source["salt_model"], cons=con)]
        records += [{**c, "role": "template"} for c in seqpair.check_template(f, r, a, b, seq)]
        if "N" not in amp:
            records.append({**seqpair.check_amplicon_structure(amp, cond.anneal_c, con), "role": "amplicon"})
        records.append({**seqpair.check_seq06(
            t["region"] + " GRCh38.p14", t["fasta_sha256"], n_bases=len(seq)), "role": "template"})
        if source["profile"] == "rt_qpcr":
            records += [{**seqpair.check_rt05(s), "role": role} for role, s in (("F", f), ("R", r))]
        if t["conversion"] != "primary":
            original = self.reference.read(t["chrom"], t["offset"] + a, t["offset"] + b)
            for role, s, foot in (("F", f, original[:len(f)]),
                                  ("R", r, seqpair.reverse_complement_raw(original[-len(r):]))):
                records.extend({**c, "role": role} for c in (
                    bisulfite.check_bis01(s), bisulfite.check_bis02(foot), bisulfite.check_bis03(s)))
        else:
            extra = spc08_extra_contigs(f, r, source.get("search_mode") or "anchored")
            if extra is not None:
                records.extend(extra)
            else:
                mito = ANNOTATIONS / "chrM.fa"
                if mito.is_file():
                    mseq = "".join(mito.read_text().splitlines()[1:]).upper()
                    records.extend({**seqpair.check_spc08(s, {"chrM": mseq}), "role": role}
                                   for role, s in (("F", f), ("R", r)))
        return {
            "records": self.provenance(records, cond, source["profile"]),
            "thermo": {"F": f_e["thermo"], "R": r_e["thermo"]},
            "melting_curves": {role: eng.melting_curve(s, t_min=20, t_max=max(95, cond.anneal_c), step=0.5)
                               for role, s in (("F", f), ("R", r))},
            "variants": variant_layer if variant_layer is not None else self.pair_variants(source, pair),
        }

    def provenance(self, records, cond, profile):
        return [{"conditions": cond.as_dict(), "param_set_sha256": self.digest,
                 "profile": profile, "engine_version": engine.ENGINE_VERSION, **r}
                for r in records]

    def write_scan_panel(self, source, selected, work):
        t = source["template"]
        panel = work / "panel.tsv"
        write_tsv(panel, ["pair_id", "primer1", "primer2", "left_start", "right_start",
                          "expected_amplicon_length", "expected_chrom", "conversion"],
                  [(p["id"], p["F"], p["R"], target["offset"] + p["start"], target["offset"] + p["end"],
                    p["end"] - p["start"], target["chrom"], target["conversion"])
                   for p in selected for target in [self.pair_source(source, p)["template"]]])
        return panel

    def scan_contracts(self, work, keys):
        contracts = {}
        for key in keys:
            contract = {row["field"]: parse_number(row["value"])
                        for row in read_tsv(work / key / "panel.contract.tsv")}
            for field in ("manifest", "index_dir", "panel"):
                contract.pop(field, None)
            contract["manifest_sha256"] = self.index_status[key]["manifest_sha256"]
            contracts[key] = contract
        return contracts

    def genome_search(self, source, selected, work, stage, data):
        bis = source["template"]["conversion"] in ("CT", "GA")
        needed = ["CT", "GA"] if bis else ["primary"]
        if not (BINARY.is_file() and os.access(BINARY, os.X_OK)) or any(
            not self.index_status[k]["ready"] for k in needed
        ):
            raise ValueError("24 kromozomluk genom taraması şu anda kullanılamıyor.")
        anchor = finite(data.get("anchor", 12), "3′ çapa", 8, 24, integer=True)
        if any(anchor > min(len(p["F"]), len(p["R"])) for p in selected):
            raise ValueError("3′ çapa en kısa primerden uzun olamaz.")
        max_mm = finite(data.get("mismatches", 3), "Uyumsuzluk", 0, 3, integer=True)
        min_amp = finite(data.get("min_product", 50), "En kısa ürün", 20, 5000, integer=True)
        max_amp = finite(data.get("max_product", 3000), "En uzun ürün", min_amp, 5000, integer=True)
        mode = search_mode(data, source)
        panel = self.write_scan_panel(source, selected, work)
        for key in needed:
            stage(f"C++ · {key if bis else 'GRCh38'} · 24 kromozom taranıyor")
            index_dir = self.indexes[key]
            run_dir = work / key
            run_dir.mkdir()
            run_command([str(BINARY), str(index_dir / "manifest.sa8.tsv"), str(index_dir),
                         str(panel), str(run_dir / "panel"), "--anchor", str(anchor),
                         "--max-mismatch", str(max_mm), "--min-amplicon", str(min_amp),
                         "--max-amplicon", str(max_amp), "--cache", "1", "--mode", mode])
        return panel, self.scan_contracts(work, needed), (anchor, min_amp, max_amp)

    def scan(self, source, selected, work, stage, data):
        panel, contracts, settings = self.genome_search(source, selected, work, stage, data)
        return self.evaluate_scan(source, selected, work, stage, panel, contracts, settings)

    def evaluate_scan(self, source, selected, work, stage, panel, contracts, settings):
        t = source["template"]
        conversion = t["conversion"]
        bis = conversion in ("CT", "GA")
        anchor, min_amp, max_amp = settings
        result = {"bisulfite": bis, "reference": "GRCh38.p14", "pairs": {}, "contracts": contracts}
        if bis:
            stage("Python · BIS-04 iplik özgüllüğü")
            run_cli("bisulfite_qc.py", "--ct-prefix", work / "CT" / "panel",
                    "--ga-prefix", work / "GA" / "panel", "--out", work, "--intended", conversion, "--panel", panel)
            result["summary"] = [normalized_row(row) for row in read_tsv(work / "bisulfite_qc.tsv")]
            result["records"] = jsonl(work / "bisulfite_qc_records.jsonl")
            return result

        stage("Python · SPC, termodinamik ve anotasyonlar")
        folder = work / "primary"
        tmlo, tmhi = (source["constraints"][k] for k in ("tm_min", "tm_max"))
        run_cli("panel_qc.py", "--prefix", folder / "panel", "--out", folder,
                "--panel", panel, "--anchor", anchor, "--tm-window", f"{tmlo}:{tmhi}",
                "--salt-model", source["salt_model"], "--profile", source["profile"],
                *condition_args(source["conditions"]))
        by_pair = {row["pair_id"]: normalized_row(row) for row in read_tsv(folder / "panel_qc.tsv")}
        groups = collections.defaultdict(list)
        for rec in jsonl(folder / "panel_qc_records.jsonl"):
            groups[rec.get("target") or rec.get("primer_id")].append(rec)
        bindings = collections.defaultdict(list)
        for row in read_tsv(folder / "panel.bindings.tsv"):
            bindings[(row["pair_id"], row["role"])].append(row)
        products = collections.defaultdict(list)
        for row in read_tsv(folder / "panel.products.tsv"):
            products[row["pair_id"]].append(normalized_row(row))
        for pair in selected:
            t = self.pair_source(source, pair)["template"]
            pid = pair["id"]
            prods = products[pid]
            intended = [p for p in prods if p.get("orientation", "F-R") not in ("F-F", "R-R")
                        and p["chrom"] == t["chrom"]
                        and p["amplicon_start"] == t["offset"] + pair["start"]
                        and p["amplicon_end_exclusive"] == t["offset"] + pair["end"]]
            result["pairs"][pid] = {
                "summary": by_pair.get(pid), "intended": intended,
                "records": groups[pid], "primers": {
                    role: groups[f"{pid}_{role}"] for role in ("F", "R")},
                "bindings": {role: pack_bindings(bindings[pid, role], anchor) for role in ("F", "R")},
                "offtarget": sorted((p for p in prods if p not in intended),
                                    key=lambda p: p["total_mismatches"])[:50],
            }
        result["records"] = groups["__panel__"]
        stage("Python · yerel varyant kapsamı ve sıralama")
        for pair in selected:
            result["pairs"][pair["id"]]["variants"] = self.pair_variants(source, pair)
        if all("design" in pair for pair in selected):
            result["ranking"] = self.rerank(selected, result["pairs"], panel, folder)
        if len(selected) > 1:
            stage("Python · multipleks çapraz ürünleri ve dimer matrisi")
            mxdir = work / "multiplex"
            run_cli("multiplex_qc.py", "--prefix", folder / "panel", "--out", mxdir,
                    "--panel", panel, "--min-amplicon", min_amp, "--max-amplicon", max_amp,
                    "--na", source["conditions"]["na"], "--mg", source["conditions"]["mg"],
                    "--dntp", source["conditions"]["dntp"], "--primer-nm", source["conditions"]["pnm"],
                    "--fragment-bp", 166 if source["profile"] == "ctdna" else 0)
            result["multiplex"] = {
                "summary": read_tsv(mxdir / "multiplex_summary.tsv"),
                "dimer": [normalized_row(row) for row in read_tsv(mxdir / "dimer_matrix.tsv")],
                "cross": [normalized_row(row) for row in read_tsv(mxdir / "cross_products.tsv")][:200],
            }
        return result

    def rerank(self, selected, scan_pairs, panel, folder):
        (folder / "design_panel.tsv").write_bytes(panel.read_bytes())
        write_tsv(folder / "design_pairs.tsv",
                  ["left_primer", "right_primer", "product_length", "tm_left_c", "tm_right_c",
                   "delta_tm_c", "pareto_rank", "uncalibrated_score"],
                  [(p["F"], p["R"], p["design"]["product_length"],
                    p["design"]["left"]["tm_c"], p["design"]["right"]["tm_c"],
                    p["design"]["delta_tm"], p["design"]["pareto_rank"],
                    p["design"]["score"]) for p in selected])
        variant_rows = []
        for pair in selected:
            layer = scan_pairs[pair["id"]]["variants"]
            for rec in layer.get("records", []):
                variant_rows.append((pair["id"], rec["id"], rec["severity"], json.dumps(rec["value"])))
        write_tsv(folder / "variant_qc.tsv", ["pair_id", "check", "severity", "value"], variant_rows)
        run_cli("rerank.py", "--design", folder, "--qc", folder, "--variants", folder)
        return [normalized_row(row) for row in read_tsv(folder / "final_ranked.tsv")]

    def single(self, action, source, pair, data):
        source = self.pair_source(source, pair)
        if action == "qc":
            return self.pair_qc(source, pair)
        if action == "variants":
            return self.pair_variants(source, pair)
        t = source["template"]
        cond = condition_model(source["conditions"])
        seq = t["sequence"][pair["start"]:pair["end"]]
        if action in ("melt", "probe") and "N" in seq:
            raise ValueError("Amplikonda belirsiz N bazları bulundu. Bu model tam dizi gerektirir.")
        if action == "melt":
            if len(seq) > 1000:
                raise ValueError("Erime profili için ürün en fazla 1000 baz olmalı.")
            melt = amplicon.melt_profile(seq, cond)
            return {"summary": melt.as_dict(), "check": amplicon.qpcr02(melt),
                    "curve": list(zip(melt.temperatures, melt.helicity, melt.derivative))}
        if action == "probe":
            if len(seq) > 600:
                raise ValueError("Prob seçimi için ürün en fazla 600 baz olmalı.")
            s = dna(data.get("probe"), "Prob", 8, 80, optional=True)
            rejects = collections.Counter()
            if not s:
                candidates, rejected = probe.enumerate_probes(
                    seq, pair["F"], pair["R"], cond, salt_model=source["salt_model"])
                rejects.update(rule for rule, _ in rejected)
                if not candidates:
                    return {"status": "no_candidate", "rejections": dict(rejects)}
                s = candidates[0].sequence
            out = probe.evaluate_probe(pair["id"] + "_P", s, pair["F"], pair["R"],
                                       amplicon=seq, cond=cond, salt_model=source["salt_model"])
            return {"result": out, "rejections": dict(rejects)}
        if action == "isoform":
            if t["conversion"] != "primary":
                raise ValueError("İzoform kapsamı dönüşmemiş insan referansı gerektirir.")
            gene = data.get("gene", "")
            if not isinstance(gene, str) or not re.fullmatch(r"[A-Za-z0-9_.-]{1,40}", gene):
                raise ValueError("GENCODE gen simgesini girin.")
            exons = ANNOTATIONS / "transcript_exons.tsv"
            if not exons.is_file():
                raise ValueError("GENCODE ekson tablosu bulunamadı.")
            idx = load_exons_for_genes(str(exons), {gene})
            a, b = t["offset"] + pair["start"], t["offset"] + pair["end"]
            return {"check": isoform.check_rt03(gene, t["chrom"], (a, a + len(pair["F"])),
                    (b - len(pair["R"]), b), idx)}
        raise ValueError("Bilinmeyen analiz türü.")

    def oligo(self, action, data, work):
        conditions = reaction(data.get("conditions", {}))
        if action == "thermo":
            s = dna(data.get("sequence"), "Primer", 8, 80, modified=True)
            return {"result": engine.ThermoEngine(cond=condition_model(conditions)).evaluate("WEB", s)}
        mutant = dna(data.get("mutant"), "Mutant prob", 8, 80, modified=True)
        wild = dna(data.get("wt"), "WT prob", 8, 80, modified=True)
        primers = [dna(data.get(k), label, 8, 80, modified=True, optional=True)
                   for k, label in (("forward", "İleri primer"), ("reverse", "Geri primer"),
                                    ("blocker", "Blocker"))]
        if bool(primers[0]) != bool(primers[1]) or (primers[2] and not primers[0]):
            raise ValueError("Dimer denetimi için ileri ve geri primer birlikte girilmeli.")
        run_cli("blocker_qc.py", f"--mutant={mutant}", f"--wt={wild}",
                f"--forward={primers[0]}", f"--reverse={primers[1]}", f"--blocker={primers[2]}",
                "--out", work, *condition_args(conditions))
        return {"records": jsonl(work / "ct05_qc_records.jsonl")}


def specific_pairs(template, candidates, work):
    conversion = template["conversion"]
    folder = work / conversion
    offset, chrom = template["offset"], template["chrom"]
    expected = {p["id"]: (chrom, offset + p["start"], offset + p["end"]) for p in candidates}
    products, intended, exact = collections.Counter(), collections.Counter(), collections.Counter()
    with (folder / "panel.products.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["pair_id"]
            products[key] += 1
            locus = row["chrom"], int(row["amplicon_start"]), int(row["amplicon_end_exclusive"])
            if (locus == expected.get(key) and int(row["total_mismatches"]) == 0
                    and row.get("orientation", "F-R") not in ("F-F", "R-R")):
                intended[key] += 1
    with (folder / "panel.bindings.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if int(row["mismatches"]) == 0:
                exact[row["pair_id"], row["role"]] += 1
    other_products = set()
    if conversion in ("CT", "GA"):
        alternate = "GA" if conversion == "CT" else "CT"
        with (work / alternate / "panel.products.tsv").open(encoding="utf-8", newline="") as handle:
            other_products = {row["pair_id"] for row in csv.DictReader(handle, delimiter="\t")}
    return [
        pair for pair in candidates if products[pair["id"]] == intended[pair["id"]] == 1
        and exact[pair["id"], "F"] == exact[pair["id"], "R"] == 1
        and pair["id"] not in other_products
    ]


def copy_selected_rows(destination, key, suffix, batches):
    count = 0
    with destination.open("w", encoding="utf-8", newline="") as out:
        writer = None
        for folder, ids in batches:
            with (folder / key / f"panel.{suffix}.tsv").open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if writer is None:
                    writer = csv.DictWriter(out, reader.fieldnames, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                for row in reader:
                    if row["pair_id"] in ids:
                        writer.writerow(row)
                        count += 1
    return count


def variants_union(coverage):
    merged = []
    for chrom, start, end in sorted(coverage):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = chrom, merged[-1][1], max(merged[-1][2], end)
        else:
            merged.append((chrom, start, end))
    return merged


def template_data(chrom, offset, sequence, conversion):
    region = f"{chrom}:{offset + 1}-{offset + len(sequence)}"
    fasta = f">{region} GRCh38.p14 {conversion}\n{sequence}\n"
    return {"chrom": chrom, "offset": offset, "sequence": sequence, "conversion": conversion,
            "region": region, "fasta_sha256": hashlib.sha256(fasta.encode()).hexdigest()}


_CONTIG_ROLES = None


def spc08_extra_contigs(f, r, mode="anchored"):
    """SPC-08 tam kapsam: chrM + ALT + yamalar + iskeleler (C++ kaba kuvvet tarama).

    Ikili ya da dizi dosyasi yoksa None doner; cagiran chrM-only Python yoluna duser.
    """
    global _CONTIG_ROLES
    if not (CONTIG_BINARY.is_file() and os.access(CONTIG_BINARY, os.X_OK)
            and EXTRA_CONTIGS.is_file() and EXTRA_CONTIGS_TSV.is_file()):
        return None
    if _CONTIG_ROLES is None:
        _CONTIG_ROLES = {row["contig"]: row["role"] for row in read_tsv(EXTRA_CONTIGS_TSV)}
    import tempfile
    with tempfile.TemporaryDirectory(prefix="oligoscout-spc08-") as tmp:
        primers = Path(tmp) / "primers.tsv"
        primers.write_text("id\tseq\n" + "".join(f"{s}\t{s}\n" for s in dict.fromkeys((f, r))))
        out = Path(tmp) / "hits.tsv"
        run_command([str(CONTIG_BINARY), str(EXTRA_CONTIGS), str(primers), str(out), "3", "2", "--all"])
        sites = read_tsv(out)
    if mode == "sensitive":   # Primer-BLAST terminal rule, same as the genome search
        return seqpair.check_spc08_scan({"F": f, "R": r}, sites, _CONTIG_ROLES,
                                        max_mismatch=3, terminal=(5, 1))
    return seqpair.check_spc08_scan({"F": f, "R": r}, sites, _CONTIG_ROLES)


def search_mode(data, source=None):
    """anchored (exact 3' anchor, default) or sensitive (Primer-BLAST terminal rule)."""
    mode = (data or {}).get("search_mode") or (source or {}).get("search_mode") or "anchored"
    if mode not in ("anchored", "sensitive"):
        raise ValueError("Arama modu 'anchored' ya da 'sensitive' olmalı.")
    return mode


def pack_bindings(rows, anchor):
    hist = collections.Counter()
    sites = []
    for row in rows:
        mm = int(row["mismatches"])
        anchored = "X" not in row["mask"][-anchor:]
        hist[mm, anchored] += 1
        site = row["site_plus_strand"]
        if row["strand"] == "-":
            site = seqpair.reverse_complement_raw(site)
        sites.append({"chrom": row["chrom"], "start": int(row["start"]), "strand": row["strand"],
                      "mm": mm, "mask": row["mask"], "site": site, "anchored": anchored})
    sites.sort(key=lambda x: (x["mm"], not x["anchored"], x["chrom"], x["start"]))
    return {"n": len(rows), "top": sites[:40],
            "hist": [{"mm": k[0], "anchored": k[1], "n": n} for k, n in sorted(hist.items())]}
