
"""Design a multi-target panel with the existing GRCh38 and thermodynamics engines."""
from __future__ import annotations

import collections
import csv
import itertools
import time
from pathlib import Path

from local_reference import parse_region
from oligoscout_thermo import multiplex as mpx


def parse_targets(raw):
    if not isinstance(raw, list) or not 2 <= len(raw) <= 8:
        raise ValueError("Multipleks tasarımı için 2–8 hedef bölge girin.")
    targets, seen, names = [], set(), set()
    for index, item in enumerate(raw, 1):
        if not isinstance(item, dict):
            raise ValueError("Her hedef için bir ad ve genom bölgesi gerekli.")
        name = item.get("name", "").strip() if isinstance(item.get("name", ""), str) else ""
        if not name or len(name) > 60 or any(ord(c) < 32 for c in name):
            raise ValueError("Hedef adı 1–60 karakter olmalı.")
        chrom, start, end = parse_region(item.get("region"))
        region = f"{chrom}:{start + 1}-{end}"
        if region in seen or name.casefold() in names:
            raise ValueError("Hedef adları ve bölgeleri birbirinden farklı olmalı.")
        seen.add(region)
        names.add(name.casefold())
        targets.append({"id": f"T{index:02d}", "name": name, "region": region})
    return targets


def remap_ids(value, lookup):
    """Preserve all evidence while assigning candidate IDs unique to each target."""
    def rename(text):
        for old, new in lookup.items():
            if text == old or text.startswith(old + "_"):
                return new + text[len(old):]
        return text
    if isinstance(value, dict):
        return {rename(key): remap_ids(item, lookup) for key, item in value.items()}
    if isinstance(value, list):
        return [remap_ids(item, lookup) for item in value]
    return rename(value) if isinstance(value, str) else value


def edge_key(a, b):
    return tuple(sorted((a, b)))


def select_pool(groups, blocked, limit=100000):
    """Choose at most one candidate per target; maximize target coverage, with a node cap."""
    ordered = sorted(groups, key=lambda group: (len(group), group[0]["target_id"] if group else ""))
    best, visits, capped = [], 0, False

    def visit(index, chosen):
        nonlocal best, visits, capped
        if len(best) == len(groups):
            return
        if visits >= limit:
            capped = True
            return
        visits += 1
        if len(chosen) > len(best):
            best = list(chosen)
        if index == len(ordered) or len(chosen) + len(ordered) - index <= len(best):
            return
        for candidate in ordered[index]:
            if edge_key(candidate["id"], candidate["id"]) in blocked:
                continue
            if all(edge_key(candidate["id"], other["id"]) not in blocked for other in chosen):
                visit(index + 1, chosen + [candidate])
        visit(index + 1, chosen)

    visit(0, [])
    return best, {"visited": visits, "limit_reached": capped}


def assess_candidates(candidates, binding_path, cond, settings, stage):
    owners = collections.defaultdict(set)
    lookup = {p["id"]: p for p in candidates}
    for p in candidates:
        for role in ("F", "R"):
            owners[p[role]].add(p["id"])
    stage("Multipleks · tüm primerlerin dimer etkileşimleri")
    dimers = mpx.dimer_matrix({seq: seq for seq in owners}, cond)
    dg = {edge_key(r["primer_a"], r["primer_b"]): r for r in dimers}
    cross_counts, cross_examples = collections.Counter(), {}
    if binding_path:
        sites, _ = mpx.load_bindings(str(binding_path))
        stage("Multipleks · hedefler arası genomik çapraz ürünler")
        for product in mpx.iter_products(sites, set(), settings["min_product"], settings["max_product"]):
            sequence_pair = frozenset((product.left_id, product.right_id))
            product_edges = set()
            for a in owners[product.left_id]:
                for b in owners[product.right_id]:
                    pa, pb = lookup[a], lookup[b]
                    if a != b and pa["target_id"] == pb["target_id"]:
                        continue
                    if sequence_pair in (frozenset((pa["F"], pa["R"])), frozenset((pb["F"], pb["R"]))):
                        continue
                    product_edges.add(edge_key(a, b))
            for key in product_edges:
                cross_counts[key] += 1
                cross_examples.setdefault(key, {
                    "chrom": product.chrom, "start": product.start + 1, "end": product.end,
                    "length": product.length, "mismatches": product.mismatches,
                })
    blocked, conflicts = set(), []
    for candidate in candidates:
        own = [dg[edge_key(candidate[x], candidate[y])] for x, y in (("F", "F"), ("F", "R"), ("R", "R"))]
        worst_any = min((row["dG37_any"] for row in own if row["dG37_any"] is not None), default=None)
        worst_3p = min((row["dG37_3prime"] for row in own if row["dG37_3prime"] is not None), default=None)
        spread = abs(candidate["design"]["left"]["tm_c"] - candidate["design"]["right"]["tm_c"])
        reasons = []
        own_key = edge_key(candidate["id"], candidate["id"])
        if cross_counts[own_key]:
            reasons.append("Tek primerle tasarlanmamış genomik ürün")
        if worst_any is not None and worst_any < settings["dimer_any"]:
            reasons.append("Çiftin kendi dimer ΔG değeri")
        if worst_3p is not None and worst_3p < settings["dimer_3p"]:
            reasons.append("Çiftin kendi 3′ dimer ΔG değeri")
        if spread > settings["tm_spread"]:
            reasons.append("Çiftin kendi Tm aralığı")
        if reasons:
            blocked.add(edge_key(candidate["id"], candidate["id"]))
            conflicts.append({"a": candidate["id"], "b": candidate["id"], "reasons": reasons,
                              "cross_products": cross_counts[own_key], "cross_example": cross_examples.get(own_key), "dG37_any": worst_any,
                              "dG37_3prime": worst_3p, "tm_spread": round(spread, 3), "product_gap": None})
    for a, b in itertools.combinations(candidates, 2):
        if a["target_id"] == b["target_id"]:
            continue
        key = edge_key(a["id"], b["id"])
        combinations = [dg[edge_key(a[x], b[y])] for x in ("F", "R") for y in ("F", "R")]
        worst_any = min((r["dG37_any"] for r in combinations if r["dG37_any"] is not None), default=None)
        worst_3p = min((r["dG37_3prime"] for r in combinations if r["dG37_3prime"] is not None), default=None)
        tms = [p["design"][role]["tm_c"] for p in (a, b) for role in ("left", "right")]
        spread = max(tms) - min(tms)
        gap = abs(a["design"]["product_length"] - b["design"]["product_length"])
        reasons = []
        if cross_counts[key]:
            reasons.append("Genomik çapraz ürün")
        if worst_any is not None and worst_any < settings["dimer_any"]:
            reasons.append("Dimer ΔG")
        if worst_3p is not None and worst_3p < settings["dimer_3p"]:
            reasons.append("3′ dimer ΔG")
        if spread > settings["tm_spread"]:
            reasons.append("Tm aralığı")
        if gap < settings["product_gap"]:
            reasons.append("Ürün boyu ayrımı")
        if reasons:
            blocked.add(key)
            conflicts.append({
                "a": a["id"], "b": b["id"], "reasons": reasons,
                "cross_products": cross_counts[key], "cross_example": cross_examples.get(key),
                "dG37_any": worst_any, "dG37_3prime": worst_3p,
                "tm_spread": round(spread, 3), "product_gap": gap,
            })
    return blocked, conflicts, dimers


def run(motor, data, stage, work):
    from local_engine import condition_model, finite
    targets = parse_targets(data.get("targets"))
    profile = data.get("profile", "sanger")
    if profile == "bisulfite":
        raise ValueError("Bu multipleks akışı dönüşmemiş GRCh38 hedefleri içindir.")
    count = finite(data.get("candidates", 3), "Hedef başına aday", 1, 4, integer=True)
    raw = data.get("pool", {})
    if not isinstance(raw, dict):
        raise ValueError("Havuz ayarları geçersiz.")
    settings = {
        "tm_spread": finite(raw.get("tm_spread", 3), "Panel Tm aralığı", .1, 15),
        "dimer_any": finite(raw.get("dimer_any", -6), "Panel dimer ΔG", -30, 0),
        "dimer_3p": finite(raw.get("dimer_3p", -5), "Panel 3′ dimer ΔG", -30, 0),
        "product_gap": finite(raw.get("product_gap", 10), "Ürün boyu farkı", 0, 500, integer=True),
    }
    # Validate every target before starting the expensive C++ searches.
    for target in targets:
        motor.read_region(target["region"])
    start = time.monotonic()
    outputs, sources, candidates, binding_rows, binding_fields = [], {}, [], [], None
    for index, target in enumerate(targets, 1):
        directory = work / target["id"]
        directory.mkdir()
        label = f"{index}/{len(targets)} · {target['name']}"
        payload = {**data, "tool": "design", "region": target["region"], "profile": profile, "top": count}
        payload.pop("reference_sequence", None)
        result = motor.design(payload, lambda text: stage(f"{label} · {text}"), directory)
        names = {p["id"]: target["id"] + "_" + p["id"] for p in result["pairs"]}
        result = remap_ids(result, names)
        for p in result["pairs"]:
            p["target_id"] = target["id"]
        outputs.append({**target, "result": result})
        sources[target["id"]] = result
        candidates.extend(result["pairs"])
        raw_bindings = directory / "selected" / "primary" / "panel.bindings.tsv"
        if result["pairs"]:
            with raw_bindings.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                binding_fields = reader.fieldnames
                for row in reader:
                    row["pair_id"] = names[row["pair_id"]]
                    binding_rows.append(row)
    source = outputs[0]["result"]
    settings["min_product"] = min(50, source["constraints"]["product_min"])
    settings["max_product"] = max(3000, source["constraints"]["product_max"])
    binding_path = None
    if binding_fields:
        binding_path = work / "panel.bindings.tsv"
        with binding_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=binding_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(binding_rows)
    blocked, conflicts, dimers = assess_candidates(
        candidates, binding_path, condition_model(source["conditions"]), settings, stage)
    stage("Multipleks · hedef başına bir çift ile uyumlu panel seçimi")
    groups = [target["result"]["pairs"] for target in outputs]
    selected, search = select_pool(groups, blocked)
    selected.sort(key=lambda p: p["target_id"])
    selected_ids = {p["id"] for p in selected}
    selected_sequences = {p[role] for p in selected for role in ("F", "R")}
    labels = collections.defaultdict(list)
    for p in selected:
        for role in ("F", "R"):
            labels[p[role]].append(p["target_id"] + "." + role)
    selected_dimers = [{**row, "primer_a": " / ".join(labels[row["primer_a"]]),
                       "primer_b": " / ".join(labels[row["primer_b"]])}
                      for row in dimers
                      if row["primer_a"] in selected_sequences and row["primer_b"] in selected_sequences]
    all_tms = [p["design"][role]["tm_c"] for p in selected for role in ("left", "right")]
    sizes = sorted(p["design"]["product_length"] for p in selected)
    missing = [t["id"] for t in targets if not any(p["target_id"] == t["id"] for p in selected)]
    return {
        "reference": "GRCh38.p14", "profile": source["profile"], "conditions": source["conditions"],
        "salt_model": source["salt_model"], "constraints": source["constraints"],
        "template": source["template"], "target_sources": sources, "targets": outputs,
        "pairs": candidates, "selected": [p["id"] for p in selected], "missing": missing,
        "complete": len(selected_ids) == len(targets), "conflicts": conflicts,
        "dimer": selected_dimers, "settings": settings, "search": search,
        "summary": {
            "targets": len(targets), "candidates": len(candidates), "selected": len(selected),
            "tm_spread": round(max(all_tms) - min(all_tms), 3) if all_tms else None,
            "product_gap": min((b - a for a, b in zip(sizes, sizes[1:])), default=None),
            "cross_products": 0 if len(selected) >= 2 else None,
        }, "elapsed_ms": round((time.monotonic() - start) * 1000),
    }
