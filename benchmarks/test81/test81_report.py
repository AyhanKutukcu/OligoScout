#!/usr/bin/env python3
"""Create compact publication-oriented Test #81 TSV, Markdown, and SVG reports."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from pathlib import Path


METHODS = [
    ("oracle", "Independent full-scan oracle"),
    ("oligoscout", "OligoScout"),
    ("bowtie1", "Bowtie 1 + contract filter"),
    ("bwa_aln", "BWA-aln + contract filter"),
    ("blastn_e1100", "BLASTN-short e=1100 + contract filter"),
]


def elapsed_seconds(text: str) -> float:
    parts = text.split(":")
    if len(parts) == 2:
        return float(parts[0]) * 60 + float(parts[1])
    if len(parts) == 3:
        return float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
    return float(text)


def parse_time(path: Path) -> tuple[float, int]:
    text = path.read_text(encoding="utf-8", errors="replace")
    elapsed_match = re.search(
        r"Elapsed \(wall clock\) time .*?:\s*([0-9:.]+)\s*$", text, re.MULTILINE
    )
    rss_match = re.search(
        r"Maximum resident set size \(kbytes\):\s*([0-9]+)", text
    )
    if not elapsed_match or not rss_match:
        raise RuntimeError(f"Cannot parse GNU time file: {path}")
    return elapsed_seconds(elapsed_match.group(1)), int(rss_match.group(1))


def method_time(directory: Path) -> tuple[float, int, str]:
    direct = directory / "times/search.time.txt"
    if direct.exists():
        elapsed, rss = parse_time(direct)
        return elapsed, rss, "single measured search process"
    chunks = sorted((directory / "times/blast_chunks").glob("*.time.txt"))
    values = [parse_time(path) for path in chunks]
    if not values:
        return 0.0, 0, "not available"
    return sum(value[0] for value in values), max(value[1] for value in values), \
        f"sum of {len(values)} sequential query chunks; RSS is chunk maximum"


def pipeline_time(directory: Path) -> tuple[float, int]:
    paths = [
        directory / "times/normalize.time.txt",
        directory / "times/sort_bindings.time.txt",
        directory / "times/assemble.time.txt",
        directory / "times/sort_products.time.txt",
    ]
    search_path = directory / "times/search.time.txt"
    if search_path.exists():
        paths.insert(0, search_path)
    else:
        paths = sorted((directory / "times/blast_chunks").glob("*.time.txt")) + paths
    values = [parse_time(path) for path in paths if path.exists()]
    return sum(value[0] for value in values), max((value[1] for value in values), default=0)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def accuracy_svg(rows: list[dict], output: Path) -> None:
    width, height = 920, 430
    left, top, chart_width, chart_height = 250, 45, 620, 310
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:serif;fill:#000}.axis{stroke:#000;stroke-width:1}.bar{fill:#555}.grid{stroke:#bbb;stroke-width:.7;stroke-dasharray:3 3}</style>',
        '<text x="460" y="25" text-anchor="middle" font-size="17">Product-set F1 against independent oracle</text>',
    ]
    for tick in (0.0, 0.25, 0.5, 0.75, 1.0):
        x = left + tick * chart_width
        lines.append(f'<line class="grid" x1="{x}" y1="{top}" x2="{x}" y2="{top + chart_height}"/>')
        lines.append(f'<text x="{x}" y="{top + chart_height + 22}" text-anchor="middle" font-size="12">{tick:.2f}</text>')
    bar_height = 38
    gap = 22
    for index, row in enumerate(rows):
        y = top + 15 + index * (bar_height + gap)
        value = float(row["product_f1"])
        lines.append(f'<text x="{left - 10}" y="{y + 25}" text-anchor="end" font-size="13">{html.escape(row["label"])}</text>')
        lines.append(f'<rect class="bar" x="{left}" y="{y}" width="{value * chart_width:.2f}" height="{bar_height}"/>')
        lines.append(f'<text x="{left + value * chart_width - 6}" y="{y + 25}" text-anchor="end" fill="white" font-size="13">{value:.6f}</text>')
    lines.append(f'<line class="axis" x1="{left}" y1="{top + chart_height}" x2="{left + chart_width}" y2="{top + chart_height}"/>')
    lines.append('</svg>')
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def runtime_svg(rows: list[dict], output: Path) -> None:
    positive = [float(row["search_wall_seconds"]) for row in rows if float(row["search_wall_seconds"]) > 0]
    maximum = max(positive, default=1.0)
    width, height = 920, 430
    left, top, chart_width, chart_height = 250, 45, 620, 310
    scale = math.log10(maximum + 1)
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:serif;fill:#000}.axis{stroke:#000;stroke-width:1}.bar{fill:#333}</style>',
        '<text x="460" y="25" text-anchor="middle" font-size="17">Measured search wall time (log-scaled bar length)</text>',
    ]
    bar_height, gap = 38, 22
    for index, row in enumerate(rows):
        y = top + 15 + index * (bar_height + gap)
        seconds = float(row["search_wall_seconds"])
        bar_width = (math.log10(seconds + 1) / scale) * chart_width if seconds > 0 else 0
        lines.append(f'<text x="{left - 10}" y="{y + 25}" text-anchor="end" font-size="13">{html.escape(row["label"])}</text>')
        lines.append(f'<rect class="bar" x="{left}" y="{y}" width="{bar_width:.2f}" height="{bar_height}"/>')
        lines.append(f'<text x="{left + 8}" y="{y + 25}" fill="white" font-size="13">{seconds:.2f} s</text>')
    lines.append(f'<line class="axis" x1="{left}" y1="{top + chart_height}" x2="{left + chart_width}" y2="{top + chart_height}"/>')
    lines.append('</svg>')
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--panel-size", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    oracle_directory = args.stage_dir / "oracle"
    for method, label in METHODS:
        directory = args.stage_dir / method
        comparison = directory / "comparison"
        binding = read_json(comparison / "bindings.summary.json")
        product = read_json(comparison / "products.summary.json")
        wall, rss, scope = method_time(directory)
        pipeline_wall, pipeline_rss = pipeline_time(directory)
        rows.append({
            "method": method,
            "label": label,
            "panel_pairs": args.panel_size,
            "binding_truth": binding["truth_count"],
            "binding_predictions": binding["prediction_count"],
            "binding_precision": binding["micro_precision"],
            "binding_recall": binding["micro_recall"],
            "binding_f1": binding["micro_f1"],
            "binding_jaccard": binding["jaccard"],
            "product_truth": product["truth_count"],
            "product_predictions": product["prediction_count"],
            "product_precision": product["micro_precision"],
            "product_recall": product["micro_recall"],
            "product_f1": product["micro_f1"],
            "product_jaccard": product["jaccard"],
            "intended_recall": product.get("intended_product_recall", ""),
            "macro_pair_f1": product["macro_pair_f1"],
            "macro_pair_f1_ci_low": product["macro_pair_f1_bootstrap_95_low"],
            "macro_pair_f1_ci_high": product["macro_pair_f1_bootstrap_95_high"],
            "search_wall_seconds": wall,
            "search_max_rss_kib": rss,
            "timing_scope": scope,
            "end_to_end_wall_seconds": pipeline_wall,
            "end_to_end_max_stage_rss_kib": pipeline_rss,
        })
    tsv = args.output_dir / "test81_results.tsv"
    with tsv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    markdown = args.output_dir / "test81_report.md"
    with markdown.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"# Test #81 publication benchmark — {args.panel_size} primer pairs\n\n")
        stream.write("Contract: 24-nt primers; exact 12-nt 3′ anchor; at most 3 mismatches per primer; 50–3000 bp products; GRCh38.p14 primary assembly.\n\n")
        stream.write("| Method | Binding P/R/F1 | Product P/R/F1 | Intended recall | Search wall | End-to-end wall | Max stage RSS |\n")
        stream.write("|---|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            stream.write(
                f"| {row['label']} | {float(row['binding_precision']):.6f} / {float(row['binding_recall']):.6f} / {float(row['binding_f1']):.6f} "
                f"| {float(row['product_precision']):.6f} / {float(row['product_recall']):.6f} / {float(row['product_f1']):.6f} "
                f"| {float(row['intended_recall']):.6f} | {float(row['search_wall_seconds']):.2f} s "
                f"| {float(row['end_to_end_wall_seconds']):.2f} s "
                f"| {int(row['end_to_end_max_stage_rss_kib']):,} KiB |\n"
            )
        stream.write("\nBowtie 1, BWA-aln, and BLASTN-short are general alignment/search baselines. Their raw hits were independently revalidated against the reference and filtered to the same primer-binding contract. They are not described as primer-native tools.\n\n")
        stream.write("The independent oracle performs a direct full-reference scan. Product-level uncertainty is reported with 10,000 fixed-seed bootstrap resamples using the primer pair—not individual hits—as the resampling unit.\n\n")
        stream.write("Timing caveat: large outputs and temporary sort files reside on the Windows E: drive through WSL DrvFS. These wall times describe this machine and storage path; accuracy/set metrics are unaffected. Search wall time is reported separately from an end-to-end sum of all sequential measured stages. BLAST search wall time is the sum of sequential query chunks and its RSS is the maximum observed chunk RSS; adaptive smaller chunks may be used for resource-heavy queries and are retained in the audit trail. OligoScout's native search stage also materializes binding and product tables; external tools first materialize raw alignments and then undergo contract normalization and product assembly. Interpret both the search-phase and end-to-end columns with these scope differences in view.\n")
    accuracy_svg(rows, args.output_dir / "test81_product_f1.svg")
    runtime_svg(rows, args.output_dir / "test81_search_runtime.svg")
    print(f"RESULT_TABLE\t{tsv}")
    print(f"REPORT\t{markdown}")
    print("TEST81_REPORT_COMPLETE\tYES")


if __name__ == "__main__":
    main()
