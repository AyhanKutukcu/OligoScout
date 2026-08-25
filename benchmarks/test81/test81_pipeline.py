#!/usr/bin/env python3
"""Normalization, product assembly, and exact comparison for Test #81."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import mmap
import random
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


BINDING_HEADER = [
    "pair_id", "primer_identity", "query_id", "chromosome", "strand",
    "binding_start", "binding_end_exclusive", "primer_length",
    "mismatch_count", "mismatch_mask", "anchor_length", "anchor_mismatches",
    "last3_mismatches", "last5_mismatches", "primer_sequence",
    "binding_sequence",
]
PRODUCT_HEADER = [
    "pair_id", "chromosome", "left_primer", "right_primer",
    "left_position", "right_position", "amplicon_start",
    "amplicon_end_exclusive", "amplicon_length", "left_mismatches",
    "right_mismatches", "total_mismatches", "left_anchor_mismatches",
    "right_anchor_mismatches", "left_last3_mismatches",
    "right_last3_mismatches", "left_mismatch_mask", "right_mismatch_mask",
]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


class IndexedFasta:
    def __init__(self, fasta: Path, fai: Path):
        self.handle = fasta.open("rb")
        self.data = mmap.mmap(self.handle.fileno(), 0, access=mmap.ACCESS_READ)
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with fai.open(encoding="utf-8") as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                self.index[fields[0]] = tuple(map(int, fields[1:5]))

    def fetch(self, name: str, start: int, end: int) -> str:
        if name not in self.index:
            raise KeyError(f"Reference sequence absent from FAI: {name}")
        length, offset, line_bases, line_bytes = self.index[name]
        if not (0 <= start <= end <= length):
            raise ValueError(f"Reference interval out of range: {name}:{start}-{end}")
        parts: list[bytes] = []
        cursor = start
        while cursor < end:
            within = cursor % line_bases
            take = min(end - cursor, line_bases - within)
            byte_offset = offset + (cursor // line_bases) * line_bytes + within
            parts.append(self.data[byte_offset:byte_offset + take])
            cursor += take
        return b"".join(parts).decode("ascii").upper()

    def close(self) -> None:
        self.data.close()
        self.handle.close()


@dataclass(frozen=True)
class Query:
    pair_id: int
    identity: str
    query_id: str
    sequence: str


def load_panel(path: Path) -> tuple[list[dict[str, str]], dict[str, Query]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise RuntimeError("Panel is empty")
    queries: dict[str, Query] = {}
    for row in rows:
        pair_id = int(row["pair_id"])
        for identity, id_field, sequence_field in (
            ("P1", "primer1_id", "primer1"),
            ("P2", "primer2_id", "primer2"),
        ):
            sequence = row[sequence_field].upper()
            if not sequence or set(sequence) - {"A", "C", "G", "T"}:
                raise RuntimeError(f"Invalid primer sequence: {row[id_field]}")
            query = Query(pair_id, identity, row[id_field], sequence)
            if query.query_id in queries:
                raise RuntimeError(f"Duplicate query ID: {query.query_id}")
            queries[query.query_id] = query
    return rows, queries


def mismatch_details(primer: str, binding: str) -> tuple[int, int]:
    mask = 0
    mismatches = 0
    for index, (expected, observed) in enumerate(zip(primer, binding)):
        if expected != observed:
            mismatches += 1
            mask |= 1 << index
    return mismatches, mask


def region_count(mask: int, length: int, region: int) -> int:
    begin = max(0, length - region)
    return sum((mask >> index) & 1 for index in range(begin, length))


def validated_binding(
    query: Query,
    chromosome: str,
    start: int,
    strand: str,
    reference: IndexedFasta,
    anchor: int,
    max_mismatches: int,
) -> list[str] | None:
    if start < 0:
        return None
    try:
        genomic = reference.fetch(chromosome, start, start + len(query.sequence))
    except (KeyError, ValueError):
        return None
    if set(genomic) - {"A", "C", "G", "T"}:
        return None
    binding = genomic if strand == "plus" else reverse_complement(genomic)
    mismatches, mask = mismatch_details(query.sequence, binding)
    anchor_mismatches = region_count(mask, len(query.sequence), anchor)
    if mismatches > max_mismatches or anchor_mismatches != 0:
        return None
    return [
        str(query.pair_id), query.identity, query.query_id, chromosome, strand,
        str(start), str(start + len(query.sequence)), str(len(query.sequence)),
        str(mismatches), f"0x{mask:016x}", str(anchor), str(anchor_mismatches),
        str(region_count(mask, len(query.sequence), 3)),
        str(region_count(mask, len(query.sequence), 5)), query.sequence, binding,
    ]


def command_make_queries(args: argparse.Namespace) -> None:
    _rows, queries = load_panel(args.panel)
    with args.fasta.open("w", encoding="ascii", newline="\n") as fasta, \
            args.fastq.open("w", encoding="ascii", newline="\n") as fastq:
        for query in queries.values():
            fasta.write(f">{query.query_id}\n{query.sequence}\n")
            fastq.write(
                f"@{query.query_id}\n{query.sequence}\n+\n{'I' * len(query.sequence)}\n"
            )
    print(f"QUERY_COUNT\t{len(queries)}")


def command_split_fasta(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    with args.input.open(encoding="ascii") as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence)))
                name = line
                sequence = []
            elif line:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence)))
    for chunk_index in range(0, len(records), args.records_per_chunk):
        chunk = records[chunk_index:chunk_index + args.records_per_chunk]
        path = args.output_dir / f"queries_{chunk_index // args.records_per_chunk:04d}.fa"
        with path.open("w", encoding="ascii", newline="\n") as output:
            for header, bases in chunk:
                output.write(f"{header}\n{bases}\n")
    print(f"FASTA_RECORDS\t{len(records)}")
    print(f"FASTA_CHUNKS\t{(len(records) + args.records_per_chunk - 1) // args.records_per_chunk}")


def sam_candidates(fields: list[str], include_xa: bool) -> Iterable[tuple[str, int, str]]:
    flag = int(fields[1])
    if not flag & 4 and fields[2] != "*":
        yield fields[2], int(fields[3]) - 1, "minus" if flag & 16 else "plus"
    if not include_xa:
        return
    for field in fields[11:]:
        if not field.startswith("XA:Z:"):
            continue
        for entry in field[5:].split(";"):
            if not entry:
                continue
            chromosome, signed_position, _cigar, _nm = entry.split(",")
            yield chromosome, abs(int(signed_position)) - 1, \
                "plus" if signed_position.startswith("+") else "minus"


def command_normalize_sam(args: argparse.Namespace) -> None:
    _rows, queries = load_panel(args.panel)
    reference = IndexedFasta(args.reference, args.fai)
    accepted = 0
    examined = 0
    try:
        with open_text(args.input) as source, args.output.open(
            "w", encoding="utf-8", newline=""
        ) as target:
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(BINDING_HEADER)
            for line in source:
                if not line or line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11 or fields[0] not in queries:
                    continue
                query = queries[fields[0]]
                for chromosome, start, strand in sam_candidates(fields, args.include_xa):
                    examined += 1
                    row = validated_binding(
                        query, chromosome, start, strand, reference,
                        args.anchor, args.max_mismatches,
                    )
                    if row is not None:
                        writer.writerow(row)
                        accepted += 1
    finally:
        reference.close()
    print(f"CANDIDATES_EXAMINED\t{examined}")
    print(f"CONTRACT_BINDINGS_WRITTEN\t{accepted}")


def command_normalize_blast(args: argparse.Namespace) -> None:
    _rows, queries = load_panel(args.panel)
    reference = IndexedFasta(args.reference, args.fai)
    examined = 0
    accepted = 0
    try:
        with open_text(args.input) as source, args.output.open(
            "w", encoding="utf-8", newline=""
        ) as target:
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(BINDING_HEADER)
            for raw in source:
                if not raw.strip():
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 15 or fields[0] not in queries:
                    continue
                qstart, qend = int(fields[2]), int(fields[3])
                sstart, send = int(fields[4]), int(fields[5])
                gaps = int(fields[9])
                if gaps:
                    continue
                query = queries[fields[0]]
                strand = "plus" if fields[6] == "plus" else "minus"
                if strand == "plus":
                    start = min(sstart, send) - 1 - (qstart - 1)
                else:
                    start = min(sstart, send) - 1 - (len(query.sequence) - qend)
                examined += 1
                row = validated_binding(
                    query, fields[1], start, strand, reference,
                    args.anchor, args.max_mismatches,
                )
                if row is not None:
                    writer.writerow(row)
                    accepted += 1
    finally:
        reference.close()
    print(f"CANDIDATES_EXAMINED\t{examined}")
    print(f"CONTRACT_BINDINGS_WRITTEN\t{accepted}")


@dataclass(frozen=True)
class Binding:
    identity: str
    start: int
    length: int
    mismatches: int
    anchor_mismatches: int
    last3_mismatches: int
    mask: str


def emit_products(
    writer: csv.writer,
    pair_id: int,
    chromosome: str,
    left: list[Binding],
    right: list[Binding],
    min_amplicon: int,
    max_amplicon: int,
) -> int:
    right.sort(key=lambda hit: hit.start)
    starts = [hit.start for hit in right]
    count = 0
    for left_hit in sorted(left, key=lambda hit: hit.start):
        minimum_right = left_hit.start + min_amplicon - (right[0].length if right else 0)
        maximum_right = left_hit.start + max_amplicon - (right[0].length if right else 0)
        begin = bisect.bisect_left(starts, minimum_right)
        end = bisect.bisect_right(starts, maximum_right)
        for right_hit in right[begin:end]:
            amplicon_end = right_hit.start + right_hit.length
            amplicon_length = amplicon_end - left_hit.start
            if not min_amplicon <= amplicon_length <= max_amplicon:
                continue
            writer.writerow([
                pair_id, chromosome, left_hit.identity, right_hit.identity,
                left_hit.start, right_hit.start, left_hit.start, amplicon_end,
                amplicon_length, left_hit.mismatches, right_hit.mismatches,
                left_hit.mismatches + right_hit.mismatches,
                left_hit.anchor_mismatches, right_hit.anchor_mismatches,
                left_hit.last3_mismatches, right_hit.last3_mismatches,
                left_hit.mask, right_hit.mask,
            ])
            count += 1
    return count


def command_assemble(args: argparse.Namespace) -> None:
    product_count = 0
    with args.bindings.open(encoding="utf-8", newline="") as source, \
            args.output.open("w", encoding="utf-8", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(PRODUCT_HEADER)
        current_key: tuple[int, str] | None = None
        group: dict[tuple[str, str], list[Binding]] = defaultdict(list)

        def flush() -> int:
            if current_key is None:
                return 0
            pair_id, chromosome = current_key
            total = 0
            total += emit_products(
                writer, pair_id, chromosome, group[("P1", "plus")],
                group[("P2", "minus")], args.min_amplicon, args.max_amplicon,
            )
            total += emit_products(
                writer, pair_id, chromosome, group[("P2", "plus")],
                group[("P1", "minus")], args.min_amplicon, args.max_amplicon,
            )
            return total

        for row in reader:
            key = (int(row["pair_id"]), row["chromosome"])
            if current_key is not None and key != current_key:
                product_count += flush()
                group = defaultdict(list)
            current_key = key
            group[(row["primer_identity"], row["strand"])].append(Binding(
                row["primer_identity"], int(row["binding_start"]),
                int(row["primer_length"]), int(row["mismatch_count"]),
                int(row["anchor_mismatches"]), int(row["last3_mismatches"]),
                row["mismatch_mask"],
            ))
        product_count += flush()
    print(f"PRODUCT_COUNT\t{product_count}")


def row_key(kind: str, row: dict[str, str]) -> tuple:
    if kind == "bindings":
        return (
            int(row["pair_id"]), row["chromosome"], row["primer_identity"],
            row["strand"], int(row["binding_start"]),
            int(row["binding_end_exclusive"]), int(row["mismatch_count"]),
            row["mismatch_mask"],
        )
    return (
        int(row["pair_id"]), row["chromosome"], int(row["amplicon_start"]),
        int(row["amplicon_end_exclusive"]), row["left_primer"],
        row["right_primer"], int(row["left_mismatches"]),
        int(row["right_mismatches"]), row["left_mismatch_mask"],
        row["right_mismatch_mask"],
    )


def key_iterator(path: Path, kind: str) -> Iterator[tuple]:
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            yield row_key(kind, row)


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def command_compare(args: argparse.Namespace) -> None:
    panel, _queries = load_panel(args.panel)
    pair_ids = [int(row["pair_id"]) for row in panel]
    truth_counts = defaultdict(int)
    prediction_counts = defaultdict(int)
    intersection_counts = defaultdict(int)
    truth = iter(key_iterator(args.truth, args.kind))
    prediction = iter(key_iterator(args.prediction, args.kind))
    left = next(truth, None)
    right = next(prediction, None)
    while left is not None or right is not None:
        if right is None or (left is not None and left < right):
            truth_counts[left[0]] += 1
            left = next(truth, None)
        elif left is None or right < left:
            prediction_counts[right[0]] += 1
            right = next(prediction, None)
        else:
            truth_counts[left[0]] += 1
            prediction_counts[right[0]] += 1
            intersection_counts[left[0]] += 1
            left = next(truth, None)
            right = next(prediction, None)
    per_pair: list[dict[str, float | int]] = []
    for pair_id in pair_ids:
        expected = truth_counts[pair_id]
        observed = prediction_counts[pair_id]
        shared = intersection_counts[pair_id]
        precision = shared / observed if observed else (1.0 if expected == 0 else 0.0)
        recall = shared / expected if expected else (1.0 if observed == 0 else 0.0)
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        per_pair.append({
            "pair_id": pair_id, "truth": expected, "prediction": observed,
            "intersection": shared, "precision": precision,
            "recall": recall, "f1": f1,
        })
    truth_total = sum(truth_counts.values())
    prediction_total = sum(prediction_counts.values())
    intersection_total = sum(intersection_counts.values())
    precision = intersection_total / prediction_total if prediction_total else 0.0
    recall = intersection_total / truth_total if truth_total else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    union = truth_total + prediction_total - intersection_total
    jaccard = intersection_total / union if union else 1.0
    pair_f1 = [float(row["f1"]) for row in per_pair]
    rng = random.Random(810081)
    bootstrap = [
        statistics.fmean(rng.choice(pair_f1) for _ in pair_f1)
        for _ in range(args.bootstrap_replicates)
    ]
    summary = {
        "kind": args.kind,
        "panel_pairs": len(pair_ids),
        "truth_count": truth_total,
        "prediction_count": prediction_total,
        "intersection_count": intersection_total,
        "truth_only_count": truth_total - intersection_total,
        "prediction_only_count": prediction_total - intersection_total,
        "micro_precision": precision,
        "micro_recall": recall,
        "micro_f1": f1,
        "jaccard": jaccard,
        "macro_pair_f1": statistics.fmean(pair_f1),
        "macro_pair_f1_bootstrap_95_low": quantile(bootstrap, 0.025),
        "macro_pair_f1_bootstrap_95_high": quantile(bootstrap, 0.975),
        "bootstrap_replicates": args.bootstrap_replicates,
        "bootstrap_unit": "primer_pair",
        "exact_set_equality": truth_total == prediction_total == intersection_total,
    }
    if args.kind == "products":
        expected = {
            int(row["pair_id"]): (
                row["expected_chromosome"],
                int(row["expected_amplicon_start"]),
                int(row["expected_amplicon_end_exclusive"]),
            )
            for row in panel
        }
        recovered: set[int] = set()
        with args.prediction.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                pair_id = int(row["pair_id"])
                target = expected.get(pair_id)
                if target == (
                    row["chromosome"], int(row["amplicon_start"]),
                    int(row["amplicon_end_exclusive"]),
                ):
                    recovered.add(pair_id)
        summary["intended_products_recovered"] = len(recovered)
        summary["intended_products_total"] = len(expected)
        summary["intended_product_recall"] = len(recovered) / len(expected)
    with args.per_pair.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(per_pair[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(per_pair)
    args.json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    with args.summary.open("w", encoding="utf-8", newline="\n") as stream:
        for key, value in summary.items():
            stream.write(f"{key}\t{value}\n")
    for key, value in summary.items():
        print(f"{key.upper()}\t{value}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    make = subparsers.add_parser("make-queries")
    make.add_argument("--panel", type=Path, required=True)
    make.add_argument("--fasta", type=Path, required=True)
    make.add_argument("--fastq", type=Path, required=True)
    make.set_defaults(function=command_make_queries)

    split = subparsers.add_parser("split-fasta")
    split.add_argument("--input", type=Path, required=True)
    split.add_argument("--output-dir", type=Path, required=True)
    split.add_argument("--records-per-chunk", type=int, default=64)
    split.set_defaults(function=command_split_fasta)

    for name, function in (("normalize-sam", command_normalize_sam),
                           ("normalize-blast", command_normalize_blast)):
        command = subparsers.add_parser(name)
        command.add_argument("--panel", type=Path, required=True)
        command.add_argument("--reference", type=Path, required=True)
        command.add_argument("--fai", type=Path, required=True)
        command.add_argument("--input", type=Path, required=True)
        command.add_argument("--output", type=Path, required=True)
        command.add_argument("--anchor", type=int, default=12)
        command.add_argument("--max-mismatches", type=int, default=3)
        if name == "normalize-sam":
            command.add_argument("--include-xa", action="store_true")
        command.set_defaults(function=function)

    assemble = subparsers.add_parser("assemble")
    assemble.add_argument("--bindings", type=Path, required=True)
    assemble.add_argument("--output", type=Path, required=True)
    assemble.add_argument("--min-amplicon", type=int, default=50)
    assemble.add_argument("--max-amplicon", type=int, default=3000)
    assemble.set_defaults(function=command_assemble)

    compare = subparsers.add_parser("compare")
    compare.add_argument("--panel", type=Path, required=True)
    compare.add_argument("--truth", type=Path, required=True)
    compare.add_argument("--prediction", type=Path, required=True)
    compare.add_argument("--kind", choices=("bindings", "products"), required=True)
    compare.add_argument("--summary", type=Path, required=True)
    compare.add_argument("--json", type=Path, required=True)
    compare.add_argument("--per-pair", type=Path, required=True)
    compare.add_argument("--bootstrap-replicates", type=int, default=10000)
    compare.set_defaults(function=command_compare)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
