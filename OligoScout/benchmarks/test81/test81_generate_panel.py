#!/usr/bin/env python3
"""Deterministically design and difficulty-stratify the Test #81 panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import mmap
import random
import subprocess
from dataclasses import dataclass, asdict
from pathlib import Path


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


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
        length, offset, line_bases, line_bytes = self.index[name]
        if not (0 <= start <= end <= length):
            raise ValueError(f"FASTA interval out of range: {name}:{start}-{end}")
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


@dataclass
class Candidate:
    source_id: int
    chromosome: str
    window_start: int
    primer1: str
    primer2: str
    expected_start: int
    expected_end: int
    expected_length: int
    primer1_anchor_count: int = 0
    primer2_anchor_count: int = 0
    difficulty_score: int = 0
    difficulty_stratum: str = ""


def make_windows(
    reference: IndexedFasta,
    count: int,
    window_length: int,
    seed: int,
) -> list[tuple[int, str, int, str]]:
    rng = random.Random(seed)
    weighted: list[str] = []
    cumulative: list[int] = []
    total = 0
    for name, (length, *_rest) in reference.index.items():
        capacity = max(0, length - window_length + 1)
        if capacity:
            weighted.append(name)
            total += capacity
            cumulative.append(total)
    windows: list[tuple[int, str, int, str]] = []
    seen: set[tuple[str, int]] = set()
    attempts = 0
    while len(windows) < count and attempts < count * 50:
        attempts += 1
        ticket = rng.randrange(total)
        chromosome_index = next(i for i, edge in enumerate(cumulative) if ticket < edge)
        chromosome = weighted[chromosome_index]
        chromosome_length = reference.index[chromosome][0]
        start = rng.randrange(chromosome_length - window_length + 1)
        key = (chromosome, start)
        if key in seen:
            continue
        sequence = reference.fetch(chromosome, start, start + window_length)
        if set(sequence) <= {"A", "C", "G", "T"}:
            seen.add(key)
            windows.append((len(windows), chromosome, start, sequence))
    if len(windows) != count:
        raise RuntimeError(f"Only {len(windows)} strict-DNA windows found out of {count}")
    return windows


def write_primer3_input(path: Path, windows: list[tuple[int, str, int, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as output:
        for source_id, chromosome, start, sequence in windows:
            output.write(f"SEQUENCE_ID={source_id}|{chromosome}|{start}\n")
            output.write(f"SEQUENCE_TEMPLATE={sequence}\n")
            output.write("SEQUENCE_TARGET=300,200\n")
            output.write("PRIMER_TASK=generic\n")
            output.write("PRIMER_PICK_LEFT_PRIMER=1\n")
            output.write("PRIMER_PICK_INTERNAL_OLIGO=0\n")
            output.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
            output.write("PRIMER_OPT_SIZE=24\nPRIMER_MIN_SIZE=24\nPRIMER_MAX_SIZE=24\n")
            output.write("PRIMER_OPT_TM=60.0\nPRIMER_MIN_TM=57.0\nPRIMER_MAX_TM=64.0\n")
            output.write("PRIMER_MIN_GC=25.0\nPRIMER_MAX_GC=75.0\n")
            output.write("PRIMER_MAX_NS_ACCEPTED=0\n")
            output.write("PRIMER_PRODUCT_SIZE_RANGE=150-500\n")
            output.write("PRIMER_NUM_RETURN=1\n")
            output.write("PRIMER_EXPLAIN_FLAG=1\n")
            output.write("=\n")


def parse_primer3(path: Path, reference: IndexedFasta) -> list[Candidate]:
    records: list[dict[str, str]] = []
    record: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            line = raw.rstrip("\n")
            if line == "=":
                records.append(record)
                record = {}
            elif "=" in line:
                key, value = line.split("=", 1)
                record[key] = value
    candidates: list[Candidate] = []
    seen: set[tuple[str, int, int, str, str]] = set()
    for row in records:
        required = {
            "SEQUENCE_ID", "PRIMER_LEFT_0", "PRIMER_RIGHT_0",
            "PRIMER_LEFT_0_SEQUENCE", "PRIMER_RIGHT_0_SEQUENCE",
            "PRIMER_PAIR_0_PRODUCT_SIZE",
        }
        if not required <= row.keys():
            continue
        source_text, chromosome, window_text = row["SEQUENCE_ID"].split("|")
        window_start = int(window_text)
        left_start, left_length = map(int, row["PRIMER_LEFT_0"].split(","))
        right_position, right_length = map(int, row["PRIMER_RIGHT_0"].split(","))
        primer1 = row["PRIMER_LEFT_0_SEQUENCE"].upper()
        primer2 = row["PRIMER_RIGHT_0_SEQUENCE"].upper()
        expected_start = window_start + left_start
        expected_end = window_start + right_position + 1
        expected_length = expected_end - expected_start
        if left_length != 24 or right_length != 24 or len(primer1) != 24 or len(primer2) != 24:
            continue
        if expected_length != int(row["PRIMER_PAIR_0_PRODUCT_SIZE"]):
            continue
        if reference.fetch(chromosome, expected_start, expected_start + 24) != primer1:
            continue
        if reverse_complement(reference.fetch(chromosome, expected_end - 24, expected_end)) != primer2:
            continue
        key = (chromosome, expected_start, expected_end, primer1, primer2)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(Candidate(
            int(source_text), chromosome, window_start, primer1, primer2,
            expected_start, expected_end, expected_length,
        ))
    return candidates


def query_jellyfish(jellyfish: Path, database: Path, sequences: list[str]) -> list[int]:
    process = subprocess.run(
        [str(jellyfish), "query", "-i", str(database)],
        input="\n".join(sequences) + "\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    counts = [int(line) for line in process.stdout.splitlines() if line.strip()]
    if len(counts) != len(sequences):
        raise RuntimeError(
            f"Jellyfish returned {len(counts)} counts for {len(sequences)} anchors"
        )
    return counts


def spaced_selection(rows: list[Candidate], count: int) -> list[Candidate]:
    if len(rows) < count:
        raise RuntimeError(f"Stratum has {len(rows)} candidates; {count} required")
    if len(rows) == count:
        return rows
    return [rows[(i * len(rows)) // count] for i in range(count)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--primer3", type=Path, required=True)
    parser.add_argument("--jellyfish", type=Path, required=True)
    parser.add_argument("--kmer-db", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--pairs", type=int, default=1024)
    parser.add_argument("--candidate-windows", type=int, default=6000)
    parser.add_argument("--seed", type=int, default=810024)
    args = parser.parse_args()
    if args.pairs % 4:
        raise SystemExit("--pairs must be divisible by four")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    reference = IndexedFasta(args.reference, args.fai)
    try:
        windows = make_windows(reference, args.candidate_windows, 800, args.seed)
        primer3_input = args.work_dir / "primer3_input.txt"
        primer3_output = args.work_dir / "primer3_output.txt"
        primer3_error = args.work_dir / "primer3_stderr.txt"
        write_primer3_input(primer3_input, windows)
        with primer3_input.open("rb") as source, primer3_output.open("wb") as target, \
                primer3_error.open("wb") as error:
            subprocess.run([str(args.primer3)], stdin=source, stdout=target,
                           stderr=error, check=True)
        candidates = parse_primer3(primer3_output, reference)
        if len(candidates) < args.pairs:
            raise RuntimeError(
                f"Primer3 produced only {len(candidates)} valid unique pairs; "
                f"need {args.pairs}. Increase --candidate-windows."
            )
        anchors: list[str] = []
        for candidate in candidates:
            anchors.extend((candidate.primer1[-12:], candidate.primer2[-12:]))
        counts = query_jellyfish(args.jellyfish, args.kmer_db, anchors)
        for index, candidate in enumerate(candidates):
            candidate.primer1_anchor_count = counts[index * 2]
            candidate.primer2_anchor_count = counts[index * 2 + 1]
            candidate.difficulty_score = max(
                candidate.primer1_anchor_count, candidate.primer2_anchor_count
            )
        candidates.sort(key=lambda row: (
            row.difficulty_score, row.primer1_anchor_count,
            row.primer2_anchor_count, row.chromosome, row.expected_start,
        ))
        names = ["Q1_LOW", "Q2_MEDIUM_LOW", "Q3_MEDIUM_HIGH", "Q4_HIGH"]
        quartiles: list[list[Candidate]] = []
        for index, name in enumerate(names):
            begin = (index * len(candidates)) // 4
            end = ((index + 1) * len(candidates)) // 4
            chosen = spaced_selection(candidates[begin:end], args.pairs // 4)
            for candidate in chosen:
                candidate.difficulty_stratum = name
            quartiles.append(chosen)
        selected: list[Candidate] = []
        for within in range(args.pairs // 4):
            for quartile in quartiles:
                selected.append(quartile[within])

        args.output.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = [
            "pair_id", "primer1_id", "primer1", "primer2_id", "primer2",
            "expected_chromosome", "expected_amplicon_start",
            "expected_amplicon_end_exclusive", "expected_amplicon_length",
            "difficulty_stratum", "primer1_anchor12_count",
            "primer2_anchor12_count", "difficulty_score", "source_window_id",
            "source_window_start",
        ]
        with args.output.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for pair_id, row in enumerate(selected):
                writer.writerow({
                    "pair_id": pair_id,
                    "primer1_id": f"test81_pair{pair_id}_primer1",
                    "primer1": row.primer1,
                    "primer2_id": f"test81_pair{pair_id}_primer2",
                    "primer2": row.primer2,
                    "expected_chromosome": row.chromosome,
                    "expected_amplicon_start": row.expected_start,
                    "expected_amplicon_end_exclusive": row.expected_end,
                    "expected_amplicon_length": row.expected_length,
                    "difficulty_stratum": row.difficulty_stratum,
                    "primer1_anchor12_count": row.primer1_anchor_count,
                    "primer2_anchor12_count": row.primer2_anchor_count,
                    "difficulty_score": row.difficulty_score,
                    "source_window_id": row.source_id,
                    "source_window_start": row.window_start,
                })
        panel_sha = hashlib.sha256(args.output.read_bytes()).hexdigest()
        metadata = {
            "schema": "test81_panel_generation_v1",
            "seed": args.seed,
            "candidate_windows_requested": args.candidate_windows,
            "candidate_pairs_valid": len(candidates),
            "selected_pairs": len(selected),
            "primer_length": 24,
            "anchor_length": 12,
            "strata": {name: args.pairs // 4 for name in names},
            "panel_sha256": panel_sha,
            "difficulty_min": min(row.difficulty_score for row in selected),
            "difficulty_max": max(row.difficulty_score for row in selected),
            "selection": "four anchor12-frequency quartiles, deterministic spaced sample, interleaved",
        }
        args.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        print(f"VALID_CANDIDATE_PAIRS\t{len(candidates)}")
        print(f"SELECTED_PAIRS\t{len(selected)}")
        print(f"PANEL_SHA256\t{panel_sha}")
        print("TEST81_PANEL_GENERATION_COMPLETE\tYES")
    finally:
        reference.close()


if __name__ == "__main__":
    main()
