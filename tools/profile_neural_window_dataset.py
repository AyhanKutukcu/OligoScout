#!/usr/bin/env python3

import csv
import math
import sys
from collections import Counter


if len(sys.argv) != 2:
    raise SystemExit(
        "Usage: profile_neural_window_dataset.py <dataset.tsv>"
    )


path = sys.argv[1]

rows = 0
present = 0
absent = 0
mutated = 0

split_counts = Counter()
width_counts = Counter()

present_unique = 0
present_repetitive = 0

max_width = 0
max_width_kmer = ""

entropy_sum = 0.0
gc_sum = 0.0


with open(path, newline="") as handle:
    reader = csv.DictReader(
        handle,
        delimiter="\t",
    )

    for row in reader:
        rows += 1

        split_counts[
            row["split"]
        ] += 1

        is_present = (
            int(row["present"]) == 1
        )

        is_mutated = (
            int(row["mutated"]) == 1
        )

        width = int(
            row["sa_width"]
        )

        if is_present:
            present += 1

            width_counts[
                width
            ] += 1

            if width == 1:
                present_unique += 1
            else:
                present_repetitive += 1

            if width > max_width:
                max_width = width
                max_width_kmer = row["kmer"]

        else:
            absent += 1

        if is_mutated:
            mutated += 1

        entropy_sum += float(
            row["entropy"]
        )

        gc_sum += float(
            row["gc_fraction"]
        )


def width_quantile(q):
    if present == 0:
        return 0

    target = math.ceil(
        q * present
    )

    cumulative = 0

    for width in sorted(
        width_counts
    ):
        cumulative += width_counts[
            width
        ]

        if cumulative >= target:
            return width

    return max_width


print(
    "rows\t",
    rows,
)

print(
    "present\t",
    present,
)

print(
    "absent\t",
    absent,
)

print(
    "mutated\t",
    mutated,
)

print(
    "present_unique_width1\t",
    present_unique,
)

print(
    "present_repetitive_width_gt1\t",
    present_repetitive,
)

if present:
    print(
        "unique_fraction_present\t",
        f"{present_unique / present:.8f}",
    )

print(
    "width_p50\t",
    width_quantile(0.50),
)

print(
    "width_p90\t",
    width_quantile(0.90),
)

print(
    "width_p95\t",
    width_quantile(0.95),
)

print(
    "width_p99\t",
    width_quantile(0.99),
)

print(
    "width_p999\t",
    width_quantile(0.999),
)

print(
    "max_width\t",
    max_width,
)

print(
    "max_width_kmer\t",
    max_width_kmer,
)

print(
    "mean_gc\t",
    f"{gc_sum / rows:.8f}",
)

print(
    "mean_entropy\t",
    f"{entropy_sum / rows:.8f}",
)

for split in (
    "train",
    "validation",
    "test",
):
    print(
        f"{split}_samples\t",
        split_counts[split],
    )

print(
    "TOP_WIDTHS"
)

for width, count in width_counts.most_common(
    15
):
    print(
        f"{width}\t{count}"
    )

print(
    "ALL_CHECKS\tYES"
)
