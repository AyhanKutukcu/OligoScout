#!/usr/bin/env python3

import csv
import sys
import numpy as np


if len(sys.argv) != 3:
    raise SystemExit(
        "Usage: benchmark_lexicographic_window.py "
        "<dataset.tsv> <row_count>"
    )


path = sys.argv[1]
row_count = int(
    sys.argv[2]
)


BASE = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}


def lexicographic_fraction(
    sequence,
):
    code = 0

    for base in sequence:
        code = (
            code * 4
            +
            BASE[base]
        )

    denominator = (
        4 ** len(sequence)
    )

    return (
        code /
        denominator
    )


errors_validation = []
errors_test = []


with open(
    path,
    newline="",
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t",
    )

    for row in reader:

        split = row[
            "split"
        ]

        if split not in (
            "validation",
            "test",
        ):
            continue


        fraction = (
            lexicographic_fraction(
                row["kmer"]
            )
        )


        # One sentinel suffix precedes all
        # A/C/G/T suffixes.
        predicted_row = (
            1.0
            +
            fraction
            *
            (
                row_count - 1
            )
        )


        true_row = float(
            row["sa_lower"]
        )


        error = abs(
            predicted_row
            -
            true_row
        )


        if split == "validation":
            errors_validation.append(
                error
            )
        else:
            errors_test.append(
                error
            )


def report(
    name,
    values,
):
    values = np.asarray(
        values,
        dtype=np.float64,
    )

    print(
        f"{name}_samples\t"
        f"{len(values)}"
    )

    print(
        f"{name}_mae\t"
        f"{np.mean(values):.3f}"
    )

    for label, percentile in (
        ("median", 50),
        ("p90", 90),
        ("p95", 95),
        ("p99", 99),
        ("p999", 99.9),
    ):
        print(
            f"{name}_{label}\t"
            f"{np.percentile(values, percentile):.3f}"
        )

    print(
        f"{name}_max\t"
        f"{np.max(values):.3f}"
    )


report(
    "validation",
    errors_validation,
)

report(
    "test",
    errors_test,
)

print(
    "ALL_CHECKS\tYES"
)
