#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path
import sys


def mean(values):
    return (
        sum(values) / len(values)
        if values
        else 0.0
    )


def load_rows(path):
    rows = []

    with Path(path).open() as handle:
        for line in handle:
            if not line.startswith("PRIMER\t"):
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) != 12:
                raise RuntimeError(
                    f"Expected 12 columns, got {len(f)}"
                )

            (
                _tag,
                length,
                primer_index,
                max_seed,
                total_seed,
                estimator_us,
                exhaustive_us,
                candidate_us,
                winner,
                exhaustive_hits,
                candidate_hits,
                equality,
            ) = f

            if equality != "YES":
                raise RuntimeError(
                    "Equality violation in benchmark."
                )

            rows.append(
                {
                    "length": int(length),
                    "primer_index": int(primer_index),
                    "max_seed": int(max_seed),
                    "estimator_us": float(estimator_us),
                    "exhaustive_us": float(exhaustive_us),
                    "candidate_us": float(candidate_us),
                }
            )

    return rows


def threshold_time(row, threshold):
    backend = (
        row["candidate_us"]
        if row["max_seed"] <= threshold
        else row["exhaustive_us"]
    )

    return (
        row["estimator_us"]
        +
        backend
    )


def choose_threshold(rows):
    thresholds = (
        [-1]
        +
        sorted(
            {
                row["max_seed"]
                for row in rows
            }
        )
    )

    best_threshold = None
    best_time = None

    for threshold in thresholds:

        value = mean(
            [
                threshold_time(
                    row,
                    threshold,
                )
                for row in rows
            ]
        )

        if (
            best_time is None
            or value < best_time
        ):
            best_time = value
            best_threshold = threshold

    return best_threshold


def evaluate_cv(rows):
    policies = defaultdict(list)

    universal_short_thresholds = []

    per_length_thresholds = defaultdict(list)

    for fold in range(5):

        train = [
            r
            for r in rows
            if (
                r["primer_index"] % 5
            )
            != fold
        ]

        test = [
            r
            for r in rows
            if (
                r["primer_index"] % 5
            )
            == fold
        ]

        # --------------------------------------------
        # One threshold trained across lengths 18–22.
        # --------------------------------------------

        train_short = [
            r
            for r in train
            if r["length"] <= 22
        ]

        short_threshold = choose_threshold(
            train_short
        )

        universal_short_thresholds.append(
            short_threshold
        )

        # --------------------------------------------
        # Separate thresholds only for 18–22.
        # --------------------------------------------

        length_threshold = {}

        for length in range(18, 23):

            train_length = [
                r
                for r in train
                if r["length"] == length
            ]

            threshold = choose_threshold(
                train_length
            )

            length_threshold[length] = threshold

            per_length_thresholds[
                length
            ].append(
                threshold
            )

        # --------------------------------------------
        # Evaluate held-out fold.
        # --------------------------------------------

        for row in test:

            length = row["length"]

            # Always exhaustive.
            policies[
                "always_exhaustive"
            ].append(
                row["exhaustive_us"]
            )

            # Always candidate.
            policies[
                "always_candidate"
            ].append(
                row["candidate_us"]
            )

            # Current implementation:
            # non-20 k3 -> exhaustive,
            # 20 nt -> estimator + threshold 12275.
            if length == 20:
                current_time = threshold_time(
                    row,
                    12275,
                )
            else:
                current_time = row[
                    "exhaustive_us"
                ]

            policies[
                "current_router"
            ].append(
                current_time
            )

            # Historical universal threshold:
            # estimator everywhere.
            policies[
                "t12275_all_lengths"
            ].append(
                threshold_time(
                    row,
                    12275,
                )
            )

            # NEW HYBRID A:
            # one learned threshold for 18–22,
            # direct Candidate for 23–35.
            if length <= 22:
                hybrid_universal = (
                    threshold_time(
                        row,
                        short_threshold,
                    )
                )
            else:
                hybrid_universal = row[
                    "candidate_us"
                ]

            policies[
                "hybrid_short_universal"
            ].append(
                hybrid_universal
            )

            # NEW HYBRID B:
            # length-specific threshold for 18–22,
            # direct Candidate for 23–35.
            if length <= 22:
                hybrid_specific = (
                    threshold_time(
                        row,
                        length_threshold[
                            length
                        ],
                    )
                )
            else:
                hybrid_specific = row[
                    "candidate_us"
                ]

            policies[
                "hybrid_short_length_specific"
            ].append(
                hybrid_specific
            )

            # Simple fixed policy:
            # 18–22 threshold 12275,
            # 23–35 Candidate directly.
            if length <= 22:
                fixed_hybrid = (
                    threshold_time(
                        row,
                        12275,
                    )
                )
            else:
                fixed_hybrid = row[
                    "candidate_us"
                ]

            policies[
                "hybrid_short_t12275"
            ].append(
                fixed_hybrid
            )

    baseline = mean(
        policies[
            "always_exhaustive"
        ]
    )

    current = mean(
        policies[
            "current_router"
        ]
    )

    print(
        "policy\tmean_us\t"
        "speedup_vs_exhaustive\t"
        "speedup_vs_current\t"
        "gain_vs_current_pct"
    )

    for name in [
        "always_exhaustive",
        "always_candidate",
        "current_router",
        "t12275_all_lengths",
        "hybrid_short_t12275",
        "hybrid_short_universal",
        "hybrid_short_length_specific",
    ]:
        value = mean(
            policies[name]
        )

        print(
            f"{name}\t"
            f"{value:.3f}\t"
            f"{baseline / value:.3f}\t"
            f"{current / value:.3f}\t"
            f"{100.0 * (current - value) / current:.3f}"
        )

    print()

    print(
        "short_universal_cv_thresholds\t"
        +
        ",".join(
            str(x)
            for x in
            universal_short_thresholds
        )
    )

    print(
        "short_universal_threshold_mean\t"
        f"{mean(universal_short_thresholds):.1f}"
    )

    for length in range(18, 23):

        values = per_length_thresholds[
            length
        ]

        print(
            f"length_{length}_cv_thresholds\t"
            +
            ",".join(
                str(x)
                for x in values
            )
        )


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: analyze_hybrid_router.py "
            "<estimator-benchmark.tsv>"
        )

    rows = load_rows(
        sys.argv[1]
    )

    print(
        f"rows\t{len(rows)}"
    )

    evaluate_cv(
        rows
    )


if __name__ == "__main__":
    main()
