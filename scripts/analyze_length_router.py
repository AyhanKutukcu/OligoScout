#!/usr/bin/env python3

from __future__ import annotations

import statistics
import sys
from collections import defaultdict
from pathlib import Path


def mean(values):
    return sum(values) / len(values) if values else 0.0


def routed_time(row, threshold):
    if row["max_seed"] <= threshold:
        return row["candidate_us"]
    return row["exhaustive_us"]


def choose_threshold(rows):
    """
    Candidate if:
        max_seed_occurrences <= threshold

    threshold = -1 means always exhaustive.

    max(observed max_seed) effectively permits
    always-candidate.
    """
    thresholds = [-1] + sorted(
        {row["max_seed"] for row in rows}
    )

    best_threshold = None
    best_mean = None

    for threshold in thresholds:
        current = mean(
            [
                routed_time(row, threshold)
                for row in rows
            ]
        )

        if (
            best_mean is None
            or current < best_mean
            or (
                current == best_mean
                and threshold < best_threshold
            )
        ):
            best_mean = current
            best_threshold = threshold

    return best_threshold, best_mean


def route_fraction(rows, threshold):
    if not rows:
        return 0.0

    routed_candidate = sum(
        row["max_seed"] <= threshold
        for row in rows
    )

    return routed_candidate / len(rows)


def decision_accuracy(rows, threshold):
    """
    Fraction where threshold router selects
    the actually faster backend.
    """
    correct = 0

    for row in rows:
        candidate_selected = (
            row["max_seed"] <= threshold
        )

        candidate_faster = (
            row["candidate_us"]
            <
            row["exhaustive_us"]
        )

        if candidate_selected == candidate_faster:
            correct += 1

    return (
        correct / len(rows)
        if rows
        else 0.0
    )


def load_rows(path):
    grouped = defaultdict(list)

    with Path(path).open() as handle:
        for line in handle:
            if not line.startswith("PRIMER\t"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 11:
                raise RuntimeError(
                    f"Unexpected PRIMER column count: "
                    f"{len(fields)}"
                )

            (
                _tag,
                length,
                primer_index,
                max_seed,
                total_seed,
                exhaustive_us,
                candidate_us,
                winner,
                exhaustive_hits,
                candidate_hits,
                equality,
            ) = fields

            if equality != "YES":
                raise RuntimeError(
                    "Input contains equality violation."
                )

            row = {
                "length": int(length),
                "primer_index": int(primer_index),
                "max_seed": int(max_seed),
                "total_seed": int(total_seed),
                "exhaustive_us": float(
                    exhaustive_us
                ),
                "candidate_us": float(
                    candidate_us
                ),
                "winner": winner,
                "exhaustive_hits": int(
                    exhaustive_hits
                ),
                "candidate_hits": int(
                    candidate_hits
                ),
            }

            grouped[row["length"]].append(row)

    return grouped


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: analyze_length_router.py "
            "<length-router.tsv>"
        )

    grouped = load_rows(sys.argv[1])

    print(
        "\t".join(
            [
                "length",
                "rows",
                "candidate_win_pct",
                "always_exhaustive_us",
                "always_candidate_us",
                "oracle_us",
                "full_threshold",
                "full_router_us",
                "full_gain_pct",
                "full_candidate_route_pct",
                "full_accuracy_pct",
                "cv_thresholds",
                "cv_threshold_median",
                "cv_router_us",
                "cv_gain_pct",
                "cv_candidate_route_pct",
                "cv_accuracy_pct",
                "cv_regret_vs_oracle_pct",
                "threshold_12275_us",
                "threshold_12275_gain_pct",
            ]
        )
    )

    for length in sorted(grouped):
        rows = grouped[length]

        exhaustive_mean = mean(
            [
                row["exhaustive_us"]
                for row in rows
            ]
        )

        candidate_mean = mean(
            [
                row["candidate_us"]
                for row in rows
            ]
        )

        oracle_mean = mean(
            [
                min(
                    row["exhaustive_us"],
                    row["candidate_us"],
                )
                for row in rows
            ]
        )

        candidate_win_pct = (
            100.0
            *
            sum(
                row["candidate_us"]
                <
                row["exhaustive_us"]
                for row in rows
            )
            /
            len(rows)
        )

        full_threshold, full_router_mean = (
            choose_threshold(rows)
        )

        full_gain_pct = (
            100.0
            *
            (
                exhaustive_mean
                -
                full_router_mean
            )
            /
            exhaustive_mean
        )

        full_route_pct = (
            100.0
            *
            route_fraction(
                rows,
                full_threshold,
            )
        )

        full_accuracy_pct = (
            100.0
            *
            decision_accuracy(
                rows,
                full_threshold,
            )
        )

        # -----------------------------------
        # 5-fold deterministic CV
        # -----------------------------------

        fold_thresholds = []
        cv_times = []
        cv_selected_candidate = 0
        cv_correct = 0
        cv_rows = 0

        for fold in range(5):
            train = [
                row
                for row in rows
                if (
                    row["primer_index"]
                    % 5
                )
                != fold
            ]

            test = [
                row
                for row in rows
                if (
                    row["primer_index"]
                    % 5
                )
                == fold
            ]

            threshold, _ = (
                choose_threshold(train)
            )

            fold_thresholds.append(
                threshold
            )

            for row in test:
                candidate_selected = (
                    row["max_seed"]
                    <=
                    threshold
                )

                candidate_faster = (
                    row["candidate_us"]
                    <
                    row["exhaustive_us"]
                )

                cv_times.append(
                    row["candidate_us"]
                    if candidate_selected
                    else row["exhaustive_us"]
                )

                cv_selected_candidate += int(
                    candidate_selected
                )

                cv_correct += int(
                    candidate_selected
                    ==
                    candidate_faster
                )

                cv_rows += 1

        cv_router_mean = mean(
            cv_times
        )

        cv_gain_pct = (
            100.0
            *
            (
                exhaustive_mean
                -
                cv_router_mean
            )
            /
            exhaustive_mean
        )

        cv_route_pct = (
            100.0
            *
            cv_selected_candidate
            /
            cv_rows
        )

        cv_accuracy_pct = (
            100.0
            *
            cv_correct
            /
            cv_rows
        )

        cv_regret_pct = (
            100.0
            *
            (
                cv_router_mean
                -
                oracle_mean
            )
            /
            oracle_mean
        )

        cv_threshold_median = (
            statistics.median(
                fold_thresholds
            )
        )

        # -----------------------------------
        # Historical 20-mer threshold check
        # -----------------------------------

        historical_threshold = 12275

        old_threshold_mean = mean(
            [
                routed_time(
                    row,
                    historical_threshold,
                )
                for row in rows
            ]
        )

        old_threshold_gain_pct = (
            100.0
            *
            (
                exhaustive_mean
                -
                old_threshold_mean
            )
            /
            exhaustive_mean
        )

        print(
            "\t".join(
                [
                    str(length),
                    str(len(rows)),
                    f"{candidate_win_pct:.3f}",
                    f"{exhaustive_mean:.3f}",
                    f"{candidate_mean:.3f}",
                    f"{oracle_mean:.3f}",
                    str(full_threshold),
                    f"{full_router_mean:.3f}",
                    f"{full_gain_pct:.3f}",
                    f"{full_route_pct:.3f}",
                    f"{full_accuracy_pct:.3f}",
                    ",".join(
                        str(x)
                        for x in fold_thresholds
                    ),
                    f"{cv_threshold_median:.1f}",
                    f"{cv_router_mean:.3f}",
                    f"{cv_gain_pct:.3f}",
                    f"{cv_route_pct:.3f}",
                    f"{cv_accuracy_pct:.3f}",
                    f"{cv_regret_pct:.3f}",
                    f"{old_threshold_mean:.3f}",
                    f"{old_threshold_gain_pct:.3f}",
                ]
            )
        )


if __name__ == "__main__":
    main()
