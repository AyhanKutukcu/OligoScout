#!/usr/bin/env python3

from __future__ import annotations

import statistics
import sys
from collections import defaultdict
from pathlib import Path


def mean(values):
    return (
        sum(values) / len(values)
        if values
        else 0.0
    )


def routed_time(
    row,
    threshold,
):
    """
    Real adaptive-router cost:

        estimator
          +
        selected backend
    """
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
    best_mean = None

    for threshold in thresholds:

        value = mean(
            [
                routed_time(
                    row,
                    threshold,
                )
                for row in rows
            ]
        )

        if (
            best_mean is None
            or value < best_mean
        ):
            best_mean = value
            best_threshold = threshold

    return (
        best_threshold,
        best_mean,
    )


def load_rows(path):
    grouped = defaultdict(list)

    with Path(path).open() as handle:

        for line in handle:

            if not line.startswith(
                "PRIMER\t"
            ):
                continue

            f = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if len(f) != 12:
                raise RuntimeError(
                    "Expected 12 PRIMER columns, "
                    f"found {len(f)}."
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
                    "Equality violation in input."
                )

            grouped[
                int(length)
            ].append(
                {
                    "length":
                        int(length),

                    "primer_index":
                        int(primer_index),

                    "max_seed":
                        int(max_seed),

                    "total_seed":
                        int(total_seed),

                    "estimator_us":
                        float(estimator_us),

                    "exhaustive_us":
                        float(exhaustive_us),

                    "candidate_us":
                        float(candidate_us),
                }
            )

    return grouped


def cv_threshold(
    rows,
):
    thresholds = []
    times = []

    candidate_routes = 0

    for fold in range(5):

        train = [
            r
            for r in rows
            if (
                r["primer_index"]
                %
                5
            )
            != fold
        ]

        test = [
            r
            for r in rows
            if (
                r["primer_index"]
                %
                5
            )
            == fold
        ]

        threshold, _ = (
            choose_threshold(
                train
            )
        )

        thresholds.append(
            threshold
        )

        for row in test:

            times.append(
                routed_time(
                    row,
                    threshold,
                )
            )

            candidate_routes += int(
                row["max_seed"]
                <=
                threshold
            )

    return {
        "thresholds":
            thresholds,

        "median":
            statistics.median(
                thresholds
            ),

        "mean_us":
            mean(times),

        "candidate_route_pct":
            (
                100.0
                *
                candidate_routes
                /
                len(rows)
            ),
    }


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: "
            "analyze_length_router_estimator.py "
            "<benchmark.tsv>"
        )

    grouped = load_rows(
        sys.argv[1]
    )

    historical_threshold = 12275

    print(
        "\t".join(
            [
                "length",
                "estimator_mean_us",
                "always_exhaustive_us",
                "always_candidate_us",
                "cv_thresholds",
                "cv_threshold_median",
                "cv_router_us",
                "cv_gain_vs_exhaustive_pct",
                "cv_gain_vs_candidate_pct",
                "cv_candidate_route_pct",
                "t12275_router_us",
                "t12275_gain_vs_exhaustive_pct",
                "t12275_vs_cv_pct",
            ]
        )
    )


    all_rows = []


    for length in sorted(grouped):

        rows = grouped[length]

        all_rows.extend(
            rows
        )

        estimator_mean = mean(
            [
                r["estimator_us"]
                for r in rows
            ]
        )

        exhaustive_mean = mean(
            [
                r["exhaustive_us"]
                for r in rows
            ]
        )

        # Fixed backend policies do NOT need
        # estimator execution.

        candidate_mean = mean(
            [
                r["candidate_us"]
                for r in rows
            ]
        )

        cv = cv_threshold(
            rows
        )

        cv_gain_exhaustive = (
            100.0
            *
            (
                exhaustive_mean
                -
                cv["mean_us"]
            )
            /
            exhaustive_mean
        )

        cv_gain_candidate = (
            100.0
            *
            (
                candidate_mean
                -
                cv["mean_us"]
            )
            /
            candidate_mean
        )

        t12275_mean = mean(
            [
                routed_time(
                    r,
                    historical_threshold,
                )
                for r in rows
            ]
        )

        t12275_gain = (
            100.0
            *
            (
                exhaustive_mean
                -
                t12275_mean
            )
            /
            exhaustive_mean
        )

        t12275_vs_cv = (
            100.0
            *
            (
                t12275_mean
                -
                cv["mean_us"]
            )
            /
            cv["mean_us"]
        )

        print(
            "\t".join(
                [
                    str(length),

                    f"{estimator_mean:.3f}",

                    f"{exhaustive_mean:.3f}",

                    f"{candidate_mean:.3f}",

                    ",".join(
                        str(x)
                        for x in
                        cv["thresholds"]
                    ),

                    f"{cv['median']:.1f}",

                    f"{cv['mean_us']:.3f}",

                    f"{cv_gain_exhaustive:.3f}",

                    f"{cv_gain_candidate:.3f}",

                    f"{cv['candidate_route_pct']:.3f}",

                    f"{t12275_mean:.3f}",

                    f"{t12275_gain:.3f}",

                    f"{t12275_vs_cv:.3f}",
                ]
            )
        )


    # ========================================================
    # Global single-threshold 5-fold CV
    # ========================================================

    global_cv = cv_threshold(
        all_rows
    )

    global_exhaustive = mean(
        [
            r["exhaustive_us"]
            for r in all_rows
        ]
    )

    global_candidate = mean(
        [
            r["candidate_us"]
            for r in all_rows
        ]
    )

    global_estimator = mean(
        [
            r["estimator_us"]
            for r in all_rows
        ]
    )

    global_12275 = mean(
        [
            routed_time(
                r,
                historical_threshold,
            )
            for r in all_rows
        ]
    )

    print()
    print(
        "GLOBAL"
        "\t"
        f"rows={len(all_rows)}"
        "\t"
        f"estimator_mean_us="
        f"{global_estimator:.3f}"
        "\t"
        f"always_exhaustive_us="
        f"{global_exhaustive:.3f}"
        "\t"
        f"always_candidate_us="
        f"{global_candidate:.3f}"
    )

    print(
        "GLOBAL_CV"
        "\t"
        "thresholds="
        +
        ",".join(
            str(x)
            for x in
            global_cv["thresholds"]
        )
        +
        "\t"
        f"median_threshold="
        f"{global_cv['median']:.1f}"
        "\t"
        f"router_us="
        f"{global_cv['mean_us']:.3f}"
        "\t"
        f"candidate_route_pct="
        f"{global_cv['candidate_route_pct']:.3f}"
    )

    print(
        "GLOBAL_T12275"
        "\t"
        f"router_us="
        f"{global_12275:.3f}"
        "\t"
        f"gain_vs_exhaustive_pct="
        f"{(
            100.0 *
            (
                global_exhaustive
                -
                global_12275
            )
            /
            global_exhaustive
        ):.3f}"
        "\t"
        f"delta_vs_global_cv_pct="
        f"{(
            100.0 *
            (
                global_12275
                -
                global_cv['mean_us']
            )
            /
            global_cv['mean_us']
        ):.3f}"
    )


if __name__ == "__main__":
    main()
