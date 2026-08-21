#!/usr/bin/env python3

import csv
import hashlib
import statistics
import sys
from pathlib import Path


def percentile(values, q):
    if not values:
        return 0.0

    values = sorted(values)

    if len(values) == 1:
        return values[0]

    x = q * (len(values) - 1)

    lo = int(x)
    hi = min(
        lo + 1,
        len(values) - 1
    )

    fraction = x - lo

    return (
        values[lo] *
        (1.0 - fraction)
        +
        values[hi] *
        fraction
    )


def load_records(path):
    lines = []

    with open(
        path,
        "r",
        encoding="utf-8"
    ) as handle:

        for line in handle:

            if line.startswith("#"):
                continue

            if not line.strip():
                continue

            lines.append(line)

    reader = csv.DictReader(
        lines,
        delimiter="\t"
    )

    records = []

    for row in reader:

        if row["all_checks"] != "YES":
            raise RuntimeError(
                f"Row {row['primer_index']} "
                f"has all_checks={row['all_checks']}"
            )

        primer = row["primer"]

        digest = hashlib.sha256(
            primer.encode("ascii")
        ).digest()

        fold = (
            int.from_bytes(
                digest[:8],
                "big"
            )
            &
            1
        )

        records.append(
            {
                "index":
                    int(
                        row["primer_index"]
                    ),

                "primer":
                    primer,

                "fold":
                    fold,

                "total":
                    int(
                        row[
                            "total_seed_occurrences"
                        ]
                    ),

                "max":
                    int(
                        row[
                            "max_seed_occurrences"
                        ]
                    ),

                "est":
                    float(
                        row[
                            "estimator_median_us"
                        ]
                    ),

                "ex":
                    float(
                        row[
                            "exhaustive_median_us"
                        ]
                    ),

                "cand":
                    float(
                        row[
                            "candidate_median_us"
                        ]
                    ),

                "hit_count":
                    int(
                        row["hit_count"]
                    ),
            }
        )

    if not records:
        raise RuntimeError(
            "No benchmark records found."
        )

    return records


def backend_winner(row):

    if row["cand"] < row["ex"]:
        return "CANDIDATE"

    return "EXHAUSTIVE"


def router_latency(
    row,
    choose_candidate
):
    backend = (
        row["cand"]
        if choose_candidate
        else row["ex"]
    )

    # Adaptive routing always pays
    # the estimator cost first.
    return (
        row["est"]
        +
        backend
    )


def policy_stats(
    rows,
    selector
):
    latencies = []

    correct = 0

    candidate_choices = 0

    wrong_candidate = 0
    wrong_exhaustive = 0

    for row in rows:

        choose_candidate = (
            selector(row)
        )

        if choose_candidate:
            candidate_choices += 1

        predicted = (
            "CANDIDATE"
            if choose_candidate
            else "EXHAUSTIVE"
        )

        winner = (
            backend_winner(row)
        )

        if predicted == winner:

            correct += 1

        elif choose_candidate:

            wrong_candidate += 1

        else:

            wrong_exhaustive += 1

        latencies.append(
            router_latency(
                row,
                choose_candidate
            )
        )

    return {
        "mean":
            statistics.fmean(
                latencies
            ),

        "median":
            statistics.median(
                latencies
            ),

        "p95":
            percentile(
                latencies,
                0.95
            ),

        "accuracy":
            correct /
            len(rows),

        "candidate_choices":
            candidate_choices,

        "wrong_candidate":
            wrong_candidate,

        "wrong_exhaustive":
            wrong_exhaustive,
    }


def baseline_stats(
    records,
    field
):
    values = [
        row[field]
        for row in records
    ]

    return {
        "mean":
            statistics.fmean(
                values
            ),

        "median":
            statistics.median(
                values
            ),

        "p95":
            percentile(
                values,
                0.95
            ),
    }


def print_baseline(
    name,
    stats
):
    print(
        f"{name}\t"
        f"mean={stats['mean']:.3f}\t"
        f"median={stats['median']:.3f}\t"
        f"p95={stats['p95']:.3f}"
    )


def optimize_1d(
    rows,
    metric
):
    thresholds = (
        [-1]
        +
        sorted(
            {
                row[metric]
                for row in rows
            }
        )
    )

    best = None

    for threshold in thresholds:

        selector = (
            lambda row,
            t=threshold,
            m=metric:
                row[m] <= t
        )

        stats = policy_stats(
            rows,
            selector
        )

        key = (
            stats["mean"],
            stats["p95"],
            stats["wrong_candidate"],
        )

        if (
            best is None
            or
            key < best["key"]
        ):
            best = {
                "threshold":
                    threshold,

                "stats":
                    stats,

                "key":
                    key,
            }

    return best


def optimize_2d(rows):

    total_thresholds = (
        [-1]
        +
        sorted(
            {
                row["total"]
                for row in rows
            }
        )
    )

    max_thresholds = (
        [-1]
        +
        sorted(
            {
                row["max"]
                for row in rows
            }
        )
    )

    best = None

    for total_threshold in total_thresholds:

        for max_threshold in max_thresholds:

            selector = (
                lambda row,
                tt=total_threshold,
                mt=max_threshold:
                    (
                        row["total"] <= tt
                        and
                        row["max"] <= mt
                    )
            )

            stats = policy_stats(
                rows,
                selector
            )

            key = (
                stats["mean"],
                stats["p95"],
                stats[
                    "wrong_candidate"
                ],
            )

            if (
                best is None
                or
                key < best["key"]
            ):
                best = {
                    "total_threshold":
                        total_threshold,

                    "max_threshold":
                        max_threshold,

                    "stats":
                        stats,

                    "key":
                        key,
                }

    return best


def summarize_combined(
    fold_models,
    exhaustive_mean
):
    latencies = []

    correct = 0
    total = 0

    candidate_choices = 0

    wrong_candidate = 0
    wrong_exhaustive = 0

    for (
        holdout,
        selector
    ) in fold_models:

        for row in holdout:

            choose_candidate = (
                selector(row)
            )

            if choose_candidate:
                candidate_choices += 1

            predicted = (
                "CANDIDATE"
                if choose_candidate
                else "EXHAUSTIVE"
            )

            winner = (
                backend_winner(row)
            )

            if predicted == winner:

                correct += 1

            elif choose_candidate:

                wrong_candidate += 1

            else:

                wrong_exhaustive += 1

            latencies.append(
                router_latency(
                    row,
                    choose_candidate
                )
            )

            total += 1

    mean_latency = (
        statistics.fmean(
            latencies
        )
    )

    return {
        "mean":
            mean_latency,

        "median":
            statistics.median(
                latencies
            ),

        "p95":
            percentile(
                latencies,
                0.95
            ),

        "speedup":
            exhaustive_mean /
            mean_latency,

        "accuracy":
            correct /
            total,

        "candidate_choices":
            candidate_choices,

        "wrong_candidate":
            wrong_candidate,

        "wrong_exhaustive":
            wrong_exhaustive,
    }


def main():

    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: "
            "analyze_sensitive_router_thresholds.py "
            "<router_raw.tsv>"
        )

    records = load_records(
        Path(
            sys.argv[1]
        )
    )

    candidate_wins = sum(
        backend_winner(row) ==
        "CANDIDATE"
        for row in records
    )

    exhaustive_wins = (
        len(records)
        -
        candidate_wins
    )

    print(
        f"records={len(records)}"
    )

    print(
        f"candidate_wins="
        f"{candidate_wins}"
    )

    print(
        f"exhaustive_wins="
        f"{exhaustive_wins}"
    )

    print(
        "candidate_win_pct="
        f"{100.0 * candidate_wins / len(records):.3f}"
    )

    print()

    exhaustive_stats = (
        baseline_stats(
            records,
            "ex"
        )
    )

    candidate_stats = (
        baseline_stats(
            records,
            "cand"
        )
    )

    estimator_stats = (
        baseline_stats(
            records,
            "est"
        )
    )

    oracle_backend = [
        min(
            row["ex"],
            row["cand"]
        )
        for row in records
    ]

    oracle_router = [
        (
            row["est"]
            +
            min(
                row["ex"],
                row["cand"]
            )
        )
        for row in records
    ]

    print(
        "=== BASELINES ==="
    )

    print_baseline(
        "always_exhaustive",
        exhaustive_stats
    )

    print_baseline(
        "always_candidate",
        candidate_stats
    )

    print_baseline(
        "estimator_only",
        estimator_stats
    )

    print_baseline(
        "oracle_backend_no_estimator",
        {
            "mean":
                statistics.fmean(
                    oracle_backend
                ),

            "median":
                statistics.median(
                    oracle_backend
                ),

            "p95":
                percentile(
                    oracle_backend,
                    0.95
                ),
        }
    )

    print_baseline(
        "oracle_router_with_estimator",
        {
            "mean":
                statistics.fmean(
                    oracle_router
                ),

            "median":
                statistics.median(
                    oracle_router
                ),

            "p95":
                percentile(
                    oracle_router,
                    0.95
                ),
        }
    )

    print()

    fold_sizes = {
        fold:
            sum(
                row["fold"] ==
                fold
                for row in records
            )
        for fold in (
            0,
            1
        )
    }

    print(
        f"fold0="
        f"{fold_sizes[0]}"
    )

    print(
        f"fold1="
        f"{fold_sizes[1]}"
    )

    print()

    families = {
        "TOTAL":
            [],

        "MAX":
            [],

        "TOTAL_AND_MAX":
            [],
    }

    for train_fold in (
        0,
        1
    ):

        holdout_fold = (
            1 -
            train_fold
        )

        train = [
            row
            for row in records
            if row["fold"] ==
                train_fold
        ]

        holdout = [
            row
            for row in records
            if row["fold"] ==
                holdout_fold
        ]

        if (
            not train
            or
            not holdout
        ):
            raise RuntimeError(
                "Hash split produced an empty "
                "train or holdout fold."
            )

        print(
            f"=== TRAIN FOLD "
            f"{train_fold} "
            f"-> HOLDOUT "
            f"{holdout_fold} ==="
        )

        # --------------------------------------------
        # TOTAL
        # --------------------------------------------

        total_fit = (
            optimize_1d(
                train,
                "total"
            )
        )

        total_threshold = (
            total_fit[
                "threshold"
            ]
        )

        total_selector = (
            lambda row,
            t=total_threshold:
                row["total"] <= t
        )

        total_hold = (
            policy_stats(
                holdout,
                total_selector
            )
        )

        families[
            "TOTAL"
        ].append(
            (
                holdout,
                total_selector
            )
        )

        print(
            "TOTAL\t"
            f"threshold="
            f"{total_threshold}\t"
            f"train_mean="
            f"{total_fit['stats']['mean']:.3f}\t"
            f"hold_mean="
            f"{total_hold['mean']:.3f}\t"
            f"hold_median="
            f"{total_hold['median']:.3f}\t"
            f"hold_p95="
            f"{total_hold['p95']:.3f}\t"
            f"accuracy="
            f"{100.0 * total_hold['accuracy']:.2f}%\t"
            f"candidate_choices="
            f"{total_hold['candidate_choices']}\t"
            f"wrong_candidate="
            f"{total_hold['wrong_candidate']}\t"
            f"wrong_exhaustive="
            f"{total_hold['wrong_exhaustive']}"
        )

        # --------------------------------------------
        # MAX
        # --------------------------------------------

        max_fit = (
            optimize_1d(
                train,
                "max"
            )
        )

        max_threshold = (
            max_fit[
                "threshold"
            ]
        )

        max_selector = (
            lambda row,
            t=max_threshold:
                row["max"] <= t
        )

        max_hold = (
            policy_stats(
                holdout,
                max_selector
            )
        )

        families[
            "MAX"
        ].append(
            (
                holdout,
                max_selector
            )
        )

        print(
            "MAX\t"
            f"threshold="
            f"{max_threshold}\t"
            f"train_mean="
            f"{max_fit['stats']['mean']:.3f}\t"
            f"hold_mean="
            f"{max_hold['mean']:.3f}\t"
            f"hold_median="
            f"{max_hold['median']:.3f}\t"
            f"hold_p95="
            f"{max_hold['p95']:.3f}\t"
            f"accuracy="
            f"{100.0 * max_hold['accuracy']:.2f}%\t"
            f"candidate_choices="
            f"{max_hold['candidate_choices']}\t"
            f"wrong_candidate="
            f"{max_hold['wrong_candidate']}\t"
            f"wrong_exhaustive="
            f"{max_hold['wrong_exhaustive']}"
        )

        # --------------------------------------------
        # TOTAL AND MAX
        # --------------------------------------------

        two_fit = (
            optimize_2d(
                train
            )
        )

        total_t = (
            two_fit[
                "total_threshold"
            ]
        )

        max_t = (
            two_fit[
                "max_threshold"
            ]
        )

        two_selector = (
            lambda row,
            tt=total_t,
            mt=max_t:
                (
                    row["total"] <= tt
                    and
                    row["max"] <= mt
                )
        )

        two_hold = (
            policy_stats(
                holdout,
                two_selector
            )
        )

        families[
            "TOTAL_AND_MAX"
        ].append(
            (
                holdout,
                two_selector
            )
        )

        print(
            "TOTAL_AND_MAX\t"
            f"total_threshold="
            f"{total_t}\t"
            f"max_threshold="
            f"{max_t}\t"
            f"train_mean="
            f"{two_fit['stats']['mean']:.3f}\t"
            f"hold_mean="
            f"{two_hold['mean']:.3f}\t"
            f"hold_median="
            f"{two_hold['median']:.3f}\t"
            f"hold_p95="
            f"{two_hold['p95']:.3f}\t"
            f"accuracy="
            f"{100.0 * two_hold['accuracy']:.2f}%\t"
            f"candidate_choices="
            f"{two_hold['candidate_choices']}\t"
            f"wrong_candidate="
            f"{two_hold['wrong_candidate']}\t"
            f"wrong_exhaustive="
            f"{two_hold['wrong_exhaustive']}"
        )

        print()

    print(
        "=== COMBINED 2-FOLD "
        "CROSS-VALIDATED PERFORMANCE ==="
    )

    for (
        family_name,
        fold_models
    ) in families.items():

        stats = (
            summarize_combined(
                fold_models,
                exhaustive_stats[
                    "mean"
                ]
            )
        )

        print(
            f"{family_name}\t"
            f"mean="
            f"{stats['mean']:.3f}\t"
            f"median="
            f"{stats['median']:.3f}\t"
            f"p95="
            f"{stats['p95']:.3f}\t"
            f"speedup_vs_exhaustive_mean="
            f"{stats['speedup']:.3f}x\t"
            f"accuracy="
            f"{100.0 * stats['accuracy']:.2f}%\t"
            f"candidate_choices="
            f"{stats['candidate_choices']}\t"
            f"wrong_candidate="
            f"{stats['wrong_candidate']}\t"
            f"wrong_exhaustive="
            f"{stats['wrong_exhaustive']}"
        )


if __name__ == "__main__":
    main()
