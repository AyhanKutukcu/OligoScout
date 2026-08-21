#!/usr/bin/env python3

import math
import statistics
import sys
from pathlib import Path


if len(sys.argv) != 2:
    raise SystemExit(
        "Usage: analyze_anchor_oracle.py <oracle.tsv>"
    )

path = Path(sys.argv[1])

if not path.exists():
    raise SystemExit(
        f"ERROR: file not found: {path}"
    )


rows = []

for line in path.read_text().splitlines():
    if not line.startswith("PAIR"):
        continue

    f = line.split()

    if len(f) < 12:
        raise RuntimeError(
            f"Malformed PAIR row: {line}"
        )

    rows.append(
        {
            "idx": int(f[1]),
            "adaptive_us": float(f[2]),
            "p1_us": float(f[3]),
            "p2_us": float(f[4]),
            "oracle_us": float(f[5]),
            "p1_max": int(f[6]),
            "p2_max": int(f[7]),
            "p1_total": int(f[8]),
            "p2_total": int(f[9]),
            "adaptive_anchor": f[10],
            "oracle_anchor": f[11],
        }
    )


if not rows:
    raise SystemExit(
        "ERROR: no PAIR rows found."
    )


def mean(values):
    return statistics.fmean(values)


def median(values):
    return statistics.median(values)


def pct_gain(old, new):
    if old == 0:
        return 0.0

    return 100.0 * (1.0 - new / old)


def forced_time_for_anchor(row, anchor):
    return (
        row["p1_us"]
        if anchor == "P1"
        else row["p2_us"]
    )


def oracle_anchor(row):
    return (
        "P1"
        if row["p1_us"] <= row["p2_us"]
        else "P2"
    )


def relative_forced_gap(row):
    low = min(
        row["p1_us"],
        row["p2_us"]
    )

    high = max(
        row["p1_us"],
        row["p2_us"]
    )

    if low == 0:
        return 0.0

    return (
        high - low
    ) / low


# ============================================================
# 1. Timing noise:
# adaptive vs forced SAME anchor
# ============================================================

same_anchor_abs_pct = []
same_anchor_delta_us = []

for r in rows:
    same_forced = forced_time_for_anchor(
        r,
        r["adaptive_anchor"]
    )

    same_anchor_delta_us.append(
        r["adaptive_us"] -
        same_forced
    )

    denominator = max(
        same_forced,
        1e-12
    )

    same_anchor_abs_pct.append(
        100.0 *
        abs(
            r["adaptive_us"] -
            same_forced
        )
        /
        denominator
    )


# ============================================================
# 2. Noise-resistant policy comparison
#
# Do NOT compare adaptive_us directly with oracle_us.
# Compare:
#
# forced timing corresponding to adaptive's chosen anchor
# vs
# min(forced P1, forced P2)
# ============================================================

current_policy_forced = []

oracle_times = []

current_correct = 0

for r in rows:
    current_policy_forced.append(
        forced_time_for_anchor(
            r,
            r["adaptive_anchor"]
        )
    )

    oracle_times.append(
        min(
            r["p1_us"],
            r["p2_us"]
        )
    )

    if (
        r["adaptive_anchor"] ==
        oracle_anchor(r)
    ):
        current_correct += 1


current_forced_mean = mean(
    current_policy_forced
)

oracle_mean = mean(
    oracle_times
)


# ============================================================
# 3. Margin-aware accuracy
# ============================================================

margin_stats = []

for threshold in (
    0.00,
    0.02,
    0.05,
    0.10,
    0.20,
):
    subset = [
        r
        for r in rows
        if relative_forced_gap(r) >= threshold
    ]

    if not subset:
        margin_stats.append(
            (
                threshold,
                0,
                0.0
            )
        )
        continue

    correct = sum(
        1
        for r in subset
        if (
            r["adaptive_anchor"] ==
            oracle_anchor(r)
        )
    )

    margin_stats.append(
        (
            threshold,
            len(subset),
            100.0 *
            correct /
            len(subset)
        )
    )


# ============================================================
# 4. Simple deterministic anchor rules
# ============================================================

def low_max(r):
    if r["p1_max"] < r["p2_max"]:
        return "P1"

    if r["p2_max"] < r["p1_max"]:
        return "P2"

    return (
        "P1"
        if r["p1_total"] <= r["p2_total"]
        else "P2"
    )


def high_max(r):
    if r["p1_max"] > r["p2_max"]:
        return "P1"

    if r["p2_max"] > r["p1_max"]:
        return "P2"

    return (
        "P1"
        if r["p1_total"] >= r["p2_total"]
        else "P2"
    )


def low_total(r):
    return (
        "P1"
        if r["p1_total"] <= r["p2_total"]
        else "P2"
    )


def high_total(r):
    return (
        "P1"
        if r["p1_total"] >= r["p2_total"]
        else "P2"
    )


def low_logsum(r):
    p1 = (
        math.log1p(r["p1_max"]) +
        math.log1p(r["p1_total"])
    )

    p2 = (
        math.log1p(r["p2_max"]) +
        math.log1p(r["p2_total"])
    )

    return (
        "P1"
        if p1 <= p2
        else "P2"
    )


def high_logsum(r):
    p1 = (
        math.log1p(r["p1_max"]) +
        math.log1p(r["p1_total"])
    )

    p2 = (
        math.log1p(r["p2_max"]) +
        math.log1p(r["p2_total"])
    )

    return (
        "P1"
        if p1 >= p2
        else "P2"
    )


rules = {
    "CURRENT_ADAPTIVE": (
        lambda r: r["adaptive_anchor"]
    ),
    "LOW_MAX": low_max,
    "HIGH_MAX": high_max,
    "LOW_TOTAL": low_total,
    "HIGH_TOTAL": high_total,
    "LOW_LOGSUM": low_logsum,
    "HIGH_LOGSUM": high_logsum,
}


rule_results = []

for name, rule in rules.items():

    selected = []

    correct = 0

    for r in rows:
        anchor = rule(r)

        selected.append(
            forced_time_for_anchor(
                r,
                anchor
            )
        )

        if anchor == oracle_anchor(r):
            correct += 1

    selected_mean = mean(
        selected
    )

    rule_results.append(
        (
            name,
            selected_mean,
            pct_gain(
                current_forced_mean,
                selected_mean
            ),
            100.0 *
            correct /
            len(rows)
        )
    )


# ============================================================
# 5. Small cross-validated score search
#
# score =
# alpha * log(max_seed + 1)
# +
# (1-alpha) * log(total_seed + 1)
#
# Train alpha/polarity on 4 folds,
# evaluate on held-out fifth.
# ============================================================

alphas = [
    -1.00 + 0.05 * i
    for i in range(61)
]


def score(r, primer, alpha):
    if primer == 1:
        m = r["p1_max"]
        t = r["p1_total"]
    else:
        m = r["p2_max"]
        t = r["p2_total"]

    return (
        alpha * math.log1p(m)
        +
        (1.0 - alpha) *
        math.log1p(t)
    )


def select_score_anchor(
    r,
    alpha,
    polarity
):
    s1 = score(
        r,
        1,
        alpha
    )

    s2 = score(
        r,
        2,
        alpha
    )

    if polarity == "LOW":
        return (
            "P1"
            if s1 <= s2
            else "P2"
        )

    return (
        "P1"
        if s1 >= s2
        else "P2"
    )


cv_selected_times = []
cv_oracle_times = []
cv_current_times = []
cv_correct = 0
cv_parameters = []


for fold in range(5):

    train = [
        r
        for r in rows
        if r["idx"] % 5 != fold
    ]

    test = [
        r
        for r in rows
        if r["idx"] % 5 == fold
    ]

    best = None

    for polarity in (
        "LOW",
        "HIGH",
    ):
        for alpha in alphas:

            times = []

            for r in train:
                a = select_score_anchor(
                    r,
                    alpha,
                    polarity
                )

                times.append(
                    forced_time_for_anchor(
                        r,
                        a
                    )
                )

            m = mean(times)

            candidate = (
                m,
                alpha,
                polarity
            )

            if (
                best is None
                or candidate[0] <
                    best[0]
            ):
                best = candidate

    _, best_alpha, best_polarity = best

    cv_parameters.append(
        (
            fold,
            best_alpha,
            best_polarity
        )
    )

    for r in test:

        selected_anchor = (
            select_score_anchor(
                r,
                best_alpha,
                best_polarity
            )
        )

        selected_time = (
            forced_time_for_anchor(
                r,
                selected_anchor
            )
        )

        current_time = (
            forced_time_for_anchor(
                r,
                r["adaptive_anchor"]
            )
        )

        oracle_time = min(
            r["p1_us"],
            r["p2_us"]
        )

        cv_selected_times.append(
            selected_time
        )

        cv_current_times.append(
            current_time
        )

        cv_oracle_times.append(
            oracle_time
        )

        if (
            selected_anchor ==
            oracle_anchor(r)
        ):
            cv_correct += 1


# ============================================================
# OUTPUT
# ============================================================

print(f"rows\t{len(rows)}")

print(
    "same_anchor_noise_median_abs_pct\t"
    f"{median(same_anchor_abs_pct):.3f}"
)

print(
    "same_anchor_noise_mean_abs_pct\t"
    f"{mean(same_anchor_abs_pct):.3f}"
)

print(
    "same_anchor_delta_mean_us\t"
    f"{mean(same_anchor_delta_us):.3f}"
)

print(
    "current_policy_forced_mean_us\t"
    f"{current_forced_mean:.3f}"
)

print(
    "oracle_forced_mean_us\t"
    f"{oracle_mean:.3f}"
)

print(
    "oracle_policy_gain_pct\t"
    f"{pct_gain(current_forced_mean, oracle_mean):.3f}"
)

print(
    "current_policy_oracle_accuracy_pct\t"
    f"{100.0 * current_correct / len(rows):.3f}"
)

print()

print(
    "# MARGIN_AWARE_ACCURACY"
)

print(
    "minimum_gap_pct\tpairs\taccuracy_pct"
)

for threshold, n, accuracy in margin_stats:
    print(
        f"{100.0 * threshold:.0f}\t"
        f"{n}\t"
        f"{accuracy:.3f}"
    )


print()

print(
    "# SIMPLE_RULES"
)

print(
    "rule\tmean_us\tgain_vs_current_pct\toracle_accuracy_pct"
)

for (
    name,
    selected_mean,
    gain,
    accuracy
) in sorted(
    rule_results,
    key=lambda x: x[1]
):
    print(
        f"{name}\t"
        f"{selected_mean:.3f}\t"
        f"{gain:.3f}\t"
        f"{accuracy:.3f}"
    )


print()

print(
    "# CROSS_VALIDATED_LOG_SCORE"
)

for fold, alpha, polarity in cv_parameters:
    print(
        f"fold_{fold}_model\t"
        f"alpha={alpha:.2f}\t"
        f"polarity={polarity}"
    )

cv_mean = mean(
    cv_selected_times
)

cv_current_mean = mean(
    cv_current_times
)

cv_oracle_mean = mean(
    cv_oracle_times
)

print(
    "cv_current_forced_mean_us\t"
    f"{cv_current_mean:.3f}"
)

print(
    "cv_model_mean_us\t"
    f"{cv_mean:.3f}"
)

print(
    "cv_oracle_mean_us\t"
    f"{cv_oracle_mean:.3f}"
)

print(
    "cv_gain_vs_current_pct\t"
    f"{pct_gain(cv_current_mean, cv_mean):.3f}"
)

print(
    "cv_remaining_gap_to_oracle_pct\t"
    f"{pct_gain(cv_mean, cv_oracle_mean):.3f}"
)

print(
    "cv_oracle_accuracy_pct\t"
    f"{100.0 * cv_correct / len(rows):.3f}"
)
