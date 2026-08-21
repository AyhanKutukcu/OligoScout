#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import statistics
import subprocess
import sys


if len(sys.argv) != 7:
    print(
        "Usage: benchmark_multiplex_strategy_router.py "
        "<benchmark> <manifest> <index_dir> "
        "<chromosome> <panel.tsv> <output_dir>"
    )
    sys.exit(2)


benchmark = Path(sys.argv[1])
manifest = Path(sys.argv[2])
index_dir = Path(sys.argv[3])
chromosome = sys.argv[4]
panel = Path(sys.argv[5])
output_dir = Path(sys.argv[6])

raw_dir = output_dir / "raw"

raw_dir.mkdir(
    parents=True,
    exist_ok=True
)


panel_sizes = [
    8,
    12,
    16,
    20,
    24,
    32,
]


outer_repetitions = 3

inner_repetitions = 10


def parse_number(
    text,
    key,
    integer=False,
):
    match = re.search(
        rf'(?m)^{re.escape(key)}[ \t]+'
        rf'([0-9]+(?:\.[0-9]+)?)',
        text
    )

    if match is None:
        raise RuntimeError(
            f"Cannot parse {key}"
        )

    if integer:
        return int(
            float(
                match.group(1)
            )
        )

    return float(
        match.group(1)
    )


records = []


for pair_count in panel_sizes:

    run_records = []


    for outer in range(
        1,
        outer_repetitions + 1
    ):

        command = [
            "taskset",
            "-c",
            "2",

            str(benchmark),

            str(manifest),

            str(index_dir),

            chromosome,

            str(panel),

            str(pair_count),

            str(inner_repetitions),
        ]


        print()
        print(
            "RUN",
            f"pairs={pair_count}",
            f"outer={outer}",
            f"inner={inner_repetitions}",
            sep="\t"
        )


        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )


        raw_path = (
            raw_dir
            /
            (
                f"chr22_pairs{pair_count}"
                f"_outer{outer}.txt"
            )
        )


        raw_path.write_text(
            completed.stdout
        )


        print(
            completed.stdout
        )


        if completed.returncode != 0:
            raise RuntimeError(
                "Production benchmark failed: "
                f"pairs={pair_count}, "
                f"outer={outer}, "
                f"exit={completed.returncode}"
            )


        actual_pairs = parse_number(
            completed.stdout,
            "pair_count",
            integer=True
        )


        if actual_pairs != pair_count:
            raise RuntimeError(
                "Requested/actual pair count mismatch: "
                f"requested={pair_count}, "
                f"actual={actual_pairs}"
            )


        if (
            re.search(
                r'(?m)^VERIFY_SHARD_EQUIVALENT'
                r'[ \t]+YES$',
                completed.stdout
            )
            is None
        ):
            raise RuntimeError(
                "Production V1/V2 equivalence failed "
                f"at pair_count={pair_count}"
            )


        baseline_ms = parse_number(
            completed.stdout,
            "baseline_join_ms"
        )


        global_ms = parse_number(
            completed.stdout,
            "global_join_ms"
        )


        reported_speedup = parse_number(
            completed.stdout,
            "join_speedup"
        )


        if global_ms <= 0.0:
            raise RuntimeError(
                "Invalid global join runtime"
            )


        computed_speedup = (
            baseline_ms
            /
            global_ms
        )


        relative_error = abs(
            computed_speedup
            -
            reported_speedup
        ) / max(
            reported_speedup,
            1e-12
        )


        if relative_error > 0.02:
            raise RuntimeError(
                "Reported/computed speedup mismatch: "
                f"reported={reported_speedup}, "
                f"computed={computed_speedup}"
            )


        run_records.append(
            {
                "outer": outer,
                "baseline_ms": baseline_ms,
                "global_ms": global_ms,
                "speedup": computed_speedup,
            }
        )


    median_baseline = statistics.median(
        row["baseline_ms"]
        for row in run_records
    )


    median_global = statistics.median(
        row["global_ms"]
        for row in run_records
    )


    median_speedup = statistics.median(
        row["speedup"]
        for row in run_records
    )


    min_speedup = min(
        row["speedup"]
        for row in run_records
    )


    max_speedup = max(
        row["speedup"]
        for row in run_records
    )


    records.append(
        {
            "pairs": pair_count,
            "median_v1_ms": median_baseline,
            "median_v2_ms": median_global,
            "median_speedup": median_speedup,
            "min_speedup": min_speedup,
            "max_speedup": max_speedup,
            "runs": run_records,
        }
    )


# ------------------------------------------------------------
# Stable empirical crossover.
#
# Primary rule:
#   median V1/V2 >= 1.05
#   and every larger measured size median >= 1.00
#
# Fallback:
#   first size from which all measured sizes > 1.00
#
# This avoids choosing a threshold from one noisy sample.
# ------------------------------------------------------------

threshold = None

selection_rule = None


for index, row in enumerate(records):

    if row["median_speedup"] < 1.05:
        continue


    if all(
        later["median_speedup"] >= 1.00
        for later in records[index:]
    ):
        threshold = row["pairs"]

        selection_rule = (
            "median>=1.05_and_all_larger>=1.00"
        )

        break


if threshold is None:

    for index, row in enumerate(records):

        if all(
            later["median_speedup"] > 1.00
            for later in records[index:]
        ):
            threshold = row["pairs"]

            selection_rule = (
                "fallback_all_from_threshold_gt_1.00"
            )

            break


if threshold is None:
    raise RuntimeError(
        "No stable V1->V2 crossover found "
        "within measured sizes."
    )


summary_path = (
    output_dir
    /
    "production_chr22_crossover.tsv"
)


with summary_path.open(
    "w",
    newline=""
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t"
    )


    writer.writerow(
        [
            "pairs",
            "median_v1_ms",
            "median_v2_ms",
            "median_speedup_v1_over_v2",
            "min_speedup",
            "max_speedup",
            "winner",
        ]
    )


    for row in records:

        winner = (
            "V2"
            if row["median_speedup"] > 1.0
            else
            "V1"
        )


        writer.writerow(
            [
                row["pairs"],
                f'{row["median_v1_ms"]:.9f}',
                f'{row["median_v2_ms"]:.9f}',
                f'{row["median_speedup"]:.9f}',
                f'{row["min_speedup"]:.9f}',
                f'{row["max_speedup"]:.9f}',
                winner,
            ]
        )


detail_path = (
    output_dir
    /
    "production_chr22_crossover_runs.tsv"
)


with detail_path.open(
    "w",
    newline=""
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t"
    )


    writer.writerow(
        [
            "pairs",
            "outer_run",
            "inner_repetitions",
            "v1_ms",
            "v2_ms",
            "speedup_v1_over_v2",
        ]
    )


    for row in records:

        for run in row["runs"]:

            writer.writerow(
                [
                    row["pairs"],
                    run["outer"],
                    inner_repetitions,
                    f'{run["baseline_ms"]:.9f}',
                    f'{run["global_ms"]:.9f}',
                    f'{run["speedup"]:.9f}',
                ]
            )


(output_dir / "measured_threshold.txt").write_text(
    f"{threshold}\n"
)


(output_dir / "selection_rule.txt").write_text(
    selection_rule
    +
    "\n"
)


print()
print(
    "================================================"
)

print(
    " PRODUCTION chr22 V1/V2 CROSSOVER"
)

print(
    "================================================"
)


for row in records:

    winner = (
        "V2"
        if row["median_speedup"] > 1.0
        else
        "V1"
    )


    print(
        "PAIRS",
        row["pairs"],
        f'V1_MS={row["median_v1_ms"]:.6f}',
        f'V2_MS={row["median_v2_ms"]:.6f}',
        f'SPEEDUP={row["median_speedup"]:.4f}',
        f'RANGE={row["min_speedup"]:.4f}'
        f'-{row["max_speedup"]:.4f}',
        f'WINNER={winner}',
        sep="\t"
    )


print(
    "MEASURED_THRESHOLD",
    threshold,
    sep="\t"
)


print(
    "SELECTION_RULE",
    selection_rule,
    sep="\t"
)


print(
    "PRODUCTION_CROSSOVER_COMPLETE",
    "YES",
    sep="\t"
)
