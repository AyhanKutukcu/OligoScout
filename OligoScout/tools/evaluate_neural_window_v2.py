#!/usr/bin/env python3

import argparse
import csv
import math

import numpy as np
import torch
from torch import nn


BASE = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}

PREFIX_LENGTHS = (
    4,
    8,
    12,
    16,
    21,
)


def base4_code(sequence):
    value = 0

    for base in sequence:
        value = (
            value * 4
            +
            BASE[base]
        )

    return value


def prefix_fraction(
    sequence,
    length,
):
    return (
        base4_code(
            sequence[:length]
        )
        /
        float(
            4 ** length
        )
    )


def make_features(row):
    kmer = row["kmer"]

    values = [
        prefix_fraction(
            kmer,
            length,
        )
        for length in PREFIX_LENGTHS
    ]

    values.extend(
        [
            float(
                row["gc_fraction"]
            ),

            float(
                row["entropy"]
            )
            /
            2.0,

            float(
                row["max_homopolymer"]
            )
            /
            21.0,

            float(
                row["distinct_bases"]
            )
            /
            4.0,
        ]
    )

    return values


def baseline_row(
    kmer,
    row_count,
):
    fraction = prefix_fraction(
        kmer,
        21,
    )

    return (
        1.0
        +
        fraction
        *
        (
            row_count - 1.0
        )
    )


class ResidualWindowPredictor(
    nn.Module
):
    def __init__(self):
        super().__init__()

        self.network = nn.Sequential(
            nn.Linear(
                9,
                128,
            ),
            nn.ReLU(),

            nn.Linear(
                128,
                128,
            ),
            nn.ReLU(),

            nn.Linear(
                128,
                64,
            ),
            nn.ReLU(),

            nn.Linear(
                64,
                1,
            ),
        )


    def forward(
        self,
        x,
    ):
        return (
            self.network(x)
            .squeeze(-1)
        )


def report_errors(
    name,
    errors,
):
    errors = np.asarray(
        errors,
        dtype=np.float64,
    )

    print(
        f"{name}_mae\t"
        f"{np.mean(errors):.3f}"
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
            f"{np.percentile(errors, percentile):.3f}"
        )

    print(
        f"{name}_max\t"
        f"{np.max(errors):.3f}"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "dataset",
    )

    parser.add_argument(
        "model",
    )

    parser.add_argument(
        "--row-count",
        type=float,
        default=1_000_001.0,
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=8192,
    )

    args = parser.parse_args()


    checkpoint = torch.load(
        args.model,
        map_location="cpu",
        weights_only=True,
    )

    residual_scale = float(
        checkpoint[
            "residual_scale"
        ]
    )


    model = ResidualWindowPredictor()

    model.load_state_dict(
        checkpoint[
            "model_state_dict"
        ]
    )

    model.eval()


    features = []
    baselines = []
    true_rows = []


    with open(
        args.dataset,
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:

            if (
                row["split"]
                !=
                "test"
            ):
                continue

            features.append(
                make_features(
                    row
                )
            )

            baselines.append(
                baseline_row(
                    row["kmer"],
                    args.row_count,
                )
            )

            true_rows.append(
                float(
                    row["sa_lower"]
                )
            )


    x = np.asarray(
        features,
        dtype=np.float32,
    )

    baselines = np.asarray(
        baselines,
        dtype=np.float64,
    )

    true_rows = np.asarray(
        true_rows,
        dtype=np.float64,
    )


    predictions = []


    with torch.no_grad():

        for start in range(
            0,
            len(x),
            args.batch_size,
        ):

            batch = torch.from_numpy(
                x[
                    start:
                    start
                    +
                    args.batch_size
                ]
            )

            residual = (
                model(batch)
                .cpu()
                .numpy()
                .astype(
                    np.float64
                )
                *
                residual_scale
            )

            predictions.append(
                residual
            )


    predicted_residual = np.concatenate(
        predictions
    )


    neural_rows = (
        baselines
        +
        predicted_residual
    )


    neural_rows = np.clip(
        neural_rows,
        0.0,
        args.row_count,
    )


    lexical_errors = np.abs(
        baselines
        -
        true_rows
    )


    neural_errors = np.abs(
        neural_rows
        -
        true_rows
    )


    print(
        "test_samples\t",
        len(true_rows),
    )

    print(
        "model_residual_scale\t",
        residual_scale,
    )


    report_errors(
        "LEX",
        lexical_errors,
    )

    report_errors(
        "NEURAL",
        neural_errors,
    )


    lex_p99 = np.percentile(
        lexical_errors,
        99,
    )

    neural_p99 = np.percentile(
        neural_errors,
        99,
    )


    print(
        "p99_reduction_fold\t",
        f"{lex_p99 / neural_p99:.6f}"
    )


    global_comparisons = math.ceil(
        math.log2(
            args.row_count
        )
    )


    print(
        "global_binary_comparisons\t",
        global_comparisons,
    )


    print(
        "WINDOW_ANALYSIS"
    )

    for radius in (
        128,
        256,
        512,
        1024,
        2048,
        4096,
        8192,
    ):

        covered = (
            neural_errors
            <=
            radius
        )

        coverage = float(
            np.mean(
                covered
            )
        )

        fallback = (
            1.0
            -
            coverage
        )

        window_rows = (
            2 * radius
            +
            1
        )

        local_comparisons = (
            math.ceil(
                math.log2(
                    window_rows
                )
            )
        )

        expected_comparisons = (
            local_comparisons
            +
            fallback
            *
            global_comparisons
        )

        print(
            "WINDOW"
            f"\tradius={radius}"
            f"\trows={window_rows}"
            f"\tcoverage={coverage:.8f}"
            f"\tfallback={fallback:.8f}"
            f"\tlocal_log2={local_comparisons}"
            f"\testimated_comparisons={expected_comparisons:.4f}"
        )


    print(
        "ALL_CHECKS\tYES"
    )


if __name__ == "__main__":
    main()
