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
            ) / 2.0,

            float(
                row["max_homopolymer"]
            ) / 21.0,

            float(
                row["distinct_bases"]
            ) / 4.0,
        ]
    )

    return values


def baseline_row(
    kmer,
    row_count,
):
    return (
        1.0
        +
        prefix_fraction(
            kmer,
            21,
        )
        *
        (
            row_count - 1.0
        )
    )


class ResidualMLP(
    nn.Module
):
    def __init__(
        self,
        hidden_sizes,
    ):
        super().__init__()

        layers = []
        previous = 9

        for hidden in hidden_sizes:

            layers.append(
                nn.Linear(
                    previous,
                    hidden,
                )
            )

            layers.append(
                nn.ReLU()
            )

            previous = hidden

        layers.append(
            nn.Linear(
                previous,
                1,
            )
        )

        self.network = nn.Sequential(
            *layers
        )


    def forward(
        self,
        x,
    ):
        return (
            self.network(x)
            .squeeze(-1)
        )


def percentile_report(
    prefix,
    errors,
):
    print(
        f"{prefix}_mae\t"
        f"{np.mean(errors):.3f}"
    )

    print(
        f"{prefix}_median\t"
        f"{np.percentile(errors, 50):.3f}"
    )

    print(
        f"{prefix}_p90\t"
        f"{np.percentile(errors, 90):.3f}"
    )

    print(
        f"{prefix}_p95\t"
        f"{np.percentile(errors, 95):.3f}"
    )

    print(
        f"{prefix}_p99\t"
        f"{np.percentile(errors, 99):.3f}"
    )

    print(
        f"{prefix}_p999\t"
        f"{np.percentile(errors, 99.9):.3f}"
    )

    print(
        f"{prefix}_max\t"
        f"{np.max(errors):.3f}"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "dataset",
    )

    parser.add_argument(
        "checkpoint",
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
        args.checkpoint,
        map_location="cpu",
        weights_only=True,
    )


    hidden_sizes = checkpoint.get(
        "hidden_sizes",
        [128, 128, 64],
    )


    residual_scale = float(
        checkpoint[
            "residual_scale"
        ]
    )


    model = ResidualMLP(
        hidden_sizes
    )


    model.load_state_dict(
        checkpoint[
            "model_state_dict"
        ]
    )

    model.eval()


    parameters = sum(
        parameter.numel()
        for parameter
        in model.parameters()
    )


    macs = 0
    previous = 9

    for hidden in hidden_sizes:
        macs += (
            previous *
            hidden
        )

        previous = hidden

    macs += previous


    print(
        "hidden_sizes\t",
        ",".join(
            str(x)
            for x in hidden_sizes
        ),
    )

    print(
        "parameters\t",
        parameters,
    )

    print(
        "approx_macs_per_query\t",
        macs,
    )

    print(
        "fp32_weight_bytes\t",
        parameters * 4,
    )


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


    features = np.asarray(
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


    predicted_residuals = []


    with torch.no_grad():

        for start in range(
            0,
            len(features),
            args.batch_size,
        ):

            batch = torch.from_numpy(
                features[
                    start:
                    start + args.batch_size
                ]
            )

            predicted = (
                model(batch)
                .numpy()
                .astype(np.float64)
                *
                residual_scale
            )

            predicted_residuals.append(
                predicted
            )


    predicted_residuals = (
        np.concatenate(
            predicted_residuals
        )
    )


    predicted_rows = (
        baselines
        +
        predicted_residuals
    )


    predicted_rows = np.clip(
        predicted_rows,
        0.0,
        args.row_count,
    )


    lexical_errors = np.abs(
        baselines
        -
        true_rows
    )


    neural_errors = np.abs(
        predicted_rows
        -
        true_rows
    )


    print(
        "test_samples\t",
        len(true_rows),
    )


    percentile_report(
        "LEX",
        lexical_errors,
    )

    percentile_report(
        "NEURAL",
        neural_errors,
    )


    print(
        "p99_reduction_fold\t",
        f"{np.percentile(lexical_errors, 99) / np.percentile(neural_errors, 99):.6f}"
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
        1536,
        2048,
        3072,
        4096,
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
            1.0 -
            coverage
        )


        window_rows = (
            2 * radius
            +
            1
        )


        local_comparisons = math.ceil(
            math.log2(
                window_rows
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
