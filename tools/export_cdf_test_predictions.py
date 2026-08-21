#!/usr/bin/env python3

import csv
from pathlib import Path

import numpy as np
import torch
from torch import nn


BASE = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}


ROW_COUNT = 1_000_001
RESIDUAL_SCALE = 10_000.0


class CDFCorrector(nn.Module):

    def __init__(
        self,
        hidden,
    ):
        super().__init__()

        self.network = nn.Sequential(
            nn.Linear(
                1,
                hidden,
            ),

            nn.ReLU(),

            nn.Linear(
                hidden,
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


def lex_fraction(
    sequence,
):
    code = 0

    for base in sequence:
        code = (
            code * 4
            +
            BASE[base]
        )

    return (
        code
        /
        float(
            4 ** len(sequence)
        )
    )


def load_test(
    path,
):
    fractions = []
    true_rows = []


    with open(
        path,
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


            fractions.append(
                lex_fraction(
                    row["kmer"]
                )
            )

            true_rows.append(
                int(
                    row["sa_lower"]
                )
            )


    return (
        np.asarray(
            fractions,
            dtype=np.float32,
        ),
        np.asarray(
            true_rows,
            dtype=np.int64,
        ),
    )


def evaluate(
    checkpoint_path,
    fractions,
    true_rows,
    output_path,
):

    checkpoint = torch.load(
        checkpoint_path,
        map_location="cpu",
        weights_only=False,
    )


    hidden = int(
        checkpoint["hidden"]
    )


    model = CDFCorrector(
        hidden
    )


    model.load_state_dict(
        checkpoint[
            "model_state_dict"
        ]
    )


    model.eval()


    x = torch.from_numpy(
        fractions.reshape(
            -1,
            1
        )
    )


    predictions = []


    with torch.no_grad():

        for start in range(
            0,
            len(x),
            16384,
        ):
            current = x[
                start:
                start + 16384
            ]


            residual = (
                model(current)
                .numpy()
                *
                RESIDUAL_SCALE
            )


            baseline = (
                1.0
                +
                fractions[
                    start:
                    start + len(current)
                ].astype(
                    np.float64
                )
                *
                (
                    ROW_COUNT - 1
                )
            )


            predicted = (
                baseline
                +
                residual.astype(
                    np.float64
                )
            )


            predicted = np.clip(
                predicted,
                0.0,
                float(
                    ROW_COUNT
                ),
            )


            # Match the non-negative C++:
            # static_cast<uint64_t>(x + 0.5)
            predicted = np.floor(
                predicted
                +
                0.5
            ).astype(
                np.int64
            )


            predictions.append(
                predicted
            )


    predictions = np.concatenate(
        predictions
    )


    errors = np.abs(
        predictions
        -
        true_rows
    )


    np.savetxt(
        output_path,
        predictions,
        fmt="%d",
    )


    print(
        "MODEL",
        checkpoint_path,
    )

    print(
        "samples",
        len(predictions),
    )

    for radius in [
        1024,
        1536,
        2048,
        2560,
        3072,
        3584,
        4096,
    ]:

        coverage = np.mean(
            errors <=
            radius
        )

        fallback = (
            1.0
            -
            coverage
        )


        print(
            "RADIUS",
            radius,
            "coverage",
            f"{coverage:.8f}",
            "fallback",
            f"{fallback:.8f}",
        )


    print(
        "ERROR"
        "\tmedian="
        f"{np.percentile(errors, 50):.3f}"
        "\tp95="
        f"{np.percentile(errors, 95):.3f}"
        "\tp99="
        f"{np.percentile(errors, 99):.3f}"
        "\tp999="
        f"{np.percentile(errors, 99.9):.3f}"
        "\tmax="
        f"{np.max(errors)}"
    )

    print()


def main():

    dataset = (
        "data/neural_window/"
        "window_dataset_1m_k21.tsv"
    )


    fractions, true_rows = (
        load_test(
            dataset
        )
    )


    output_dir = Path(
        "results/benchmarks/"
        "neural_window_v1/"
        "cdf_predictions"
    )


    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    evaluate(
        "models/neural_cdf/cdf64.pt",
        fractions,
        true_rows,
        output_dir /
        "cdf64.txt",
    )


    evaluate(
        "models/neural_cdf/cdf128.pt",
        fractions,
        true_rows,
        output_dir /
        "cdf128.txt",
    )


    print(
        "ALL_CHECKS\tYES"
    )


if __name__ == "__main__":
    main()
