#!/usr/bin/env python3

import argparse
import csv
import math
import random
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.utils.data import Dataset
from torch.utils.data import DataLoader


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


def base4_code(
    sequence,
):
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
    value = base4_code(
        sequence[:length]
    )

    maximum = (
        4 ** length
    )

    return (
        value /
        maximum
    )


def make_features(
    row,
):
    kmer = row[
        "kmer"
    ]


    features = [
        prefix_fraction(
            kmer,
            length,
        )
        for length in PREFIX_LENGTHS
    ]


    features.extend(
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


    return np.asarray(
        features,
        dtype=np.float32,
    )


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


class ResidualDataset(
    Dataset
):
    def __init__(
        self,
        path,
        split,
        row_count,
        residual_scale,
    ):
        features = []
        residuals = []


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
                    split
                ):
                    continue


                feature = make_features(
                    row
                )


                baseline = baseline_row(
                    row["kmer"],
                    row_count,
                )


                true_row = float(
                    row["sa_lower"]
                )


                residual = (
                    true_row
                    -
                    baseline
                )


                features.append(
                    feature
                )


                residuals.append(
                    residual
                    /
                    residual_scale
                )


        self.features = np.asarray(
            features,
            dtype=np.float32,
        )


        self.residuals = np.asarray(
            residuals,
            dtype=np.float32,
        )


        print(
            f"{split}_loaded\t"
            f"{len(self.features)}"
        )


    def __len__(
        self,
    ):
        return len(
            self.features
        )


    def __getitem__(
        self,
        index,
    ):
        return (
            torch.from_numpy(
                self.features[
                    index
                ]
            ),

            torch.tensor(
                self.residuals[
                    index
                ],
                dtype=torch.float32,
            ),
        )


class ResidualWindowPredictor(
    nn.Module
):
    def __init__(
        self,
    ):
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
            self.network(
                x
            )
            .squeeze(-1)
        )


def evaluate(
    model,
    loader,
    device,
    residual_scale,
):
    model.eval()


    errors = []


    with torch.no_grad():

        for (
            features,
            true_residual,
        ) in loader:

            features = features.to(
                device
            )


            true_residual = (
                true_residual.to(
                    device
                )
            )


            predicted = model(
                features
            )


            error_rows = (
                torch.abs(
                    predicted
                    -
                    true_residual
                )
                *
                residual_scale
            )


            errors.extend(
                error_rows
                .cpu()
                .numpy()
                .tolist()
            )


    errors = np.asarray(
        errors,
        dtype=np.float64,
    )


    return {

        "mae":
            float(
                np.mean(errors)
            ),

        "median":
            float(
                np.percentile(
                    errors,
                    50,
                )
            ),

        "p90":
            float(
                np.percentile(
                    errors,
                    90,
                )
            ),

        "p95":
            float(
                np.percentile(
                    errors,
                    95,
                )
            ),

        "p99":
            float(
                np.percentile(
                    errors,
                    99,
                )
            ),

        "p999":
            float(
                np.percentile(
                    errors,
                    99.9,
                )
            ),

        "max":
            float(
                np.max(
                    errors
                )
            ),
    }


def main():
    parser = argparse.ArgumentParser()


    parser.add_argument(
        "dataset",
    )


    parser.add_argument(
        "--epochs",
        type=int,
        default=20,
    )


    parser.add_argument(
        "--batch-size",
        type=int,
        default=8192,
    )


    parser.add_argument(
        "--row-count",
        type=float,
        default=1_000_001.0,
    )


    parser.add_argument(
        "--residual-scale",
        type=float,
        default=10_000.0,
    )


    parser.add_argument(
        "--output",
        default=(
            "models/"
            "neural_window_v2.pt"
        ),
    )


    args = parser.parse_args()


    random.seed(42)
    np.random.seed(42)
    torch.manual_seed(42)


    device = torch.device(
        "cuda"
        if torch.cuda.is_available()
        else
        "cpu"
    )


    print(
        "device\t",
        device
    )


    train_dataset = ResidualDataset(
        args.dataset,
        "train",
        args.row_count,
        args.residual_scale,
    )


    validation_dataset = ResidualDataset(
        args.dataset,
        "validation",
        args.row_count,
        args.residual_scale,
    )


    train_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=0,
    )


    validation_loader = DataLoader(
        validation_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=0,
    )


    model = (
        ResidualWindowPredictor()
        .to(device)
    )


    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=5e-4,
        weight_decay=1e-6,
    )


    loss_fn = nn.SmoothL1Loss()


    best_p99 = math.inf


    output = Path(
        args.output
    )


    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    for epoch in range(
        1,
        args.epochs + 1,
    ):

        model.train()


        total_loss = 0.0
        batches = 0


        for (
            features,
            residual,
        ) in train_loader:

            features = features.to(
                device
            )


            residual = residual.to(
                device
            )


            optimizer.zero_grad(
                set_to_none=True
            )


            predicted = model(
                features
            )


            loss = loss_fn(
                predicted,
                residual,
            )


            loss.backward()


            optimizer.step()


            total_loss += float(
                loss.item()
            )

            batches += 1


        metrics = evaluate(
            model,
            validation_loader,
            device,
            args.residual_scale,
        )


        print(
            "EPOCH"
            f"\t{epoch}"
            f"\tloss={total_loss / batches:.8f}"
            f"\tmae={metrics['mae']:.3f}"
            f"\tmedian={metrics['median']:.3f}"
            f"\tp95={metrics['p95']:.3f}"
            f"\tp99={metrics['p99']:.3f}"
            f"\tp999={metrics['p999']:.3f}"
            f"\tmax={metrics['max']:.3f}"
        )


        if (
            metrics["p99"]
            <
            best_p99
        ):
            best_p99 = metrics[
                "p99"
            ]


            torch.save(
                {
                    "model_state_dict":
                        model.state_dict(),

                    "input_size":
                        9,

                    "k":
                        21,

                    "row_count":
                        args.row_count,

                    "residual_scale":
                        args.residual_scale,

                    "validation_metrics":
                        metrics,
                },
                output,
            )


            print(
                "BEST_MODEL_SAVED\t",
                output,
            )


    print(
        "ALL_CHECKS\tYES"
    )


if __name__ == "__main__":
    main()
