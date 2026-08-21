#!/usr/bin/env python3

import argparse
import csv
import math
import random
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset


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


def load_dataset(
    path,
    row_count,
    residual_scale,
):
    data = {
        "train_x": [],
        "train_y": [],
        "validation_x": [],
        "validation_y": [],
    }

    with open(
        path,
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:

            split = row["split"]

            if split not in (
                "train",
                "validation",
            ):
                continue

            features = make_features(
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
            ) / residual_scale

            if split == "train":

                data["train_x"].append(
                    features
                )

                data["train_y"].append(
                    residual
                )

            else:

                data["validation_x"].append(
                    features
                )

                data["validation_y"].append(
                    residual
                )

    for key in data:
        data[key] = np.asarray(
            data[key],
            dtype=np.float32,
        )

    print(
        "train_loaded\t",
        len(data["train_x"]),
    )

    print(
        "validation_loaded\t",
        len(
            data["validation_x"]
        ),
    )

    return data


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


def parameter_count(
    model,
):
    return sum(
        parameter.numel()
        for parameter
        in model.parameters()
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

        for features, truth in loader:

            features = features.to(
                device
            )

            truth = truth.to(
                device
            )

            prediction = model(
                features
            )

            error = (
                torch.abs(
                    prediction - truth
                )
                *
                residual_scale
            )

            errors.extend(
                error.cpu()
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
                np.max(errors)
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
        default=8,
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
        "--output-dir",
        default="models/neural_window_sweep",
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
        device,
    )

    data = load_dataset(
        args.dataset,
        args.row_count,
        args.residual_scale,
    )

    train_dataset = TensorDataset(
        torch.from_numpy(
            data["train_x"]
        ),
        torch.from_numpy(
            data["train_y"]
        ),
    )

    validation_dataset = TensorDataset(
        torch.from_numpy(
            data["validation_x"]
        ),
        torch.from_numpy(
            data["validation_y"]
        ),
    )

    validation_loader = DataLoader(
        validation_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=0,
    )

    architectures = [
        (
            "linear",
            [],
        ),

        (
            "tiny16",
            [16],
        ),

        (
            "tiny32",
            [32],
        ),

        (
            "tiny32x16",
            [32, 16],
        ),

        (
            "tiny64x32",
            [64, 32],
        ),

        (
            "v2_128x128x64",
            [128, 128, 64],
        ),
    ]

    output_dir = Path(
        args.output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "SWEEP_START"
    )

    for (
        name,
        hidden_sizes,
    ) in architectures:

        torch.manual_seed(42)

        model = ResidualMLP(
            hidden_sizes
        ).to(
            device
        )

        params = parameter_count(
            model
        )

        print(
            "MODEL_START"
            f"\t{name}"
            f"\tparams={params}"
            f"\thidden={hidden_sizes}"
        )

        train_loader = DataLoader(
            train_dataset,
            batch_size=args.batch_size,
            shuffle=True,
            num_workers=0,
            generator=(
                torch.Generator()
                .manual_seed(42)
            ),
        )

        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=5e-4,
            weight_decay=1e-6,
        )

        loss_fn = nn.SmoothL1Loss()

        best_p99 = math.inf
        best_metrics = None

        output_path = (
            output_dir
            /
            f"{name}.pt"
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
                truth,
            ) in train_loader:

                features = features.to(
                    device
                )

                truth = truth.to(
                    device
                )

                optimizer.zero_grad(
                    set_to_none=True
                )

                prediction = model(
                    features
                )

                loss = loss_fn(
                    prediction,
                    truth,
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
                "MODEL_EPOCH"
                f"\t{name}"
                f"\t{epoch}"
                f"\tloss={total_loss / batches:.8f}"
                f"\tmae={metrics['mae']:.3f}"
                f"\tmedian={metrics['median']:.3f}"
                f"\tp95={metrics['p95']:.3f}"
                f"\tp99={metrics['p99']:.3f}"
                f"\tp999={metrics['p999']:.3f}"
            )

            if (
                metrics["p99"]
                <
                best_p99
            ):
                best_p99 = metrics[
                    "p99"
                ]

                best_metrics = metrics

                torch.save(
                    {
                        "model_state_dict":
                            model.state_dict(),

                        "hidden_sizes":
                            hidden_sizes,

                        "parameter_count":
                            params,

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
                    output_path,
                )

        print(
            "MODEL_RESULT"
            f"\t{name}"
            f"\tparams={params}"
            f"\tmae={best_metrics['mae']:.3f}"
            f"\tmedian={best_metrics['median']:.3f}"
            f"\tp95={best_metrics['p95']:.3f}"
            f"\tp99={best_metrics['p99']:.3f}"
            f"\tp999={best_metrics['p999']:.3f}"
            f"\tmax={best_metrics['max']:.3f}"
        )

    print(
        "ALL_CHECKS\tYES"
    )


if __name__ == "__main__":
    main()
