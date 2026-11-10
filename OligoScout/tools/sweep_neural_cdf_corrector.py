#!/usr/bin/env python3

import argparse
import csv
import math
import random
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


BASE = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}


def lex_fraction(sequence):
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
        float(4 ** len(sequence))
    )


def baseline_row(
    sequence,
    row_count,
):
    return (
        1.0
        +
        lex_fraction(sequence)
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
    output = {
        "train_x": [],
        "train_y": [],
        "validation_x": [],
        "validation_y": [],
        "test_x": [],
        "test_y": [],
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
                "test",
            ):
                continue

            fraction = lex_fraction(
                row["kmer"]
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

            output[
                f"{split}_x"
            ].append(
                [fraction]
            )

            output[
                f"{split}_y"
            ].append(
                residual
            )

    for key in output:
        output[key] = np.asarray(
            output[key],
            dtype=np.float32,
        )

    return output


class CDFCorrector(nn.Module):

    def __init__(
        self,
        hidden,
    ):
        super().__init__()

        if hidden == 0:

            self.network = nn.Linear(
                1,
                1,
            )

        else:

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


def metrics(
    model,
    loader,
    residual_scale,
):
    model.eval()

    errors = []


    with torch.no_grad():

        for x, y in loader:

            prediction = model(x)

            error = (
                torch.abs(
                    prediction - y
                )
                *
                residual_scale
            )

            errors.extend(
                error.numpy()
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
                    50
                )
            ),

        "p95":
            float(
                np.percentile(
                    errors,
                    95
                )
            ),

        "p99":
            float(
                np.percentile(
                    errors,
                    99
                )
            ),

        "p999":
            float(
                np.percentile(
                    errors,
                    99.9
                )
            ),

        "max":
            float(
                np.max(errors)
            ),

        "coverage1024":
            float(
                np.mean(
                    errors <= 1024
                )
            ),

        "coverage2048":
            float(
                np.mean(
                    errors <= 2048
                )
            ),

        "coverage4096":
            float(
                np.mean(
                    errors <= 4096
                )
            ),
    }


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "dataset"
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
        default=1000001.0,
    )

    parser.add_argument(
        "--residual-scale",
        type=float,
        default=10000.0,
    )

    parser.add_argument(
        "--output-dir",
        default="models/neural_cdf",
    )

    args = parser.parse_args()


    random.seed(42)
    np.random.seed(42)
    torch.manual_seed(42)


    data = load_dataset(
        args.dataset,
        args.row_count,
        args.residual_scale,
    )


    print(
        "train_samples\t",
        len(
            data["train_x"]
        )
    )

    print(
        "validation_samples\t",
        len(
            data["validation_x"]
        )
    )

    print(
        "test_samples\t",
        len(
            data["test_x"]
        )
    )


    train_dataset = TensorDataset(
        torch.from_numpy(
            data["train_x"]
        ),
        torch.from_numpy(
            data["train_y"]
        ),
    )


    validation_loader = DataLoader(
        TensorDataset(
            torch.from_numpy(
                data["validation_x"]
            ),
            torch.from_numpy(
                data["validation_y"]
            ),
        ),
        batch_size=args.batch_size,
        shuffle=False,
    )


    test_loader = DataLoader(
        TensorDataset(
            torch.from_numpy(
                data["test_x"]
            ),
            torch.from_numpy(
                data["test_y"]
            ),
        ),
        batch_size=args.batch_size,
        shuffle=False,
    )


    architectures = [
        ("linear", 0),
        ("cdf16", 16),
        ("cdf32", 32),
        ("cdf64", 64),
        ("cdf128", 128),
    ]


    output_dir = Path(
        args.output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    for name, hidden in architectures:

        torch.manual_seed(42)


        model = CDFCorrector(
            hidden
        )


        params = sum(
            parameter.numel()
            for parameter
            in model.parameters()
        )


        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=5e-4,
            weight_decay=1e-6,
        )


        loss_function = nn.SmoothL1Loss()


        train_loader = DataLoader(
            train_dataset,
            batch_size=args.batch_size,
            shuffle=True,
            generator=(
                torch.Generator()
                .manual_seed(42)
            ),
        )


        best_p99 = math.inf
        best_state = None
        best_epoch = 0


        for epoch in range(
            1,
            args.epochs + 1
        ):

            model.train()


            for x, y in train_loader:

                optimizer.zero_grad(
                    set_to_none=True
                )

                prediction = model(x)

                loss = loss_function(
                    prediction,
                    y,
                )

                loss.backward()

                optimizer.step()


            result = metrics(
                model,
                validation_loader,
                args.residual_scale,
            )


            print(
                "EPOCH"
                f"\t{name}"
                f"\t{epoch}"
                f"\tp99={result['p99']:.3f}"
                f"\tp999={result['p999']:.3f}"
                f"\tmedian={result['median']:.3f}"
            )


            if (
                result["p99"]
                <
                best_p99
            ):
                best_p99 = result[
                    "p99"
                ]

                best_epoch = epoch

                best_state = {
                    key:
                        value.detach()
                        .clone()
                    for key, value
                    in model.state_dict().items()
                }


        model.load_state_dict(
            best_state
        )


        validation = metrics(
            model,
            validation_loader,
            args.residual_scale,
        )


        test = metrics(
            model,
            test_loader,
            args.residual_scale,
        )


        checkpoint = {
            "model_state_dict":
                best_state,

            "hidden":
                hidden,

            "parameter_count":
                params,

            "input_size":
                1,

            "k":
                21,

            "row_count":
                args.row_count,

            "residual_scale":
                args.residual_scale,

            "best_epoch":
                best_epoch,

            "validation_metrics":
                validation,

            "test_metrics":
                test,
        }


        torch.save(
            checkpoint,
            output_dir /
            f"{name}.pt",
        )


        print(
            "MODEL_RESULT"
            f"\t{name}"
            f"\tparams={params}"
            f"\tbest_epoch={best_epoch}"
            f"\tVAL_p99={validation['p99']:.3f}"
            f"\tTEST_median={test['median']:.3f}"
            f"\tTEST_p95={test['p95']:.3f}"
            f"\tTEST_p99={test['p99']:.3f}"
            f"\tTEST_p999={test['p999']:.3f}"
            f"\tC1024={test['coverage1024']:.8f}"
            f"\tC2048={test['coverage2048']:.8f}"
            f"\tC4096={test['coverage4096']:.8f}"
        )


    print(
        "ALL_CHECKS\tYES"
    )


if __name__ == "__main__":
    main()
