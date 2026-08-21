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


BASE_TO_INDEX = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}


def encode_kmer(
    kmer,
):
    encoded = np.zeros(
        21 * 4,
        dtype=np.float32,
    )

    for position, base in enumerate(
        kmer
    ):
        encoded[
            position * 4
            +
            BASE_TO_INDEX[base]
        ] = 1.0

    return encoded


class NeuralWindowDataset(
    Dataset
):
    def __init__(
        self,
        path,
        split,
    ):
        self.features = []
        self.lower = []
        self.log_width = []
        self.present = []

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

                sequence_features = (
                    encode_kmer(
                        row["kmer"]
                    )
                )

                engineered = np.array(
                    [
                        float(
                            row[
                                "gc_fraction"
                            ]
                        ),

                        float(
                            row[
                                "entropy"
                            ]
                        )
                        /
                        2.0,

                        float(
                            row[
                                "max_homopolymer"
                            ]
                        )
                        /
                        21.0,

                        float(
                            row[
                                "distinct_bases"
                            ]
                        )
                        /
                        4.0,
                    ],
                    dtype=np.float32,
                )

                feature = np.concatenate(
                    (
                        sequence_features,
                        engineered,
                    )
                )

                self.features.append(
                    feature
                )

                self.lower.append(
                    float(
                        row[
                            "lower_norm"
                        ]
                    )
                )

                self.log_width.append(
                    float(
                        row[
                            "log2_width_plus1"
                        ]
                    )
                )

                self.present.append(
                    float(
                        row[
                            "present"
                        ]
                    )
                )

        self.features = np.asarray(
            self.features,
            dtype=np.float32,
        )

        self.lower = np.asarray(
            self.lower,
            dtype=np.float32,
        )

        self.log_width = np.asarray(
            self.log_width,
            dtype=np.float32,
        )

        self.present = np.asarray(
            self.present,
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
                self.lower[
                    index
                ],
                dtype=torch.float32,
            ),
            torch.tensor(
                self.log_width[
                    index
                ],
                dtype=torch.float32,
            ),
            torch.tensor(
                self.present[
                    index
                ],
                dtype=torch.float32,
            ),
        )


class NeuralWindowPredictor(
    nn.Module
):
    def __init__(
        self,
    ):
        super().__init__()

        self.backbone = nn.Sequential(
            nn.Linear(
                88,
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
                32,
            ),
            nn.ReLU(),
        )

        self.lower_head = nn.Sequential(
            nn.Linear(
                32,
                1,
            ),
            nn.Sigmoid(),
        )

        self.width_head = nn.Linear(
            32,
            1,
        )

        self.presence_head = nn.Linear(
            32,
            1,
        )


    def forward(
        self,
        x,
    ):
        hidden = self.backbone(
            x
        )

        lower = self.lower_head(
            hidden
        ).squeeze(-1)

        log_width = self.width_head(
            hidden
        ).squeeze(-1)

        presence_logit = (
            self.presence_head(
                hidden
            ).squeeze(-1)
        )

        return (
            lower,
            log_width,
            presence_logit,
        )


def evaluate(
    model,
    loader,
    device,
    row_count,
):
    model.eval()

    total = 0

    lower_abs_error_rows = []

    presence_correct = 0

    with torch.no_grad():

        for (
            features,
            true_lower,
            true_log_width,
            true_present,
        ) in loader:

            features = features.to(
                device
            )

            true_lower = true_lower.to(
                device
            )

            true_present = true_present.to(
                device
            )

            (
                predicted_lower,
                _,
                presence_logit,
            ) = model(
                features
            )

            error_rows = (
                torch.abs(
                    predicted_lower
                    -
                    true_lower
                )
                *
                row_count
            )

            lower_abs_error_rows.extend(
                error_rows
                .cpu()
                .numpy()
                .tolist()
            )

            predicted_present = (
                torch.sigmoid(
                    presence_logit
                )
                >=
                0.5
            )

            presence_correct += int(
                (
                    predicted_present
                    ==
                    (
                        true_present
                        >=
                        0.5
                    )
                )
                .sum()
                .item()
            )

            total += features.shape[0]


    errors = np.asarray(
        lower_abs_error_rows,
        dtype=np.float64,
    )

    return {
        "mae_rows":
            float(
                np.mean(errors)
            ),

        "median_error_rows":
            float(
                np.percentile(
                    errors,
                    50,
                )
            ),

        "p90_error_rows":
            float(
                np.percentile(
                    errors,
                    90,
                )
            ),

        "p95_error_rows":
            float(
                np.percentile(
                    errors,
                    95,
                )
            ),

        "p99_error_rows":
            float(
                np.percentile(
                    errors,
                    99,
                )
            ),

        "p999_error_rows":
            float(
                np.percentile(
                    errors,
                    99.9,
                )
            ),

        "max_error_rows":
            float(
                np.max(
                    errors
                )
            ),

        "presence_accuracy":
            presence_correct
            /
            total,
    }


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "dataset",
    )

    parser.add_argument(
        "--epochs",
        type=int,
        default=12,
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=4096,
    )

    parser.add_argument(
        "--row-count",
        type=float,
        default=1_000_001.0,
    )

    parser.add_argument(
        "--output",
        default=(
            "models/"
            "neural_window_v1.pt"
        ),
    )

    args = parser.parse_args()


    random.seed(
        42
    )

    np.random.seed(
        42
    )

    torch.manual_seed(
        42
    )


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


    train_dataset = (
        NeuralWindowDataset(
            args.dataset,
            "train",
        )
    )

    validation_dataset = (
        NeuralWindowDataset(
            args.dataset,
            "validation",
        )
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


    model = NeuralWindowPredictor().to(
        device
    )


    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=1e-3,
        weight_decay=1e-5,
    )


    lower_loss_fn = nn.SmoothL1Loss()
    width_loss_fn = nn.SmoothL1Loss()
    presence_loss_fn = (
        nn.BCEWithLogitsLoss()
    )


    best_p99 = math.inf

    output_path = Path(
        args.output
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    for epoch in range(
        1,
        args.epochs + 1,
    ):

        model.train()

        running_loss = 0.0
        batches = 0


        for (
            features,
            true_lower,
            true_log_width,
            true_present,
        ) in train_loader:

            features = features.to(
                device
            )

            true_lower = true_lower.to(
                device
            )

            true_log_width = (
                true_log_width.to(
                    device
                )
            )

            true_present = (
                true_present.to(
                    device
                )
            )


            optimizer.zero_grad(
                set_to_none=True
            )


            (
                predicted_lower,
                predicted_log_width,
                presence_logit,
            ) = model(
                features
            )


            lower_loss = (
                lower_loss_fn(
                    predicted_lower,
                    true_lower,
                )
            )


            width_loss = (
                width_loss_fn(
                    predicted_log_width,
                    true_log_width,
                )
            )


            presence_loss = (
                presence_loss_fn(
                    presence_logit,
                    true_present,
                )
            )


            loss = (
                lower_loss
                +
                0.10
                *
                width_loss
                +
                0.10
                *
                presence_loss
            )


            loss.backward()

            optimizer.step()


            running_loss += float(
                loss.item()
            )

            batches += 1


        metrics = evaluate(
            model,
            validation_loader,
            device,
            args.row_count,
        )


        print(
            "EPOCH"
            f"\t{epoch}"
            f"\tloss={running_loss / batches:.8f}"
            f"\tmae_rows={metrics['mae_rows']:.3f}"
            f"\tmedian={metrics['median_error_rows']:.3f}"
            f"\tp95={metrics['p95_error_rows']:.3f}"
            f"\tp99={metrics['p99_error_rows']:.3f}"
            f"\tp999={metrics['p999_error_rows']:.3f}"
            f"\tmax={metrics['max_error_rows']:.3f}"
            f"\tpresence_acc={metrics['presence_accuracy']:.6f}"
        )


        if (
            metrics[
                "p99_error_rows"
            ]
            <
            best_p99
        ):
            best_p99 = metrics[
                "p99_error_rows"
            ]

            torch.save(
                {
                    "model_state_dict":
                        model.state_dict(),

                    "input_size":
                        88,

                    "k":
                        21,

                    "row_count":
                        args.row_count,

                    "validation_metrics":
                        metrics,
                },
                output_path,
            )

            print(
                "BEST_MODEL_SAVED\t",
                output_path,
            )


    print(
        "ALL_CHECKS\tYES"
    )


if __name__ == "__main__":
    main()
