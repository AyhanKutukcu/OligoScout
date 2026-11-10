#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import numpy as np
import torch


BASE = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}


def prefix_fraction(
    sequence,
    length,
):
    value = 0

    for base in sequence[:length]:
        value = (
            value * 4
            +
            BASE[base]
        )

    return (
        value /
        float(
            4 ** length
        )
    )


def sequence_features(
    sequence,
):
    prefix_lengths = (
        4,
        8,
        12,
        16,
        21,
    )

    values = [
        prefix_fraction(
            sequence,
            length,
        )
        for length in prefix_lengths
    ]

    gc = sum(
        base in ("G", "C")
        for base in sequence
    ) / 21.0

    counts = [
        sequence.count(base)
        for base in "ACGT"
    ]

    entropy = 0.0

    for count in counts:

        if count == 0:
            continue

        probability = (
            count /
            21.0
        )

        entropy -= (
            probability
            *
            math.log2(
                probability
            )
        )

    maximum_run = 1
    current_run = 1

    for i in range(
        1,
        len(sequence),
    ):

        if (
            sequence[i]
            ==
            sequence[i - 1]
        ):
            current_run += 1
            maximum_run = max(
                maximum_run,
                current_run,
            )

        else:
            current_run = 1

    distinct = len(
        set(sequence)
    )

    values.extend(
        [
            gc,
            entropy / 2.0,
            maximum_run / 21.0,
            distinct / 4.0,
        ]
    )

    return np.asarray(
        values,
        dtype=np.float32,
    )


class Model(torch.nn.Module):
    def __init__(self):
        super().__init__()

        self.network = torch.nn.Sequential(
            torch.nn.Linear(9, 64),
            torch.nn.ReLU(),

            torch.nn.Linear(64, 32),
            torch.nn.ReLU(),

            torch.nn.Linear(32, 1),
        )

    def forward(self, x):
        return (
            self.network(x)
            .squeeze(-1)
        )


def float_literal(value):
    value = np.float32(
        value
    )

    text = format(
        float(value),
        ".9g",
    )

    if (
        "." not in text
        and
        "e" not in text.lower()
    ):
        text += ".0"

    return (
        text
        +
        "f"
    )


def write_array(
    handle,
    name,
    values,
):
    values = (
        values
        .detach()
        .cpu()
        .numpy()
        .astype(np.float32)
        .reshape(-1)
    )

    handle.write(
        f"inline constexpr "
        f"std::array<float, {len(values)}> "
        f"{name}{{{{\n"
    )

    for i, value in enumerate(
        values
    ):
        handle.write(
            "    "
            +
            float_literal(value)
        )

        if (
            i + 1
            !=
            len(values)
        ):
            handle.write(",")

        handle.write("\n")

    handle.write(
        "}};\n\n"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "checkpoint",
    )

    parser.add_argument(
        "output_header",
    )

    parser.add_argument(
        "golden_header",
    )

    args = parser.parse_args()


    checkpoint = torch.load(
        args.checkpoint,
        map_location="cpu",
        weights_only=True,
    )


    hidden_sizes = checkpoint[
        "hidden_sizes"
    ]

    if (
        list(hidden_sizes)
        !=
        [64, 32]
    ):
        raise RuntimeError(
            f"Unexpected architecture: "
            f"{hidden_sizes}"
        )


    model = Model()

    model.load_state_dict(
        checkpoint[
            "model_state_dict"
        ]
    )

    model.eval()


    output = Path(
        args.output_header
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    state = checkpoint[
        "model_state_dict"
    ]


    with output.open("w") as handle:

        handle.write(
            "#pragma once\n\n"
            "#include <array>\n"
            "#include <cstddef>\n\n"
            "namespace primerpair::"
            "neural_window_weights {\n\n"
        )

        handle.write(
            "inline constexpr "
            "std::size_t input_size = 9;\n"
            "inline constexpr "
            "std::size_t hidden1_size = 64;\n"
            "inline constexpr "
            "std::size_t hidden2_size = 32;\n"
            "inline constexpr "
            "float residual_scale = 10000.0f;\n\n"
        )

        write_array(
            handle,
            "layer1_weight",
            state[
                "network.0.weight"
            ],
        )

        write_array(
            handle,
            "layer1_bias",
            state[
                "network.0.bias"
            ],
        )

        write_array(
            handle,
            "layer2_weight",
            state[
                "network.2.weight"
            ],
        )

        write_array(
            handle,
            "layer2_bias",
            state[
                "network.2.bias"
            ],
        )

        write_array(
            handle,
            "output_weight",
            state[
                "network.4.weight"
            ],
        )

        write_array(
            handle,
            "output_bias",
            state[
                "network.4.bias"
            ],
        )

        handle.write(
            "}  // namespace "
            "primerpair::neural_window_weights\n"
        )


    golden_sequences = [
        "AAAAAAAAAAAAAAAAAAAAA",
        "ACGTACGTACGTACGTACGTA",
        "CGCGGAAGCAAAGTGACTTCC",
        "GAAATATAGGTATCAACGGAG",
        "CTGAATGGAATTCCTCCGATC",
        "CAAATGACGATGTCCTTGGGT",
        "GGGTTTTTTTTACACACACGT",
        "TGCATGCATGCATGCATGCAT",
    ]


    features = np.stack(
        [
            sequence_features(
                sequence
            )
            for sequence
            in golden_sequences
        ]
    )


    with torch.no_grad():

        predicted = (
            model(
                torch.from_numpy(
                    features
                )
            )
            .numpy()
            *
            float(
                checkpoint[
                    "residual_scale"
                ]
            )
        )


    golden = Path(
        args.golden_header
    )

    golden.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    with golden.open("w") as handle:

        handle.write(
            "#pragma once\n\n"
            "#include <array>\n"
            "#include <string_view>\n\n"
            "namespace primerpair::"
            "neural_window_golden {\n\n"
        )

        handle.write(
            "inline constexpr "
            "std::array<std::string_view, 8> "
            "kmers{{\n"
        )

        for sequence in golden_sequences:
            handle.write(
                f'    "{sequence}",\n'
            )

        handle.write(
            "}};\n\n"
        )

        handle.write(
            "inline constexpr "
            "std::array<float, 8> "
            "residual_rows{{\n"
        )

        for value in predicted:
            handle.write(
                "    "
                +
                float_literal(
                    value
                )
                +
                ",\n"
            )

        handle.write(
            "}};\n\n"
            "}  // namespace "
            "primerpair::neural_window_golden\n"
        )


    parameter_count = sum(
        tensor.numel()
        for tensor
        in state.values()
    )


    print(
        "hidden_sizes\t64,32"
    )

    print(
        "parameters\t",
        parameter_count,
    )

    print(
        "fp32_bytes\t",
        parameter_count * 4,
    )

    print(
        "weights_header\t",
        output,
    )

    print(
        "golden_header\t",
        golden,
    )

    print(
        "ALL_CHECKS\tYES"
    )


if __name__ == "__main__":
    main()
