#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import torch


MODELS = [
    (
        "cdf64_v1",
        "models/neural_cdf/cdf64.pt",
        64,
    ),
    (
        "cdf128_v1",
        "models/neural_cdf/cdf128.pt",
        128,
    ),
]


def cpp_float(value):
    value = float(
        np.float32(value)
    )

    text = format(
        value,
        ".9g",
    )

    if (
        "." not in text
        and
        "e" not in text.lower()
    ):
        text += ".0"

    return text + "f"


def format_array(
    values,
    indent="        ",
    per_line=6,
):
    values = list(
        values
    )

    lines = []


    for start in range(
        0,
        len(values),
        per_line,
    ):
        chunk = values[
            start:
            start + per_line
        ]


        lines.append(
            indent
            +
            ", ".join(
                cpp_float(value)
                for value
                in chunk
            )
        )


    return ",\n".join(
        lines
    )


def export_model(
    namespace,
    checkpoint_path,
    expected_hidden,
):
    checkpoint = torch.load(
        checkpoint_path,
        map_location="cpu",
        weights_only=False,
    )


    hidden = int(
        checkpoint["hidden"]
    )


    if (
        hidden !=
        expected_hidden
    ):
        raise RuntimeError(
            f"{checkpoint_path}: "
            f"hidden={hidden}, "
            f"expected={expected_hidden}"
        )


    state = checkpoint[
        "model_state_dict"
    ]


    hidden_weight = (
        state["network.0.weight"]
        .detach()
        .cpu()
        .numpy()
        .reshape(-1)
    )


    hidden_bias = (
        state["network.0.bias"]
        .detach()
        .cpu()
        .numpy()
        .reshape(-1)
    )


    output_weight = (
        state["network.2.weight"]
        .detach()
        .cpu()
        .numpy()
        .reshape(-1)
    )


    output_bias = float(
        state["network.2.bias"]
        .detach()
        .cpu()
        .numpy()
        .reshape(-1)[0]
    )


    if (
        len(hidden_weight) !=
        hidden
        or
        len(hidden_bias) !=
        hidden
        or
        len(output_weight) !=
        hidden
    ):
        raise RuntimeError(
            f"{checkpoint_path}: "
            "unexpected tensor shape."
        )


    parameter_count = int(
        checkpoint[
            "parameter_count"
        ]
    )


    residual_scale = float(
        checkpoint[
            "residual_scale"
        ]
    )


    return f'''
namespace {namespace} {{

inline constexpr
std::size_t hidden_size =
    {hidden};


inline constexpr
std::size_t parameter_count =
    {parameter_count};


inline constexpr
float residual_scale =
    {cpp_float(residual_scale)};


inline constexpr
std::array<float, hidden_size>
hidden_weight{{
{format_array(hidden_weight)}
}};


inline constexpr
std::array<float, hidden_size>
hidden_bias{{
{format_array(hidden_bias)}
}};


inline constexpr
std::array<float, hidden_size>
output_weight{{
{format_array(output_weight)}
}};


inline constexpr
float output_bias =
    {cpp_float(output_bias)};

}}  // namespace {namespace}

'''


def main():
    output = Path(
        "include/primerpair/generated/"
        "neural_cdf_v1_weights.hpp"
    )


    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    parts = [
        "#pragma once\n\n",
        "#include <array>\n",
        "#include <cstddef>\n\n",
        "namespace primerpair::generated {\n\n",
    ]


    for (
        namespace,
        checkpoint_path,
        hidden,
    ) in MODELS:
        checkpoint = Path(
            checkpoint_path
        )


        if not checkpoint.exists():
            raise FileNotFoundError(
                checkpoint
            )


        parts.append(
            export_model(
                namespace,
                checkpoint,
                hidden,
            )
        )


    parts.append(
        "}  // namespace primerpair::generated\n"
    )


    output.write_text(
        "".join(parts)
    )


    print(
        "output\t",
        output
    )

    print(
        "cdf64_params\t193"
    )

    print(
        "cdf128_params\t385"
    )

    print(
        "ALL_CHECKS\tYES"
    )


if __name__ == "__main__":
    main()
