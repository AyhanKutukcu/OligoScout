#pragma once

#include <array>
#include <string_view>

namespace primerpair::neural_window_golden {

inline constexpr std::array<std::string_view, 8> kmers{{
    "AAAAAAAAAAAAAAAAAAAAA",
    "ACGTACGTACGTACGTACGTA",
    "CGCGGAAGCAAAGTGACTTCC",
    "GAAATATAGGTATCAACGGAG",
    "CTGAATGGAATTCCTCCGATC",
    "CAAATGACGATGTCCTTGGGT",
    "GGGTTTTTTTTACACACACGT",
    "TGCATGCATGCATGCATGCAT",
}};

inline constexpr std::array<float, 8> residual_rows{{
    -1588.58716f,
    3407.802f,
    1227.55383f,
    195.017761f,
    1310.65967f,
    -435.600281f,
    -3277.45728f,
    -1419.14185f,
}};

}  // namespace primerpair::neural_window_golden
