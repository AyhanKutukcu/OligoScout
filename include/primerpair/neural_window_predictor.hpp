#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace primerpair {

class NeuralWindowPredictor {
public:

    static constexpr
    std::size_t kmer_length = 21;

    static constexpr
    std::size_t parameter_count = 2753;

    static constexpr
    std::size_t fp32_weight_bytes =
        parameter_count *
        sizeof(float);


    [[nodiscard]]
    float predict_residual_rows(
        std::string_view kmer
    ) const;


    [[nodiscard]]
    std::uint64_t predict_row(
        std::string_view kmer,
        std::uint64_t row_count
    ) const;


    /*
     * Throughput-oriented inference.
     *
     * AVX2 builds evaluate eight independent
     * 21-mers in parallel. Non-AVX2 builds
     * transparently fall back to scalar inference.
     */
    [[nodiscard]]
    std::array<std::uint64_t, 8>
    predict_rows_batch8(
        const std::array<std::string_view, 8>& kmers,
        std::uint64_t row_count
    ) const;


private:

    struct PreparedInput {
        std::array<float, 9>
            features{};

        std::uint64_t
            full_code{0};
    };


    [[nodiscard]]
    static
    PreparedInput prepare_input(
        std::string_view kmer
    );


    [[nodiscard]]
    float predict_residual_rows_from_features(
        const std::array<float, 9>& features
    ) const;
};

}  // namespace primerpair
