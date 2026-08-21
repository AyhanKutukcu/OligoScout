#include "primerpair/neural_window_predictor.hpp"

#include "primerpair/generated/neural_window_tiny64x32_v3_weights.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string_view>

#if defined(PRIMERPAIR_HAVE_AVX2_BACKEND) && PRIMERPAIR_HAVE_AVX2_BACKEND
#include <immintrin.h>
#endif

namespace primerpair {

namespace {

[[nodiscard]]
std::uint64_t base_value(
    const char base
) {
    switch (base) {

        case 'A':
            return 0;

        case 'C':
            return 1;

        case 'G':
            return 2;

        case 'T':
            return 3;

        default:
            throw std::invalid_argument(
                "Neural window predictor requires "
                "canonical A/C/G/T."
            );
    }
}


/*
 * For a nucleotide count c in a 21-mer:
 *
 *     contribution =
 *     -(c/21) * log2(c/21)
 *
 * Precomputed once offline.
 *
 * This removes std::log2() from every query.
 */
inline constexpr
std::array<double, 22>
entropy_contribution{
    0.0,
    0.20915797251327431,
    0.32307784978845333,
    0.40105070315108626,
    0.45567950910071620,
    0.49294983997414238,
    0.51638712058788683,
    0.52832083357371873,
    0.53040663724905157,
    0.52388246628704915,
    0.50970920375780859,
    0.48865446883600444,
    0.46134566974720237,
    0.42830524572808026,
    0.38997500048077083,
    0.34673344797874411,
    0.29890851259334122,
    0.24678704218967404,
    0.19062207543124116,
    0.13063848939849149,
    0.067037455134664789,
    0.0
};


[[nodiscard]]
double prefix_fraction(
    const std::uint64_t code,
    const unsigned length
) {
    const std::uint64_t denominator =
        std::uint64_t{1}
        <<
        (
            2U *
            length
        );


    return
        static_cast<double>(
            code
        )
        /
        static_cast<double>(
            denominator
        );
}

}  // namespace


NeuralWindowPredictor::PreparedInput
NeuralWindowPredictor::prepare_input(
    const std::string_view kmer
) {
    if (
        kmer.size() !=
        kmer_length
    ) {
        throw std::invalid_argument(
            "Neural window predictor requires "
            "exactly 21 bases."
        );
    }


    PreparedInput result;


    std::array<std::size_t, 4>
        counts{};


    std::uint64_t code = 0;


    std::size_t current_run = 1;
    std::size_t maximum_run = 1;


    for (
        std::size_t i = 0;
        i < kmer_length;
        ++i
    ) {
        const std::uint64_t value =
            base_value(
                kmer[i]
            );


        code =
            (
                code <<
                2U
            )
            |
            value;


        ++counts[
            static_cast<std::size_t>(
                value
            )
        ];


        if (
            i != 0
            &&
            kmer[i] ==
            kmer[
                i - 1
            ]
        ) {
            ++current_run;

            maximum_run =
                std::max(
                    maximum_run,
                    current_run
                );

        } else {

            current_run = 1;
        }


        const std::size_t length =
            i + 1;


        switch (
            length
        ) {

            case 4:
                result.features[0] =
                    static_cast<float>(
                        prefix_fraction(
                            code,
                            4
                        )
                    );
                break;


            case 8:
                result.features[1] =
                    static_cast<float>(
                        prefix_fraction(
                            code,
                            8
                        )
                    );
                break;


            case 12:
                result.features[2] =
                    static_cast<float>(
                        prefix_fraction(
                            code,
                            12
                        )
                    );
                break;


            case 16:
                result.features[3] =
                    static_cast<float>(
                        prefix_fraction(
                            code,
                            16
                        )
                    );
                break;


            case 21:
                result.features[4] =
                    static_cast<float>(
                        prefix_fraction(
                            code,
                            21
                        )
                    );
                break;


            default:
                break;
        }
    }


    result.full_code =
        code;


    const std::size_t gc =
        counts[1]
        +
        counts[2];


    result.features[5] =
        static_cast<float>(
            static_cast<double>(
                gc
            )
            /
            21.0
        );


    double entropy = 0.0;

    std::size_t distinct = 0;


    for (
        const std::size_t count :
        counts
    ) {
        entropy +=
            entropy_contribution[
                count
            ];


        if (
            count != 0
        ) {
            ++distinct;
        }
    }


    result.features[6] =
        static_cast<float>(
            entropy /
            2.0
        );


    result.features[7] =
        static_cast<float>(
            static_cast<double>(
                maximum_run
            )
            /
            21.0
        );


    result.features[8] =
        static_cast<float>(
            static_cast<double>(
                distinct
            )
            /
            4.0
        );


    return result;
}


float
NeuralWindowPredictor::
predict_residual_rows_from_features(
    const std::array<float, 9>& features
) const {
    using namespace
        neural_window_weights;


    std::array<float, 64>
        hidden1{};


    for (
        std::size_t output = 0;
        output < 64;
        ++output
    ) {
        const float* weight =
            layer1_weight.data()
            +
            output *
            9;


        float sum =
            layer1_bias[
                output
            ];


        for (
            std::size_t input = 0;
            input < 9;
            ++input
        ) {
            sum +=
                weight[input]
                *
                features[input];
        }


        hidden1[output] =
            std::max(
                0.0f,
                sum
            );
    }


    std::array<float, 32>
        hidden2{};


    for (
        std::size_t output = 0;
        output < 32;
        ++output
    ) {
        const float* weight =
            layer2_weight.data()
            +
            output *
            64;


        float sum =
            layer2_bias[
                output
            ];


        for (
            std::size_t input = 0;
            input < 64;
            ++input
        ) {
            sum +=
                weight[input]
                *
                hidden1[input];
        }


        hidden2[output] =
            std::max(
                0.0f,
                sum
            );
    }


    float output =
        output_bias[0];


    for (
        std::size_t input = 0;
        input < 32;
        ++input
    ) {
        output +=
            output_weight[
                input
            ]
            *
            hidden2[
                input
            ];
    }


    return
        output
        *
        residual_scale;
}


float
NeuralWindowPredictor::predict_residual_rows(
    const std::string_view kmer
) const {
    const PreparedInput prepared =
        prepare_input(
            kmer
        );


    return
        predict_residual_rows_from_features(
            prepared.features
        );
}


std::uint64_t
NeuralWindowPredictor::predict_row(
    const std::string_view kmer,
    const std::uint64_t row_count
) const {
    if (
        row_count == 0
    ) {
        return 0;
    }


    const PreparedInput prepared =
        prepare_input(
            kmer
        );


    constexpr
    double denominator =
        static_cast<double>(
            std::uint64_t{1}
            <<
            42U
        );


    const double fraction =
        static_cast<double>(
            prepared.full_code
        )
        /
        denominator;


    const double baseline =
        1.0
        +
        fraction
        *
        static_cast<double>(
            row_count - 1
        );


    double predicted =
        baseline
        +
        static_cast<double>(
            predict_residual_rows_from_features(
                prepared.features
            )
        );


    predicted =
        std::clamp(
            predicted,
            0.0,
            static_cast<double>(
                row_count
            )
        );


    /*
     * Equivalent to llround() for the clamped,
     * non-negative prediction domain.
     */
    return
        static_cast<std::uint64_t>(
            predicted +
            0.5
        );
}


std::array<std::uint64_t, 8>
NeuralWindowPredictor::predict_rows_batch8(
    const std::array<std::string_view, 8>& kmers,
    const std::uint64_t row_count
) const {
    std::array<std::uint64_t, 8>
        result{};


    if (
        row_count == 0
    ) {
        return result;
    }


#if defined(PRIMERPAIR_HAVE_AVX2_BACKEND) && PRIMERPAIR_HAVE_AVX2_BACKEND

    using namespace
        neural_window_weights;


    /*
     * Feature extraction is still scalar because
     * each DNA string is only 21 bases.
     *
     * The expensive dense network is evaluated
     * across eight independent queries in AVX2 lanes.
     */
    std::array<PreparedInput, 8>
        prepared{};


    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        prepared[lane] =
            prepare_input(
                kmers[lane]
            );
    }


    /*
     * Transpose:
     *
     * feature_by_input[input][lane]
     *
     * so one AVX2 load contains the same feature
     * from eight independent queries.
     */
    std::array<
        std::array<float, 8>,
        9
    > feature_by_input{};


    for (
        std::size_t input = 0;
        input < 9;
        ++input
    ) {
        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            feature_by_input[input][lane] =
                prepared[lane]
                    .features[input];
        }
    }


    alignas(32)
    __m256 hidden1[64]{};


    for (
        std::size_t output = 0;
        output < 64;
        ++output
    ) {
        __m256 sum =
            _mm256_set1_ps(
                layer1_bias[output]
            );


        const float* weight =
            layer1_weight.data()
            +
            output *
            9;


        for (
            std::size_t input = 0;
            input < 9;
            ++input
        ) {
            const __m256 feature =
                _mm256_loadu_ps(
                    feature_by_input[
                        input
                    ].data()
                );


            const __m256 w =
                _mm256_set1_ps(
                    weight[input]
                );


            sum =
                _mm256_add_ps(
                    sum,
                    _mm256_mul_ps(
                        feature,
                        w
                    )
                );
        }


        hidden1[output] =
            _mm256_max_ps(
                sum,
                _mm256_setzero_ps()
            );
    }


    alignas(32)
    __m256 hidden2[32]{};


    for (
        std::size_t output = 0;
        output < 32;
        ++output
    ) {
        __m256 sum =
            _mm256_set1_ps(
                layer2_bias[output]
            );


        const float* weight =
            layer2_weight.data()
            +
            output *
            64;


        for (
            std::size_t input = 0;
            input < 64;
            ++input
        ) {
            const __m256 w =
                _mm256_set1_ps(
                    weight[input]
                );


            sum =
                _mm256_add_ps(
                    sum,
                    _mm256_mul_ps(
                        hidden1[input],
                        w
                    )
                );
        }


        hidden2[output] =
            _mm256_max_ps(
                sum,
                _mm256_setzero_ps()
            );
    }


    __m256 output =
        _mm256_set1_ps(
            output_bias[0]
        );


    for (
        std::size_t input = 0;
        input < 32;
        ++input
    ) {
        const __m256 w =
            _mm256_set1_ps(
                output_weight[input]
            );


        output =
            _mm256_add_ps(
                output,
                _mm256_mul_ps(
                    hidden2[input],
                    w
                )
            );
    }


    output =
        _mm256_mul_ps(
            output,
            _mm256_set1_ps(
                residual_scale
            )
        );


    std::array<float, 8>
        residuals{};


    _mm256_storeu_ps(
        residuals.data(),
        output
    );


    constexpr
    double denominator =
        static_cast<double>(
            std::uint64_t{1}
            <<
            42U
        );


    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        const double fraction =
            static_cast<double>(
                prepared[lane]
                    .full_code
            )
            /
            denominator;


        const double baseline =
            1.0
            +
            fraction
            *
            static_cast<double>(
                row_count - 1
            );


        double predicted =
            baseline
            +
            static_cast<double>(
                residuals[lane]
            );


        predicted =
            std::clamp(
                predicted,
                0.0,
                static_cast<double>(
                    row_count
                )
            );


        result[lane] =
            static_cast<std::uint64_t>(
                predicted +
                0.5
            );
    }


#else

    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        result[lane] =
            predict_row(
                kmers[lane],
                row_count
            );
    }

#endif


    return result;
}


}  // namespace primerpair
