#include "primerpair/generated/neural_cdf_v1_weights.hpp"
#include "primerpair/ip_bwt_index.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <immintrin.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifndef __AVX2__
#error "benchmark_cdf_actual requires -mavx2"
#endif

namespace {

using Clock =
    std::chrono::steady_clock;


constexpr
std::array<std::uint64_t, 5>
radii{
    2048,
    2560,
    3072,
    3584,
    4096
};


struct Result {
    double seconds{0.0};
    std::uint64_t checksum{0};
};


struct PredictionVerification {
    std::uint64_t exact{0};
    std::uint64_t within_one{0};
    std::uint64_t over_one{0};
    std::uint64_t max_difference{0};
};


struct SearchVerification {
    std::array<std::uint64_t, 5>
        fallbacks{};
};


std::string
make_reference(
    const std::size_t length
) {
    constexpr char alphabet[] = {
        'A',
        'C',
        'G',
        'T'
    };


    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;


    std::string reference;

    reference.reserve(
        length
    );


    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL
            +
            1442695040888963407ULL;


        reference.push_back(
            alphabet[
                static_cast<std::size_t>(
                    (
                        state >>
                        32
                    )
                    &
                    3ULL
                )
            ]
        );
    }


    const std::string motif =
        "ACGTACGTACGTACGT"
        "AAAAAAAACCCCCCCC"
        "GGGGGGGGTTTTTTTT"
        "ACACACACGTGTGTGT";


    for (
        std::size_t block_start = 8192;
        block_start + 512 < reference.size();
        block_start += 16384
    ) {
        for (
            std::size_t offset = 0;
            offset < 512;
            ++offset
        ) {
            reference[
                block_start +
                offset
            ] =
                motif[
                    offset %
                    motif.size()
                ];
        }
    }


    return reference;
}


std::vector<std::string>
load_test_queries(
    const std::string& path
) {
    std::ifstream input(
        path
    );


    if (!input) {
        throw std::runtime_error(
            "Could not open dataset."
        );
    }


    std::vector<std::string>
        queries;


    std::string line;


    std::getline(
        input,
        line
    );


    while (
        std::getline(
            input,
            line
        )
    ) {
        const std::size_t t1 =
            line.find('\t');


        if (
            t1 ==
            std::string::npos
        ) {
            continue;
        }


        const std::size_t t2 =
            line.find(
                '\t',
                t1 + 1
            );


        if (
            t2 ==
            std::string::npos
        ) {
            continue;
        }


        const std::string_view split(
            line.data() +
                t1 +
                1,
            t2 -
                t1 -
                1
        );


        if (
            split !=
            "test"
        ) {
            continue;
        }


        const std::size_t t3 =
            line.find(
                '\t',
                t2 + 1
            );


        if (
            t3 ==
            std::string::npos
        ) {
            continue;
        }


        queries.emplace_back(
            line.substr(
                t2 + 1,
                t3 -
                    t2 -
                    1
            )
        );
    }


    return queries;
}


std::vector<std::uint64_t>
load_predictions(
    const std::string& path
) {
    std::ifstream input(
        path
    );


    if (!input) {
        throw std::runtime_error(
            "Could not open prediction file."
        );
    }


    std::vector<std::uint64_t>
        predictions;


    std::uint64_t value = 0;


    while (
        input >>
        value
    ) {
        predictions.push_back(
            value
        );
    }


    return predictions;
}


std::uint64_t
checksum_interval(
    const primerpair::Interval interval
) {
    return
        interval.begin *
        0x9E3779B185EBCA87ULL
        +
        interval.end;
}


std::uint64_t
encode_21mer(
    const std::string_view kmer
) {
    if (
        kmer.size() !=
        21
    ) {
        throw std::invalid_argument(
            "CDF predictor requires 21-mer."
        );
    }


    std::uint64_t code = 0;


    for (
        const char base :
        kmer
    ) {
        std::uint64_t value = 0;


        switch (
            base
        ) {
            case 'A':
                value = 0;
                break;

            case 'C':
                value = 1;
                break;

            case 'G':
                value = 2;
                break;

            case 'T':
                value = 3;
                break;

            default:
                throw std::invalid_argument(
                    "CDF predictor accepts "
                    "only A/C/G/T."
                );
        }


        code =
            (
                code << 2U
            )
            |
            value;
    }


    return code;
}


template<std::size_t Hidden>
std::array<std::uint64_t, 8>
predict_batch8(
    const std::array<std::string_view, 8>& kmers,
    const std::uint64_t row_count,
    const std::array<float, Hidden>& hidden_weight,
    const std::array<float, Hidden>& hidden_bias,
    const std::array<float, Hidden>& output_weight,
    const float output_bias,
    const float residual_scale
) {
    std::array<std::uint64_t, 8>
        result{};


    constexpr
    double denominator =
        static_cast<double>(
            std::uint64_t{1}
            <<
            42U
        );


    alignas(32)
    float fractions[8]{};


    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        fractions[lane] =
            static_cast<float>(
                static_cast<double>(
                    encode_21mer(
                        kmers[lane]
                    )
                )
                /
                denominator
            );
    }


    const __m256 x =
        _mm256_load_ps(
            fractions
        );


    const __m256 zero =
        _mm256_setzero_ps();


    __m256 output =
        _mm256_set1_ps(
            output_bias
        );


    /*
     * One-dimensional hidden layer.
     *
     * Eight independent queries occupy
     * the eight AVX2 lanes.
     */
    for (
        std::size_t hidden = 0;
        hidden < Hidden;
        ++hidden
    ) {
        __m256 activation =
            _mm256_add_ps(
                _mm256_mul_ps(
                    x,
                    _mm256_set1_ps(
                        hidden_weight[
                            hidden
                        ]
                    )
                ),
                _mm256_set1_ps(
                    hidden_bias[
                        hidden
                    ]
                )
            );


        activation =
            _mm256_max_ps(
                activation,
                zero
            );


        output =
            _mm256_add_ps(
                output,
                _mm256_mul_ps(
                    activation,
                    _mm256_set1_ps(
                        output_weight[
                            hidden
                        ]
                    )
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


    alignas(32)
    float residuals[8]{};


    _mm256_store_ps(
        residuals,
        output
    );


    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        /*
         * Python model uses a float32 input
         * fraction. Reproduce it exactly here.
         */
        const double baseline =
            1.0
            +
            static_cast<double>(
                fractions[lane]
            )
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


    return result;
}


template<std::size_t Hidden>
std::uint64_t
predict_one(
    const std::string_view kmer,
    const std::uint64_t row_count,
    const std::array<float, Hidden>& hidden_weight,
    const std::array<float, Hidden>& hidden_bias,
    const std::array<float, Hidden>& output_weight,
    const float output_bias,
    const float residual_scale
) {
    constexpr
    double denominator =
        static_cast<double>(
            std::uint64_t{1}
            <<
            42U
        );


    const float fraction =
        static_cast<float>(
            static_cast<double>(
                encode_21mer(
                    kmer
                )
            )
            /
            denominator
        );


    float output =
        output_bias;


    for (
        std::size_t hidden = 0;
        hidden < Hidden;
        ++hidden
    ) {
        float activation =
            fraction *
            hidden_weight[
                hidden
            ]
            +
            hidden_bias[
                hidden
            ];


        activation =
            std::max(
                0.0f,
                activation
            );


        output +=
            activation *
            output_weight[
                hidden
            ];
    }


    const double baseline =
        1.0
        +
        static_cast<double>(
            fraction
        )
        *
        static_cast<double>(
            row_count - 1
        );


    double predicted =
        baseline
        +
        static_cast<double>(
            output *
            residual_scale
        );


    predicted =
        std::clamp(
            predicted,
            0.0,
            static_cast<double>(
                row_count
            )
        );


    return
        static_cast<std::uint64_t>(
            predicted +
            0.5
        );
}


template<std::size_t Hidden>
PredictionVerification
verify_predictions(
    const std::vector<std::string>& queries,
    const std::vector<std::uint64_t>& python_predictions,
    const std::uint64_t row_count,
    const std::array<float, Hidden>& hidden_weight,
    const std::array<float, Hidden>& hidden_bias,
    const std::array<float, Hidden>& output_weight,
    const float output_bias,
    const float residual_scale
) {
    if (
        python_predictions.size() !=
        queries.size()
    ) {
        throw std::runtime_error(
            "Prediction count mismatch."
        );
    }


    PredictionVerification stats;


    const std::size_t full =
        (
            queries.size()
            /
            8
        )
        *
        8;


    std::size_t i = 0;


    for (
        ;
        i < full;
        i += 8
    ) {
        const std::array<std::string_view, 8>
        batch{
            queries[i + 0],
            queries[i + 1],
            queries[i + 2],
            queries[i + 3],
            queries[i + 4],
            queries[i + 5],
            queries[i + 6],
            queries[i + 7]
        };


        const auto cpp =
            predict_batch8(
                batch,
                row_count,
                hidden_weight,
                hidden_bias,
                output_weight,
                output_bias,
                residual_scale
            );


        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            const std::uint64_t expected =
                python_predictions[
                    i + lane
                ];


            const std::uint64_t observed =
                cpp[lane];


            const std::uint64_t difference =
                (
                    expected >
                    observed
                )
                ?
                expected -
                    observed
                :
                observed -
                    expected;


            stats.max_difference =
                std::max(
                    stats.max_difference,
                    difference
                );


            if (
                difference == 0
            ) {
                ++stats.exact;

            } else if (
                difference <= 1
            ) {
                ++stats.within_one;

            } else {

                ++stats.over_one;
            }
        }
    }


    for (
        ;
        i < queries.size();
        ++i
    ) {
        const std::uint64_t observed =
            predict_one(
                queries[i],
                row_count,
                hidden_weight,
                hidden_bias,
                output_weight,
                output_bias,
                residual_scale
            );


        const std::uint64_t expected =
            python_predictions[i];


        const std::uint64_t difference =
            (
                expected >
                observed
            )
            ?
            expected -
                observed
            :
            observed -
                expected;


        stats.max_difference =
            std::max(
                stats.max_difference,
                difference
            );


        if (
            difference == 0
        ) {
            ++stats.exact;

        } else if (
            difference <= 1
        ) {
            ++stats.within_one;

        } else {

            ++stats.over_one;
        }
    }


    /*
     * A larger discrepancy would indicate
     * an export/inference bug, not harmless
     * FP rounding.
     */
    if (
        stats.max_difference >
        2
    ) {
        throw std::runtime_error(
            "CDF C++/Python prediction "
            "difference exceeds 2 rows."
        );
    }


    return stats;
}


template<std::size_t Hidden>
SearchVerification
verify_search(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string>& queries,
    const std::array<float, Hidden>& hidden_weight,
    const std::array<float, Hidden>& hidden_bias,
    const std::array<float, Hidden>& output_weight,
    const float output_bias,
    const float residual_scale
) {
    SearchVerification stats;


    const std::size_t full =
        (
            queries.size()
            /
            8
        )
        *
        8;


    std::size_t i = 0;


    for (
        ;
        i < full;
        i += 8
    ) {
        const std::array<std::string_view, 8>
        batch{
            queries[i + 0],
            queries[i + 1],
            queries[i + 2],
            queries[i + 3],
            queries[i + 4],
            queries[i + 5],
            queries[i + 6],
            queries[i + 7]
        };


        const auto predictions =
            predict_batch8(
                batch,
                index.row_count(),
                hidden_weight,
                hidden_bias,
                output_weight,
                output_bias,
                residual_scale
            );


        for (
            std::size_t r = 0;
            r < radii.size();
            ++r
        ) {
            const auto observed =
                index
                .exact_prefix_search_certified_window_batch8(
                    batch,
                    predictions,
                    radii[r]
                );


            for (
                std::size_t lane = 0;
                lane < 8;
                ++lane
            ) {
                const auto expected =
                    index.exact_search(
                        batch[lane]
                    );


                if (
                    observed[lane]
                        .interval.begin !=
                        expected.begin
                    ||
                    observed[lane]
                        .interval.end !=
                        expected.end
                ) {
                    throw std::runtime_error(
                        "CDF certified interval mismatch."
                    );
                }


                if (
                    observed[lane]
                        .used_global_fallback
                ) {
                    ++stats.fallbacks[r];
                }
            }
        }
    }


    for (
        ;
        i < queries.size();
        ++i
    ) {
        const auto prediction =
            predict_one(
                queries[i],
                index.row_count(),
                hidden_weight,
                hidden_bias,
                output_weight,
                output_bias,
                residual_scale
            );


        const auto expected =
            index.exact_search(
                queries[i]
            );


        for (
            std::size_t r = 0;
            r < radii.size();
            ++r
        ) {
            const auto observed =
                index
                .exact_prefix_search_certified_window(
                    queries[i],
                    prediction,
                    radii[r]
                );


            if (
                observed.interval.begin !=
                    expected.begin
                ||
                observed.interval.end !=
                    expected.end
            ) {
                throw std::runtime_error(
                    "CDF scalar-tail interval mismatch."
                );
            }


            if (
                observed.used_global_fallback
            ) {
                ++stats.fallbacks[r];
            }
        }
    }


    return stats;
}


Result
benchmark_batch_binary(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string>& queries,
    const std::size_t repetitions
) {
    std::uint64_t checksum = 0;


    const std::size_t full =
        (
            queries.size()
            /
            8
        )
        *
        8;


    const auto start =
        Clock::now();


    for (
        std::size_t repetition = 0;
        repetition < repetitions;
        ++repetition
    ) {
        std::size_t i = 0;


        for (
            ;
            i < full;
            i += 8
        ) {
            const std::array<std::string_view, 8>
            batch{
                queries[i + 0],
                queries[i + 1],
                queries[i + 2],
                queries[i + 3],
                queries[i + 4],
                queries[i + 5],
                queries[i + 6],
                queries[i + 7]
            };


            const auto intervals =
                index.exact_prefix_search_batch8(
                    batch
                );


            for (
                const auto interval :
                intervals
            ) {
                checksum +=
                    checksum_interval(
                        interval
                    );
            }
        }


        for (
            ;
            i < queries.size();
            ++i
        ) {
            checksum +=
                checksum_interval(
                    index.exact_search(
                        queries[i]
                    )
                );
        }
    }


    const auto end =
        Clock::now();


    return {
        std::chrono::duration<double>(
            end -
            start
        ).count(),
        checksum
    };
}


template<std::size_t Hidden>
Result
benchmark_predictor(
    const std::vector<std::string>& queries,
    const std::uint64_t row_count,
    const std::size_t repetitions,
    const std::array<float, Hidden>& hidden_weight,
    const std::array<float, Hidden>& hidden_bias,
    const std::array<float, Hidden>& output_weight,
    const float output_bias,
    const float residual_scale
) {
    std::uint64_t checksum = 0;


    const std::size_t full =
        (
            queries.size()
            /
            8
        )
        *
        8;


    const auto start =
        Clock::now();


    for (
        std::size_t repetition = 0;
        repetition < repetitions;
        ++repetition
    ) {
        std::size_t i = 0;


        for (
            ;
            i < full;
            i += 8
        ) {
            const std::array<std::string_view, 8>
            batch{
                queries[i + 0],
                queries[i + 1],
                queries[i + 2],
                queries[i + 3],
                queries[i + 4],
                queries[i + 5],
                queries[i + 6],
                queries[i + 7]
            };


            const auto predictions =
                predict_batch8(
                    batch,
                    row_count,
                    hidden_weight,
                    hidden_bias,
                    output_weight,
                    output_bias,
                    residual_scale
                );


            for (
                const auto prediction :
                predictions
            ) {
                checksum +=
                    prediction;
            }
        }


        for (
            ;
            i < queries.size();
            ++i
        ) {
            checksum +=
                predict_one(
                    queries[i],
                    row_count,
                    hidden_weight,
                    hidden_bias,
                    output_weight,
                    output_bias,
                    residual_scale
                );
        }
    }


    const auto end =
        Clock::now();


    return {
        std::chrono::duration<double>(
            end -
            start
        ).count(),
        checksum
    };
}


template<std::size_t Hidden>
Result
benchmark_actual(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string>& queries,
    const std::uint64_t radius,
    const std::size_t repetitions,
    const std::array<float, Hidden>& hidden_weight,
    const std::array<float, Hidden>& hidden_bias,
    const std::array<float, Hidden>& output_weight,
    const float output_bias,
    const float residual_scale
) {
    std::uint64_t checksum = 0;


    const std::size_t full =
        (
            queries.size()
            /
            8
        )
        *
        8;


    const auto start =
        Clock::now();


    for (
        std::size_t repetition = 0;
        repetition < repetitions;
        ++repetition
    ) {
        std::size_t i = 0;


        for (
            ;
            i < full;
            i += 8
        ) {
            const std::array<std::string_view, 8>
            batch{
                queries[i + 0],
                queries[i + 1],
                queries[i + 2],
                queries[i + 3],
                queries[i + 4],
                queries[i + 5],
                queries[i + 6],
                queries[i + 7]
            };


            const auto predictions =
                predict_batch8(
                    batch,
                    index.row_count(),
                    hidden_weight,
                    hidden_bias,
                    output_weight,
                    output_bias,
                    residual_scale
                );


            const auto intervals =
                index
                .exact_prefix_search_certified_window_batch8(
                    batch,
                    predictions,
                    radius
                );


            for (
                const auto& result :
                intervals
            ) {
                checksum +=
                    checksum_interval(
                        result.interval
                    );
            }
        }


        for (
            ;
            i < queries.size();
            ++i
        ) {
            const auto prediction =
                predict_one(
                    queries[i],
                    index.row_count(),
                    hidden_weight,
                    hidden_bias,
                    output_weight,
                    output_bias,
                    residual_scale
                );


            checksum +=
                checksum_interval(
                    index
                    .exact_prefix_search_certified_window(
                        queries[i],
                        prediction,
                        radius
                    )
                    .interval
                );
        }
    }


    const auto end =
        Clock::now();


    return {
        std::chrono::duration<double>(
            end -
            start
        ).count(),
        checksum
    };
}

}  // namespace


int
main(
    const int argc,
    char** argv
) {
    try {
        if (
            argc !=
            3
        ) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <dataset.tsv> <repetitions>\n";

            return 2;
        }


        const auto queries =
            load_test_queries(
                argv[1]
            );


        const std::size_t repetitions =
            static_cast<std::size_t>(
                std::stoull(
                    argv[2]
                )
            );


        if (
            queries.size() !=
            82314
        ) {
            throw std::runtime_error(
                "Unexpected test-query count."
            );
        }


        const auto python64 =
            load_predictions(
                "results/benchmarks/"
                "neural_window_v1/"
                "cdf_predictions/"
                "cdf64.txt"
            );


        const auto python128 =
            load_predictions(
                "results/benchmarks/"
                "neural_window_v1/"
                "cdf_predictions/"
                "cdf128.txt"
            );


        const std::string reference =
            make_reference(
                1000000
            );


        primerpair::IPBWTIndex index(
            reference,
            21
        );


        namespace cdf64 =
            primerpair::generated::cdf64_v1;

        namespace cdf128 =
            primerpair::generated::cdf128_v1;


        const auto prediction64 =
            verify_predictions(
                queries,
                python64,
                index.row_count(),
                cdf64::hidden_weight,
                cdf64::hidden_bias,
                cdf64::output_weight,
                cdf64::output_bias,
                cdf64::residual_scale
            );


        const auto prediction128 =
            verify_predictions(
                queries,
                python128,
                index.row_count(),
                cdf128::hidden_weight,
                cdf128::hidden_bias,
                cdf128::output_weight,
                cdf128::output_bias,
                cdf128::residual_scale
            );


        const auto search64 =
            verify_search(
                index,
                queries,
                cdf64::hidden_weight,
                cdf64::hidden_bias,
                cdf64::output_weight,
                cdf64::output_bias,
                cdf64::residual_scale
            );


        const auto search128 =
            verify_search(
                index,
                queries,
                cdf128::hidden_weight,
                cdf128::hidden_bias,
                cdf128::output_weight,
                cdf128::output_bias,
                cdf128::residual_scale
            );


        std::cout
            << "queries\t"
            << queries.size()
            << '\n';

        std::cout
            << "ip_bwt_rows\t"
            << index.row_count()
            << '\n';

        std::cout
            << "backend\tAVX2_BATCH8\n";


        std::cout
            << "VERIFY_PREDICTOR"
            << '\t'
            << "model=cdf64"
            << '\t'
            << "exact="
            << prediction64.exact
            << '\t'
            << "within_one="
            << prediction64.within_one
            << '\t'
            << "over_one="
            << prediction64.over_one
            << '\t'
            << "max_diff="
            << prediction64.max_difference
            << '\n';


        std::cout
            << "VERIFY_PREDICTOR"
            << '\t'
            << "model=cdf128"
            << '\t'
            << "exact="
            << prediction128.exact
            << '\t'
            << "within_one="
            << prediction128.within_one
            << '\t'
            << "over_one="
            << prediction128.over_one
            << '\t'
            << "max_diff="
            << prediction128.max_difference
            << '\n';


        for (
            std::size_t r = 0;
            r < radii.size();
            ++r
        ) {
            std::cout
                << "VERIFY_SEARCH"
                << '\t'
                << "model=cdf64"
                << '\t'
                << "radius="
                << radii[r]
                << '\t'
                << "fallbacks="
                << search64.fallbacks[r]
                << '\t'
                << "fallback_rate="
                << std::fixed
                << std::setprecision(8)
                << (
                    static_cast<double>(
                        search64.fallbacks[r]
                    )
                    /
                    static_cast<double>(
                        queries.size()
                    )
                )
                << '\n';


            std::cout
                << "VERIFY_SEARCH"
                << '\t'
                << "model=cdf128"
                << '\t'
                << "radius="
                << radii[r]
                << '\t'
                << "fallbacks="
                << search128.fallbacks[r]
                << '\t'
                << "fallback_rate="
                << (
                    static_cast<double>(
                        search128.fallbacks[r]
                    )
                    /
                    static_cast<double>(
                        queries.size()
                    )
                )
                << '\n';
        }


        const double total_queries =
            static_cast<double>(
                queries.size()
            )
            *
            static_cast<double>(
                repetitions
            );


        const auto binary =
            benchmark_batch_binary(
                index,
                queries,
                repetitions
            );


        const auto predictor64 =
            benchmark_predictor(
                queries,
                index.row_count(),
                repetitions,
                cdf64::hidden_weight,
                cdf64::hidden_bias,
                cdf64::output_weight,
                cdf64::output_bias,
                cdf64::residual_scale
            );


        const auto predictor128 =
            benchmark_predictor(
                queries,
                index.row_count(),
                repetitions,
                cdf128::hidden_weight,
                cdf128::hidden_bias,
                cdf128::output_weight,
                cdf128::output_bias,
                cdf128::residual_scale
            );


        const double binary_ns =
            binary.seconds *
            1.0e9 /
            total_queries;


        const double predictor64_ns =
            predictor64.seconds *
            1.0e9 /
            total_queries;


        const double predictor128_ns =
            predictor128.seconds *
            1.0e9 /
            total_queries;


        std::cout
            << "BATCH8_BINARY"
            << '\t'
            << "ns_per_query="
            << binary_ns
            << '\t'
            << "qps="
            << total_queries /
                binary.seconds
            << '\t'
            << "checksum="
            << binary.checksum
            << '\n';


        std::cout
            << "CDF_PREDICTOR"
            << '\t'
            << "model=cdf64"
            << '\t'
            << "params="
            << cdf64::parameter_count
            << '\t'
            << "ns_per_query="
            << predictor64_ns
            << '\t'
            << "qps="
            << total_queries /
                predictor64.seconds
            << '\t'
            << "checksum="
            << predictor64.checksum
            << '\n';


        std::cout
            << "CDF_PREDICTOR"
            << '\t'
            << "model=cdf128"
            << '\t'
            << "params="
            << cdf128::parameter_count
            << '\t'
            << "ns_per_query="
            << predictor128_ns
            << '\t'
            << "qps="
            << total_queries /
                predictor128.seconds
            << '\t'
            << "checksum="
            << predictor128.checksum
            << '\n';


        for (
            const auto radius :
            radii
        ) {
            const auto actual64 =
                benchmark_actual(
                    index,
                    queries,
                    radius,
                    repetitions,
                    cdf64::hidden_weight,
                    cdf64::hidden_bias,
                    cdf64::output_weight,
                    cdf64::output_bias,
                    cdf64::residual_scale
                );


            const auto actual128 =
                benchmark_actual(
                    index,
                    queries,
                    radius,
                    repetitions,
                    cdf128::hidden_weight,
                    cdf128::hidden_bias,
                    cdf128::output_weight,
                    cdf128::output_bias,
                    cdf128::residual_scale
                );


            if (
                actual64.checksum !=
                    binary.checksum
                ||
                actual128.checksum !=
                    binary.checksum
            ) {
                throw std::runtime_error(
                    "Actual CDF checksum mismatch."
                );
            }


            const double actual64_ns =
                actual64.seconds *
                1.0e9 /
                total_queries;


            const double actual128_ns =
                actual128.seconds *
                1.0e9 /
                total_queries;


            std::cout
                << "ACTUAL_CDF"
                << '\t'
                << "model=cdf64"
                << '\t'
                << "radius="
                << radius
                << '\t'
                << "ns_per_query="
                << actual64_ns
                << '\t'
                << "qps="
                << total_queries /
                    actual64.seconds
                << '\t'
                << "speedup_vs_batch_binary="
                << binary_ns /
                    actual64_ns
                << '\t'
                << "checksum="
                << actual64.checksum
                << '\n';


            std::cout
                << "ACTUAL_CDF"
                << '\t'
                << "model=cdf128"
                << '\t'
                << "radius="
                << radius
                << '\t'
                << "ns_per_query="
                << actual128_ns
                << '\t'
                << "qps="
                << total_queries /
                    actual128.seconds
                << '\t'
                << "speedup_vs_batch_binary="
                << binary_ns /
                    actual128_ns
                << '\t'
                << "checksum="
                << actual128.checksum
                << '\n';
        }


        std::cout
            << "ALL_CHECKS\tYES\n";


        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';

        return 1;
    }
}
