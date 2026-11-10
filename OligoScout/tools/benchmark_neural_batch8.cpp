#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/neural_window_predictor.hpp"

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Clock =
    std::chrono::steady_clock;


struct Result {
    double seconds{0.0};
    std::uint64_t checksum{0};
};


struct VerificationStats {
    std::uint64_t queries{0};
    std::array<std::uint64_t, 3>
        fallbacks{};
};


constexpr
std::array<std::uint64_t, 3>
radii{
    1024,
    1536,
    2048
};


std::string make_reference(
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


    /*
     * Must remain identical to the synthetic
     * neural-dataset generator.
     */
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


    /*
     * Skip TSV header.
     */
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
            "Could not open prediction file: "
            +
            path
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



std::uint64_t checksum_interval(
    const primerpair::Interval interval
) {
    return
        interval.begin *
        0x9E3779B185EBCA87ULL
        +
        interval.end;
}


void check_equal(
    const primerpair::Interval expected,
    const primerpair::Interval observed,
    const std::string_view label,
    const std::string_view query
) {
    if (
        expected.begin !=
            observed.begin
        ||
        expected.end !=
            observed.end
    ) {
        std::cerr
            << "INTERVAL_MISMATCH"
            << '\t'
            << label
            << '\t'
            << query
            << '\t'
            << expected.begin
            << '\t'
            << expected.end
            << '\t'
            << observed.begin
            << '\t'
            << observed.end
            << '\n';


        throw std::runtime_error(
            "Exact interval mismatch."
        );
    }
}


VerificationStats verify_all(
    const primerpair::IPBWTIndex& index,
    const primerpair::NeuralWindowPredictor& predictor,
    const std::vector<std::string>& queries
) {
    VerificationStats stats;


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


        const auto binary_batch =
            index.exact_prefix_search_batch8(
                batch
            );


        const auto predictions =
            predictor.predict_rows_batch8(
                batch,
                index.row_count()
            );


        std::array<
            std::array<
                primerpair::IPBWTCertifiedWindowResult,
                8
            >,
            3
        > neural{};


        for (
            std::size_t r = 0;
            r < radii.size();
            ++r
        ) {
            neural[r] =
                index
                .exact_prefix_search_certified_window_batch8(
                    batch,
                    predictions,
                    radii[r]
                );
        }


        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            const auto expected =
                index.exact_search(
                    batch[lane]
                );


            check_equal(
                expected,
                binary_batch[lane],
                "batch_binary",
                batch[lane]
            );


            for (
                std::size_t r = 0;
                r < radii.size();
                ++r
            ) {
                check_equal(
                    expected,
                    neural[r][lane].interval,
                    "batch_neural",
                    batch[lane]
                );


                if (
                    neural[r][lane]
                        .used_global_fallback
                ) {
                    ++stats.fallbacks[r];
                }
            }


            ++stats.queries;
        }
    }


    /*
     * Dataset size is not divisible by 8.
     * Verify the final lanes using the scalar
     * exact/certified paths.
     */
    for (
        ;
        i < queries.size();
        ++i
    ) {
        const auto expected =
            index.exact_search(
                queries[i]
            );


        const auto prediction =
            predictor.predict_row(
                queries[i],
                index.row_count()
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


            check_equal(
                expected,
                observed.interval,
                "scalar_tail_neural",
                queries[i]
            );


            if (
                observed.used_global_fallback
            ) {
                ++stats.fallbacks[r];
            }
        }


        ++stats.queries;
    }


    return stats;
}


void warm_up(
    const primerpair::IPBWTIndex& index,
    const primerpair::NeuralWindowPredictor& predictor,
    const std::vector<std::string>& queries
) {
    const std::size_t limit =
        std::min<std::size_t>(
            queries.size(),
            4096
        );


    std::uint64_t sink = 0;


    for (
        std::size_t i = 0;
        i < limit;
        ++i
    ) {
        sink +=
            checksum_interval(
                index.exact_search(
                    queries[i]
                )
            );
    }


    const std::size_t full =
        (
            limit /
            8
        )
        *
        8;


    for (
        std::size_t i = 0;
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


        const auto binary =
            index.exact_prefix_search_batch8(
                batch
            );


        const auto predictions =
            predictor.predict_rows_batch8(
                batch,
                index.row_count()
            );


        const auto neural =
            index
            .exact_prefix_search_certified_window_batch8(
                batch,
                predictions,
                1536
            );


        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            sink +=
                checksum_interval(
                    binary[lane]
                );

            sink +=
                checksum_interval(
                    neural[lane].interval
                );

            sink +=
                predictions[lane];
        }
    }


    /*
     * Prevent optimizer from discarding warm-up work.
     */
    if (
        sink ==
        0xFFFFFFFFFFFFFFFFULL
    ) {
        std::cerr
            << "warmup_sink\t"
            << sink
            << '\n';
    }
}


Result benchmark_scalar_binary(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string>& queries,
    const std::size_t repetitions
) {
    std::uint64_t checksum = 0;


    const auto start =
        Clock::now();


    for (
        std::size_t repetition = 0;
        repetition < repetitions;
        ++repetition
    ) {
        for (
            const auto& query :
            queries
        ) {
            checksum +=
                checksum_interval(
                    index.exact_search(
                        query
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


Result benchmark_batch_binary(
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


Result benchmark_batch_predictor(
    const primerpair::NeuralWindowPredictor& predictor,
    const std::vector<std::string>& queries,
    const std::uint64_t row_count,
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


            const auto predictions =
                predictor.predict_rows_batch8(
                    batch,
                    row_count
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
                predictor.predict_row(
                    queries[i],
                    row_count
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


Result benchmark_batch_neural(
    const primerpair::IPBWTIndex& index,
    const primerpair::NeuralWindowPredictor& predictor,
    const std::vector<std::string>& queries,
    const std::uint64_t radius,
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


            const auto predictions =
                predictor.predict_rows_batch8(
                    batch,
                    index.row_count()
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
                predictor.predict_row(
                    queries[i],
                    index.row_count()
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



std::vector<std::uint64_t>
precompute_neural_predictions(
    const primerpair::NeuralWindowPredictor& predictor,
    const std::vector<std::string>& queries,
    const std::uint64_t row_count
) {
    std::vector<std::uint64_t>
        predictions(
            queries.size()
        );


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


        const auto batch_predictions =
            predictor.predict_rows_batch8(
                batch,
                row_count
            );


        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            predictions[
                i + lane
            ] =
                batch_predictions[
                    lane
                ];
        }
    }


    for (
        ;
        i < queries.size();
        ++i
    ) {
        predictions[i] =
            predictor.predict_row(
                queries[i],
                row_count
            );
    }


    return predictions;
}


Result benchmark_precomputed_neural_window(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string>& queries,
    const std::vector<std::uint64_t>& predictions,
    const std::uint64_t radius,
    const std::size_t repetitions
) {
    if (
        predictions.size() !=
        queries.size()
    ) {
        throw std::runtime_error(
            "Prediction/query size mismatch."
        );
    }


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


            const std::array<std::uint64_t, 8>
            batch_predictions{
                predictions[i + 0],
                predictions[i + 1],
                predictions[i + 2],
                predictions[i + 3],
                predictions[i + 4],
                predictions[i + 5],
                predictions[i + 6],
                predictions[i + 7]
            };


            const auto intervals =
                index
                .exact_prefix_search_certified_window_batch8(
                    batch,
                    batch_predictions,
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
            checksum +=
                checksum_interval(
                    index
                    .exact_prefix_search_certified_window(
                        queries[i],
                        predictions[i],
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



Result benchmark_oracle_window(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string>& queries,
    const std::vector<std::uint64_t>& oracle_lowers,
    const std::uint64_t radius,
    const std::size_t repetitions
) {
    if (
        oracle_lowers.size() !=
        queries.size()
    ) {
        throw std::runtime_error(
            "Oracle/query size mismatch."
        );
    }


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


            const std::array<std::uint64_t, 8>
            predictions{
                oracle_lowers[i + 0],
                oracle_lowers[i + 1],
                oracle_lowers[i + 2],
                oracle_lowers[i + 3],
                oracle_lowers[i + 4],
                oracle_lowers[i + 5],
                oracle_lowers[i + 6],
                oracle_lowers[i + 7]
            };


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
            checksum +=
                checksum_interval(
                    index
                    .exact_prefix_search_certified_window(
                        queries[i],
                        oracle_lowers[i],
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


int main(
    const int argc,
    char** argv
) {
    using namespace primerpair;

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


        const std::string reference =
            make_reference(
                1000000
            );


        IPBWTIndex index(
            reference,
            21
        );


        NeuralWindowPredictor predictor;


        /*
         * Oracle lower-bounds are generated outside
         * every timed benchmark.
         *
         * This represents a hypothetical perfect,
         * zero-cost predictor.
         */
        std::vector<std::uint64_t>
            oracle_lowers;

        oracle_lowers.reserve(
            queries.size()
        );


        for (
            const auto& query :
            queries
        ) {
            oracle_lowers.push_back(
                index.exact_search(
                    query
                ).begin
            );
        }


        const auto precomputed_predictions =
            precompute_neural_predictions(
                predictor,
                queries,
                index.row_count()
            );


        const auto cdf64_predictions =
            load_predictions(
                "results/benchmarks/"
                "neural_window_v1/"
                "cdf_predictions/"
                "cdf64.txt"
            );


        const auto cdf128_predictions =
            load_predictions(
                "results/benchmarks/"
                "neural_window_v1/"
                "cdf_predictions/"
                "cdf128.txt"
            );


        if (
            cdf64_predictions.size() !=
                queries.size()
            ||
            cdf128_predictions.size() !=
                queries.size()
        ) {
            throw std::runtime_error(
                "CDF prediction count mismatch."
            );
        }


        std::cout
            << "queries\t"
            << queries.size()
            << '\n';

        std::cout
            << "ip_bwt_rows\t"
            << index.row_count()
            << '\n';

        std::cout
            << "neural_parameters\t"
            << NeuralWindowPredictor::parameter_count
            << '\n';


#if defined(PRIMERPAIR_HAVE_AVX2_BACKEND) && PRIMERPAIR_HAVE_AVX2_BACKEND

        std::cout
            << "backend\tAVX2_BATCH8\n";

#else

        std::cout
            << "backend\tSCALAR_FALLBACK\n";

#endif


        /*
         * Full exactness verification happens outside
         * the timed benchmark.
         */
        const auto verification =
            verify_all(
                index,
                predictor,
                queries
            );


        for (
            std::size_t r = 0;
            r < radii.size();
            ++r
        ) {
            const double fallback_rate =
                static_cast<double>(
                    verification.fallbacks[r]
                )
                /
                static_cast<double>(
                    verification.queries
                );


            std::cout
                << "VERIFY"
                << '\t'
                << "radius="
                << radii[r]
                << '\t'
                << "fallbacks="
                << verification.fallbacks[r]
                << '\t'
                << "fallback_rate="
                << std::fixed
                << std::setprecision(8)
                << fallback_rate
                << '\n';
        }


        warm_up(
            index,
            predictor,
            queries
        );


        const double total_queries =
            static_cast<double>(
                queries.size()
            )
            *
            static_cast<double>(
                repetitions
            );


        const auto scalar_binary =
            benchmark_scalar_binary(
                index,
                queries,
                repetitions
            );


        const auto batch_binary =
            benchmark_batch_binary(
                index,
                queries,
                repetitions
            );


        const auto predictor_only =
            benchmark_batch_predictor(
                predictor,
                queries,
                index.row_count(),
                repetitions
            );


        if (
            scalar_binary.checksum !=
            batch_binary.checksum
        ) {
            throw std::runtime_error(
                "Scalar/batch binary checksums differ."
            );
        }


        const double scalar_ns =
            scalar_binary.seconds *
            1.0e9 /
            total_queries;


        const double batch_binary_ns =
            batch_binary.seconds *
            1.0e9 /
            total_queries;


        const double predictor_ns =
            predictor_only.seconds *
            1.0e9 /
            total_queries;


        std::cout
            << "SCALAR_BINARY"
            << '\t'
            << "ns_per_query="
            << scalar_ns
            << '\t'
            << "qps="
            << total_queries /
                scalar_binary.seconds
            << '\t'
            << "checksum="
            << scalar_binary.checksum
            << '\n';


        std::cout
            << "BATCH8_BINARY"
            << '\t'
            << "ns_per_query="
            << batch_binary_ns
            << '\t'
            << "qps="
            << total_queries /
                batch_binary.seconds
            << '\t'
            << "speedup_vs_scalar="
            << scalar_ns /
                batch_binary_ns
            << '\t'
            << "checksum="
            << batch_binary.checksum
            << '\n';


        std::cout
            << "BATCH8_PREDICTOR"
            << '\t'
            << "ns_per_query="
            << predictor_ns
            << '\t'
            << "qps="
            << total_queries /
                predictor_only.seconds
            << '\t'
            << "checksum="
            << predictor_only.checksum
            << '\n';


        /*
         * Perfect-predictor lower bound.
         *
         * radius=0 is the absolute search-side ceiling.
         * Larger radii show how much the required
         * neural error window itself costs.
         */
        constexpr
        std::array<std::uint64_t, 6>
        oracle_radii{
            0,
            32,
            128,
            512,
            1024,
            1536
        };


        for (
            const auto radius :
            oracle_radii
        ) {
            const auto oracle =
                benchmark_oracle_window(
                    index,
                    queries,
                    oracle_lowers,
                    radius,
                    repetitions
                );


            if (
                oracle.checksum !=
                scalar_binary.checksum
            ) {
                throw std::runtime_error(
                    "Oracle/binary checksums differ."
                );
            }


            const double oracle_ns =
                oracle.seconds *
                1.0e9 /
                total_queries;


            std::cout
                << "ORACLE_WINDOW"
                << '\t'
                << "radius="
                << radius
                << '\t'
                << "ns_per_query="
                << oracle_ns
                << '\t'
                << "qps="
                << total_queries /
                    oracle.seconds
                << '\t'
                << "speedup_vs_scalar="
                << scalar_ns /
                    oracle_ns
                << '\t'
                << "speedup_vs_batch_binary="
                << batch_binary_ns /
                    oracle_ns
                << '\t'
                << "predictor_budget_ns="
                << batch_binary_ns -
                    oracle_ns
                << '\t'
                << "checksum="
                << oracle.checksum
                << '\n';
        }


        /*
         * Real model predictions, but prediction
         * computation excluded from timing.
         *
         * This gives the exact inference budget
         * available to beat Batch8 binary.
         */
        for (
            std::size_t r = 0;
            r < radii.size();
            ++r
        ) {
            const auto search_only =
                benchmark_precomputed_neural_window(
                    index,
                    queries,
                    precomputed_predictions,
                    radii[r],
                    repetitions
                );


            if (
                search_only.checksum !=
                scalar_binary.checksum
            ) {
                throw std::runtime_error(
                    "Precomputed neural/binary "
                    "checksums differ."
                );
            }


            const double search_ns =
                search_only.seconds *
                1.0e9 /
                total_queries;


            const double budget =
                batch_binary_ns -
                search_ns;


            std::cout
                << "PRECOMPUTED_NEURAL"
                << '\t'
                << "radius="
                << radii[r]
                << '\t'
                << "search_ns="
                << search_ns
                << '\t'
                << "qps="
                << total_queries /
                    search_only.seconds
                << '\t'
                << "speedup_vs_batch_binary="
                << batch_binary_ns /
                    search_ns
                << '\t'
                << "predictor_budget_ns="
                << budget
                << '\t'
                << "current_predictor_ns="
                << predictor_ns
                << '\t'
                << "required_predictor_speedup="
                << (
                    budget >
                    0.0
                    ?
                    predictor_ns /
                        budget
                    :
                    -1.0
                )
                << '\t'
                << "checksum="
                << search_only.checksum
                << '\n';
        }


        constexpr
        std::array<std::uint64_t, 5>
        cdf_radii{
            2048,
            2560,
            3072,
            3584,
            4096
        };


        for (
            const auto radius :
            cdf_radii
        ) {
            for (
                const auto& item :
                std::array{
                    std::pair{
                        std::string_view(
                            "cdf64"
                        ),
                        &cdf64_predictions
                    },
                    std::pair{
                        std::string_view(
                            "cdf128"
                        ),
                        &cdf128_predictions
                    }
                }
            ) {
                const auto search =
                    benchmark_precomputed_neural_window(
                        index,
                        queries,
                        *item.second,
                        radius,
                        repetitions
                    );


                if (
                    search.checksum !=
                    scalar_binary.checksum
                ) {
                    throw std::runtime_error(
                        "CDF/binary checksum differs."
                    );
                }


                const double search_ns =
                    search.seconds *
                    1.0e9 /
                    total_queries;


                const double budget =
                    batch_binary_ns -
                    search_ns;


                std::cout
                    << "PRECOMPUTED_CDF"
                    << '\t'
                    << "model="
                    << item.first
                    << '\t'
                    << "radius="
                    << radius
                    << '\t'
                    << "search_ns="
                    << search_ns
                    << '\t'
                    << "speedup_vs_batch_binary="
                    << batch_binary_ns /
                        search_ns
                    << '\t'
                    << "predictor_budget_ns="
                    << budget
                    << '\t'
                    << "checksum="
                    << search.checksum
                    << '\n';
            }
        }


        for (
            std::size_t r = 0;
            r < radii.size();
            ++r
        ) {
            const auto neural =
                benchmark_batch_neural(
                    index,
                    predictor,
                    queries,
                    radii[r],
                    repetitions
                );


            if (
                neural.checksum !=
                scalar_binary.checksum
            ) {
                throw std::runtime_error(
                    "Neural/binary checksums differ."
                );
            }


            const double neural_ns =
                neural.seconds *
                1.0e9 /
                total_queries;


            const double fallback_rate =
                static_cast<double>(
                    verification.fallbacks[r]
                )
                /
                static_cast<double>(
                    verification.queries
                );


            std::cout
                << "BATCH8_NEURAL"
                << '\t'
                << "radius="
                << radii[r]
                << '\t'
                << "ns_per_query="
                << neural_ns
                << '\t'
                << "qps="
                << total_queries /
                    neural.seconds
                << '\t'
                << "speedup_vs_scalar="
                << scalar_ns /
                    neural_ns
                << '\t'
                << "speedup_vs_batch_binary="
                << batch_binary_ns /
                    neural_ns
                << '\t'
                << "fallback_rate="
                << fallback_rate
                << '\t'
                << "checksum="
                << neural.checksum
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
