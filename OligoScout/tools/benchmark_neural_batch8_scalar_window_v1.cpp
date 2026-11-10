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
    std::ifstream input(path);

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
        const auto t1 =
            line.find('\t');

        if (
            t1 ==
            std::string::npos
        ) {
            continue;
        }


        const auto t2 =
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


        const auto t3 =
            line.find(
                '\t',
                t2 + 1
            );


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


std::uint64_t checksum_interval(
    const primerpair::Interval interval
) {
    return
        interval.begin *
        0x9E3779B185EBCA87ULL
        +
        interval.end;
}


Result benchmark_binary(
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


            for (
                std::size_t lane = 0;
                lane < 8;
                ++lane
            ) {
                checksum +=
                    checksum_interval(
                        index
                        .exact_prefix_search_certified_window(
                            batch[lane],
                            predictions[lane],
                            radius
                        )
                        .interval
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

}  // namespace


int main(
    const int argc,
    char** argv
) {
    using namespace primerpair;

    try {
        if (
            argc != 3
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


        const std::string reference =
            make_reference(
                1000000
            );


        IPBWTIndex index(
            reference,
            21
        );


        NeuralWindowPredictor predictor;


        const double total_queries =
            static_cast<double>(
                queries.size()
            )
            *
            static_cast<double>(
                repetitions
            );


        std::cout
            << "queries\t"
            << queries.size()
            << '\n';


#if defined(PRIMERPAIR_HAVE_AVX2_BACKEND) && PRIMERPAIR_HAVE_AVX2_BACKEND

        std::cout
            << "backend\tAVX2_BATCH8\n";

#else

        std::cout
            << "backend\tSCALAR_FALLBACK\n";

#endif


        const auto binary =
            benchmark_binary(
                index,
                queries,
                repetitions
            );


        const double binary_ns =
            binary.seconds *
            1.0e9 /
            total_queries;


        std::cout
            << "BINARY"
            << '\t'
            << "ns_per_query="
            << std::fixed
            << std::setprecision(3)
            << binary_ns
            << '\t'
            << "qps="
            << total_queries /
                binary.seconds
            << '\t'
            << "checksum="
            << binary.checksum
            << '\n';


        const auto predictor_only =
            benchmark_batch_predictor(
                predictor,
                queries,
                index.row_count(),
                repetitions
            );


        const double predictor_ns =
            predictor_only.seconds *
            1.0e9 /
            total_queries;


        std::cout
            << "BATCH8_PREDICTOR"
            << '\t'
            << "ns_per_query="
            << predictor_ns
            << '\t'
            << "qps="
            << total_queries /
                predictor_only.seconds
            << '\n';


        constexpr
        std::uint64_t radius =
            1536;


        const auto neural =
            benchmark_batch_neural(
                index,
                predictor,
                queries,
                radius,
                repetitions
            );


        if (
            neural.checksum !=
            binary.checksum
        ) {
            throw std::runtime_error(
                "Batch8 neural checksum "
                "differs from binary."
            );
        }


        const double neural_ns =
            neural.seconds *
            1.0e9 /
            total_queries;


        std::cout
            << "BATCH8_NEURAL"
            << '\t'
            << "radius="
            << radius
            << '\t'
            << "ns_per_query="
            << neural_ns
            << '\t'
            << "qps="
            << total_queries /
                neural.seconds
            << '\t'
            << "speedup_vs_binary="
            << binary_ns /
                neural_ns
            << '\t'
            << "checksum="
            << neural.checksum
            << '\n';


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
