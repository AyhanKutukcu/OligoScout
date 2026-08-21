#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/neural_window_predictor.hpp"

#include <algorithm>
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


struct TimingResult {
    double seconds{0.0};
    std::uint64_t checksum{0};
};


struct NeuralStats {
    std::uint64_t queries{0};
    std::uint64_t fallbacks{0};
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
            reference.at(
                block_start + offset
            ) =
                motif.at(
                    offset %
                    motif.size()
                );
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
     * Header.
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
        const std::size_t first_tab =
            line.find('\t');


        if (
            first_tab ==
            std::string::npos
        ) {
            continue;
        }


        const std::size_t second_tab =
            line.find(
                '\t',
                first_tab + 1
            );


        if (
            second_tab ==
            std::string::npos
        ) {
            continue;
        }


        const std::string_view split(
            line.data() +
                first_tab +
                1,
            second_tab -
                first_tab -
                1
        );


        if (
            split !=
            "test"
        ) {
            continue;
        }


        const std::size_t third_tab =
            line.find(
                '\t',
                second_tab + 1
            );


        if (
            third_tab ==
            std::string::npos
        ) {
            continue;
        }


        queries.emplace_back(
            line.substr(
                second_tab + 1,
                third_tab -
                    second_tab -
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


TimingResult benchmark_binary(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string>& queries,
    const std::size_t repetitions
) {
    std::uint64_t checksum = 0;


    for (
        std::size_t i = 0;
        i <
        std::min<std::size_t>(
            queries.size(),
            5000
        );
        ++i
    ) {
        checksum ^=
            checksum_interval(
                index.exact_search(
                    queries.at(i)
                )
            );
    }


    checksum = 0;


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


TimingResult benchmark_predictor_only(
    const primerpair::NeuralWindowPredictor& predictor,
    const std::vector<std::string>& queries,
    const std::uint64_t row_count,
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
                predictor.predict_row(
                    query,
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


TimingResult benchmark_neural(
    const primerpair::IPBWTIndex& index,
    const primerpair::NeuralWindowPredictor& predictor,
    const std::vector<std::string>& queries,
    const std::uint64_t radius,
    const std::size_t repetitions
) {
    std::uint64_t checksum = 0;


    for (
        std::size_t i = 0;
        i <
        std::min<std::size_t>(
            queries.size(),
            5000
        );
        ++i
    ) {
        const auto prediction =
            predictor.predict_row(
                queries.at(i),
                index.row_count()
            );


        checksum ^=
            checksum_interval(
                index
                .exact_prefix_search_certified_window(
                    queries.at(i),
                    prediction,
                    radius
                )
                .interval
            );
    }


    checksum = 0;


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
            const auto prediction =
                predictor.predict_row(
                    query,
                    index.row_count()
                );


            checksum +=
                checksum_interval(
                    index
                    .exact_prefix_search_certified_window(
                        query,
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


NeuralStats verify_radius(
    const primerpair::IPBWTIndex& index,
    const primerpair::NeuralWindowPredictor& predictor,
    const std::vector<std::string>& queries,
    const std::uint64_t radius
) {
    NeuralStats stats;


    for (
        const auto& query :
        queries
    ) {
        const auto expected =
            index.exact_search(
                query
            );


        const auto prediction =
            predictor.predict_row(
                query,
                index.row_count()
            );


        const auto observed =
            index
            .exact_prefix_search_certified_window(
                query,
                prediction,
                radius
            );


        if (
            expected.begin !=
                observed.interval.begin
            ||
            expected.end !=
                observed.interval.end
        ) {
            std::cerr
                << "CERTIFIED_MISMATCH\t"
                << query
                << '\t'
                << expected.begin
                << '\t'
                << expected.end
                << '\t'
                << observed.interval.begin
                << '\t'
                << observed.interval.end
                << '\n';

            throw std::runtime_error(
                "Certified neural search "
                "differs from binary IP-BWT."
            );
        }


        ++stats.queries;


        if (
            observed.used_global_fallback
        ) {
            ++stats.fallbacks;
        }
    }


    return stats;
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


        const std::string dataset =
            argv[1];


        const std::size_t repetitions =
            static_cast<std::size_t>(
                std::stoull(
                    argv[2]
                )
            );


        const auto queries =
            load_test_queries(
                dataset
            );


        std::cout
            << "test_queries\t"
            << queries.size()
            << '\n';


        if (
            queries.size() !=
            82314
        ) {
            throw std::runtime_error(
                "Unexpected test query count."
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


        NeuralWindowPredictor
            predictor;


        std::cout
            << "ip_bwt_rows\t"
            << index.row_count()
            << '\n';

        std::cout
            << "neural_parameters\t"
            << NeuralWindowPredictor::parameter_count
            << '\n';

        std::cout
            << "neural_weight_bytes\t"
            << NeuralWindowPredictor::fp32_weight_bytes
            << '\n';


        const std::uint64_t radii[] = {
            1024,
            1536,
            2048
        };


        std::vector<NeuralStats>
            verified_stats;


        for (
            const auto radius :
            radii
        ) {
            const auto stats =
                verify_radius(
                    index,
                    predictor,
                    queries,
                    radius
                );


            verified_stats.push_back(
                stats
            );


            const double fallback_rate =
                static_cast<double>(
                    stats.fallbacks
                )
                /
                static_cast<double>(
                    stats.queries
                );


            std::cout
                << "VERIFY"
                << '\t'
                << "radius="
                << radius
                << '\t'
                << "fallbacks="
                << stats.fallbacks
                << '\t'
                << "fallback_rate="
                << std::fixed
                << std::setprecision(8)
                << fallback_rate
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
            << "seconds="
            << binary.seconds
            << '\t'
            << "ns_per_query="
            << binary_ns
            << '\t'
            << "queries_per_sec="
            << (
                total_queries /
                binary.seconds
            )
            << '\t'
            << "checksum="
            << binary.checksum
            << '\n';


        const auto predictor_only =
            benchmark_predictor_only(
                predictor,
                queries,
                index.row_count(),
                repetitions
            );


        std::cout
            << "PREDICTOR_ONLY"
            << '\t'
            << "seconds="
            << predictor_only.seconds
            << '\t'
            << "ns_per_query="
            << (
                predictor_only.seconds *
                1.0e9 /
                total_queries
            )
            << '\t'
            << "checksum="
            << predictor_only.checksum
            << '\n';


        for (
            std::size_t i = 0;
            i < 3;
            ++i
        ) {
            const auto radius =
                radii[i];


            const auto neural =
                benchmark_neural(
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
                    "Timed neural/binary "
                    "checksums differ."
                );
            }


            const double neural_ns =
                neural.seconds *
                1.0e9 /
                total_queries;


            const double speedup =
                binary_ns /
                neural_ns;


            const double fallback_rate =
                static_cast<double>(
                    verified_stats.at(i).fallbacks
                )
                /
                static_cast<double>(
                    verified_stats.at(i).queries
                );


            std::cout
                << "NEURAL"
                << '\t'
                << "radius="
                << radius
                << '\t'
                << "seconds="
                << neural.seconds
                << '\t'
                << "ns_per_query="
                << neural_ns
                << '\t'
                << "queries_per_sec="
                << (
                    total_queries /
                    neural.seconds
                )
                << '\t'
                << "speedup_vs_binary="
                << speedup
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
