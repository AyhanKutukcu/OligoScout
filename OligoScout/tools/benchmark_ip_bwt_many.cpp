#include "primerpair/ip_bwt_index.hpp"

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


Result
benchmark_scalar(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string_view>& queries,
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
            const auto query :
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


Result
benchmark_direct_batch8(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string_view>& queries,
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


Result
benchmark_many(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string_view>& queries,
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
        const auto intervals =
            index.exact_prefix_search_many(
                queries
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


void
warm_up(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string_view>& queries
) {
    std::uint64_t sink = 0;


    for (
        std::size_t repetition = 0;
        repetition < 2;
        ++repetition
    ) {
        const auto intervals =
            index.exact_prefix_search_many(
                queries
            );


        for (
            const auto interval :
            intervals
        ) {
            sink +=
                checksum_interval(
                    interval
                );
        }
    }


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


}  // namespace


int
main(
    const int argc,
    char** argv
) {
    try {
        if (
            argc != 3
            &&
            argc != 4
        ) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <dataset.tsv> <repetitions> [query_count]\n";

            return 2;
        }


        auto owned =
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
            argc == 4
        ) {
            const std::size_t query_count =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[3]
                    )
                );


            if (
                query_count == 0
                ||
                query_count >
                    owned.size()
            ) {
                throw std::runtime_error(
                    "Invalid query_count."
                );
            }


            owned.resize(
                query_count
            );
        }


        std::vector<std::string_view>
            queries;

        queries.reserve(
            owned.size()
        );


        for (
            const auto& query :
            owned
        ) {
            queries.emplace_back(
                query
            );
        }


        const std::string reference =
            make_reference(
                1000000
            );


        primerpair::IPBWTIndex index(
            reference,
            21
        );


        /*
         * Full correctness check before timing.
         */
        const auto many_check =
            index.exact_prefix_search_many(
                queries
            );


        if (
            many_check.size() !=
            queries.size()
        ) {
            throw std::runtime_error(
                "Many result count mismatch."
            );
        }


        for (
            std::size_t i = 0;
            i < queries.size();
            ++i
        ) {
            const auto expected =
                index.exact_search(
                    queries[i]
                );


            if (
                expected.begin !=
                    many_check[i].begin
                ||
                expected.end !=
                    many_check[i].end
            ) {
                throw std::runtime_error(
                    "Many correctness mismatch."
                );
            }
        }


        warm_up(
            index,
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


        const auto scalar =
            benchmark_scalar(
                index,
                queries,
                repetitions
            );


        const auto batch8 =
            benchmark_direct_batch8(
                index,
                queries,
                repetitions
            );


        const auto many =
            benchmark_many(
                index,
                queries,
                repetitions
            );


        if (
            scalar.checksum !=
                batch8.checksum
            ||
            scalar.checksum !=
                many.checksum
        ) {
            throw std::runtime_error(
                "Benchmark checksum mismatch."
            );
        }


        const double scalar_ns =
            scalar.seconds *
            1.0e9 /
            total_queries;


        const double batch8_ns =
            batch8.seconds *
            1.0e9 /
            total_queries;


        const double many_ns =
            many.seconds *
            1.0e9 /
            total_queries;


        std::cout
            << std::fixed
            << std::setprecision(8);


        std::cout
            << "queries\t"
            << queries.size()
            << '\n';


        std::cout
            << "SCALAR"
            << '\t'
            << "ns_per_query="
            << scalar_ns
            << '\t'
            << "qps="
            << total_queries /
                scalar.seconds
            << '\t'
            << "checksum="
            << scalar.checksum
            << '\n';


        std::cout
            << "DIRECT_BATCH8"
            << '\t'
            << "ns_per_query="
            << batch8_ns
            << '\t'
            << "qps="
            << total_queries /
                batch8.seconds
            << '\t'
            << "speedup_vs_scalar="
            << scalar_ns /
                batch8_ns
            << '\t'
            << "checksum="
            << batch8.checksum
            << '\n';


        std::cout
            << "MANY_API"
            << '\t'
            << "ns_per_query="
            << many_ns
            << '\t'
            << "qps="
            << total_queries /
                many.seconds
            << '\t'
            << "speedup_vs_scalar="
            << scalar_ns /
                many_ns
            << '\t'
            << "relative_to_direct_batch8="
            << batch8_ns /
                many_ns
            << '\t'
            << "overhead_percent="
            << (
                (
                    many_ns /
                    batch8_ns
                )
                -
                1.0
            )
            *
            100.0
            << '\t'
            << "checksum="
            << many.checksum
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
