#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/ip_bwt_rmi.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Clock =
    std::chrono::steady_clock;


struct LookupStats {
    std::uint64_t lookups{0};
    std::uint64_t local{0};
    std::uint64_t fallback{0};
};


struct TimingResult {
    double seconds{0.0};
    std::uint64_t checksum{0};
};


std::string make_reference(
    const std::size_t length
) {
    if (length < 1000) {
        throw std::invalid_argument(
            "Benchmark reference must be at least 1000 bp."
        );
    }

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


    /*
     * Deterministic pseudo-random background.
     */
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

        const std::size_t index =
            static_cast<std::size_t>(
                (
                    state >>
                    32
                )
                &
                3ULL
            );

        reference.push_back(
            alphabet[index]
        );
    }


    /*
     * Add deterministic low-complexity/repetitive
     * islands so the benchmark is not purely random.
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


char mutate_base(
    const char base
) {
    switch (base) {

        case 'A':
            return 'C';

        case 'C':
            return 'G';

        case 'G':
            return 'T';

        case 'T':
            return 'A';

        default:
            throw std::runtime_error(
                "Unexpected nucleotide."
            );
    }
}


std::vector<std::string>
make_queries(
    const std::string& reference,
    const std::size_t query_count
) {
    std::vector<std::string>
        queries;

    queries.reserve(
        query_count
    );


    std::uint64_t state =
        0xD1B54A32D192ED03ULL;


    for (
        std::size_t i = 0;
        i < query_count;
        ++i
    ) {
        const std::size_t length =
            18 +
            (
                i %
                18
            );


        state =
            state *
            2862933555777941757ULL
            +
            3037000493ULL;


        const std::size_t maximum_start =
            reference.size() -
            length;


        const std::size_t start =
            static_cast<std::size_t>(
                state %
                static_cast<std::uint64_t>(
                    maximum_start + 1
                )
            );


        std::string query =
            reference.substr(
                start,
                length
            );


        /*
         * Half present, half one-mismatch queries.
         */
        if (
            (
                i &
                std::size_t{1}
            )
            != 0
        ) {
            const std::size_t mutation_position =
                length /
                2;

            query.at(
                mutation_position
            ) =
                mutate_base(
                    query.at(
                        mutation_position
                    )
                );
        }


        queries.push_back(
            std::move(
                query
            )
        );
    }


    return queries;
}


TimingResult benchmark_binary(
    const primerpair::IPBWTIndex& index,
    const std::vector<std::string>& queries,
    const std::size_t repetitions
) {
    std::uint64_t checksum = 0;


    /*
     * Warmup.
     */
    const std::size_t warmup =
        std::min<std::size_t>(
            queries.size(),
            5000
        );


    for (
        std::size_t i = 0;
        i < warmup;
        ++i
    ) {
        const auto interval =
            index.exact_search(
                queries.at(i)
            );

        checksum ^=
            interval.begin +
            (
                interval.end <<
                1U
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
            const auto interval =
                index.exact_search(
                    query
                );

            checksum +=
                interval.begin *
                0x9E3779B185EBCA87ULL
                +
                interval.end;
        }
    }


    const auto end =
        Clock::now();


    const double seconds =
        std::chrono::duration<double>(
            end -
            start
        ).count();


    return {
        seconds,
        checksum
    };
}


TimingResult benchmark_rmi(
    const primerpair::IPBWTRMI& rmi,
    const std::vector<std::string>& queries,
    const std::size_t repetitions
) {
    std::uint64_t checksum = 0;


    const std::size_t warmup =
        std::min<std::size_t>(
            queries.size(),
            5000
        );


    for (
        std::size_t i = 0;
        i < warmup;
        ++i
    ) {
        const auto interval =
            rmi.exact_search(
                queries.at(i)
            );

        checksum ^=
            interval.begin +
            (
                interval.end <<
                1U
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
            const auto interval =
                rmi.exact_search(
                    query
                );

            checksum +=
                interval.begin *
                0x9E3779B185EBCA87ULL
                +
                interval.end;
        }
    }


    const auto end =
        Clock::now();


    const double seconds =
        std::chrono::duration<double>(
            end -
            start
        ).count();


    return {
        seconds,
        checksum
    };
}


TimingResult benchmark_rmi_galloping(
    const primerpair::IPBWTRMI& rmi,
    const std::vector<std::string>& queries,
    const std::size_t repetitions
) {
    std::uint64_t checksum = 0;


    const std::size_t warmup =
        std::min<std::size_t>(
            queries.size(),
            5000
        );


    for (
        std::size_t i = 0;
        i < warmup;
        ++i
    ) {
        const auto interval =
            rmi.exact_search_galloping(
                queries.at(i)
            );

        checksum ^=
            interval.begin +
            (
                interval.end <<
                1U
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
            const auto interval =
                rmi.exact_search_galloping(
                    query
                );

            checksum +=
                interval.begin *
                0x9E3779B185EBCA87ULL
                +
                interval.end;
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


TimingResult benchmark_rmi_boundary(
    const primerpair::IPBWTRMI& rmi,
    const std::vector<std::string>& queries,
    const std::size_t repetitions
) {
    std::uint64_t checksum = 0;


    const std::size_t warmup =
        std::min<std::size_t>(
            queries.size(),
            5000
        );


    for (
        std::size_t i = 0;
        i < warmup;
        ++i
    ) {
        const auto interval =
            rmi.exact_search_boundary_routed(
                queries.at(i)
            );

        checksum ^=
            interval.begin +
            (
                interval.end <<
                1U
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
            const auto interval =
                rmi.exact_search_boundary_routed(
                    query
                );

            checksum +=
                interval.begin *
                0x9E3779B185EBCA87ULL
                +
                interval.end;
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


/*
 * Mirrors IPBWTIndex/IPBWTRMI exact-search
 * control flow, but records RMI fallback use.
 *
 * This pass is NOT used for latency timing.
 */
primerpair::Interval
rmi_search_with_stats(
    const primerpair::IPBWTIndex& index,
    const primerpair::IPBWTRMI& rmi,
    const std::string_view query,
    LookupStats& stats
) {
    const std::size_t k =
        index.chunk_length();


    std::vector<std::string_view>
        chunks;


    for (
        std::size_t position = 0;
        position < query.size();
        position += k
    ) {
        const std::size_t remaining =
            query.size() -
            position;

        const std::size_t width =
            std::min(
                k,
                remaining
            );

        chunks.push_back(
            query.substr(
                position,
                width
            )
        );
    }


    std::uint64_t low = 0;
    std::uint64_t high =
        index.row_count();


    for (
        auto it = chunks.rbegin();
        it != chunks.rend();
        ++it
    ) {
        const std::string_view chunk =
            *it;


        if (
            chunk.size() <
            k
        ) {
            std::string low_chunk(
                chunk
            );

            low_chunk.push_back(
                '$'
            );

            low_chunk.append(
                k -
                chunk.size() -
                1,
                'A'
            );


            std::string high_chunk(
                chunk
            );

            high_chunk.append(
                k -
                chunk.size(),
                'T'
            );


            const auto low_result =
                rmi.lower_bound_diagnostics(
                    low_chunk,
                    low
                );

            const auto high_result =
                rmi.lower_bound_diagnostics(
                    high_chunk,
                    high
                );


            stats.lookups += 2;

            stats.fallback +=
                static_cast<std::uint64_t>(
                    low_result.used_global_fallback
                )
                +
                static_cast<std::uint64_t>(
                    high_result.used_global_fallback
                );


            stats.local +=
                static_cast<std::uint64_t>(
                    !low_result.used_global_fallback
                )
                +
                static_cast<std::uint64_t>(
                    !high_result.used_global_fallback
                );


            low =
                low_result.position;

            high =
                high_result.position;

        } else {

            const auto low_result =
                rmi.lower_bound_diagnostics(
                    chunk,
                    low
                );

            const auto high_result =
                rmi.lower_bound_diagnostics(
                    chunk,
                    high
                );


            stats.lookups += 2;

            stats.fallback +=
                static_cast<std::uint64_t>(
                    low_result.used_global_fallback
                )
                +
                static_cast<std::uint64_t>(
                    high_result.used_global_fallback
                );


            stats.local +=
                static_cast<std::uint64_t>(
                    !low_result.used_global_fallback
                )
                +
                static_cast<std::uint64_t>(
                    !high_result.used_global_fallback
                );


            low =
                low_result.position;

            high =
                high_result.position;
        }


        if (
            low >=
            high
        ) {
            return {
                low,
                low
            };
        }
    }


    return {
        low,
        high
    };
}


void verify_backend(
    const primerpair::IPBWTIndex& index,
    const primerpair::IPBWTRMI& rmi,
    const std::vector<std::string>& queries
) {
    for (
        std::size_t i = 0;
        i < queries.size();
        ++i
    ) {
        const auto binary =
            index.exact_search(
                queries.at(i)
            );

        const auto learned =
            rmi.exact_search(
                queries.at(i)
            );


        if (
            binary.begin !=
                learned.begin
            ||
            binary.end !=
                learned.end
        ) {
            std::cerr
                << "BACKEND_MISMATCH\t"
                << i
                << '\t'
                << queries.at(i)
                << '\t'
                << binary.begin
                << '\t'
                << binary.end
                << '\t'
                << learned.begin
                << '\t'
                << learned.end
                << '\n';

            throw std::runtime_error(
                "Binary IP-BWT and RMI differ."
            );
        }
    }
}

}  // namespace


int main(
    const int argc,
    char** argv
) {
    using namespace primerpair;

    try {
        std::size_t reference_length =
            1'000'000;

        std::size_t query_count =
            100'000;

        std::size_t repetitions =
            5;


        if (argc >= 2) {
            reference_length =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[1]
                    )
                );
        }

        if (argc >= 3) {
            query_count =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[2]
                    )
                );
        }

        if (argc >= 4) {
            repetitions =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[3]
                    )
                );
        }


        if (
            query_count == 0 ||
            repetitions == 0
        ) {
            throw std::invalid_argument(
                "Query count and repetitions "
                "must be non-zero."
            );
        }


        std::cout
            << std::fixed
            << std::setprecision(3);


        std::cout
            << "reference_length\t"
            << reference_length
            << '\n';

        std::cout
            << "query_count\t"
            << query_count
            << '\n';

        std::cout
            << "repetitions\t"
            << repetitions
            << '\n';


        const std::string reference =
            make_reference(
                reference_length
            );


        const auto queries =
            make_queries(
                reference,
                query_count
            );


        std::cout
            << "building_ip_bwt\tYES\n";


        const auto build_start =
            Clock::now();


        IPBWTIndex index(
            reference,
            21
        );


        const auto build_end =
            Clock::now();


        const double build_seconds =
            std::chrono::duration<double>(
                build_end -
                build_start
            ).count();


        std::cout
            << "ip_bwt_build_seconds\t"
            << build_seconds
            << '\n';

        std::cout
            << "ip_bwt_rows\t"
            << index.row_count()
            << '\n';

        std::cout
            << "ip_bwt_compact_storage_bytes\t"
            << index.compact_storage_bytes()
            << '\n';


        const TimingResult binary =
            benchmark_binary(
                index,
                queries,
                repetitions
            );


        const double total_searches =
            static_cast<double>(
                query_count
            )
            *
            static_cast<double>(
                repetitions
            );


        const double binary_ns =
            (
                binary.seconds *
                1.0e9
            )
            /
            total_searches;


        const double binary_qps =
            total_searches /
            binary.seconds;


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
            << binary_qps
            << '\t'
            << "checksum="
            << binary.checksum
            << '\n';


        const std::size_t leaf_counts[] = {
            64,
            256,
            1024,
            4096
        };


        for (
            const std::size_t leaf_count :
            leaf_counts
        ) {
            const auto rmi_build_start =
                Clock::now();


            IPBWTRMI rmi(
                index,
                leaf_count,
                8
            );


            const auto rmi_build_end =
                Clock::now();


            const double rmi_build_seconds =
                std::chrono::duration<double>(
                    rmi_build_end -
                    rmi_build_start
                ).count();


            /*
             * Full correctness check BEFORE timing.
             */
            verify_backend(
                index,
                rmi,
                queries
            );


            const TimingResult learned =
                benchmark_rmi(
                    rmi,
                    queries,
                    repetitions
                );


            const TimingResult galloping_learned =
                benchmark_rmi_galloping(
                    rmi,
                    queries,
                    repetitions
                );


            const TimingResult boundary_learned =
                benchmark_rmi_boundary(
                    rmi,
                    queries,
                    repetitions
                );


            if (
                galloping_learned.checksum !=
                    binary.checksum
            ) {
                throw std::runtime_error(
                    "Timed galloping-RMI checksum differs."
                );
            }


            if (
                learned.checksum !=
                binary.checksum
            ) {
                throw std::runtime_error(
                    "Timed binary/RMI checksums differ."
                );
            }


            if (
                boundary_learned.checksum !=
                    binary.checksum
            ) {
                throw std::runtime_error(
                    "Timed boundary-RMI checksum differs."
                );
            }


            LookupStats stats;


            for (
                const auto& query :
                queries
            ) {
                const auto interval =
                    rmi_search_with_stats(
                        index,
                        rmi,
                        query,
                        stats
                    );

                const auto expected =
                    index.exact_search(
                        query
                    );

                if (
                    interval.begin !=
                        expected.begin
                    ||
                    interval.end !=
                        expected.end
                ) {
                    throw std::runtime_error(
                        "Diagnostic RMI path differs "
                        "from binary IP-BWT."
                    );
                }
            }


            const double learned_ns =
                (
                    learned.seconds *
                    1.0e9
                )
                /
                total_searches;


            const double learned_qps =
                total_searches /
                learned.seconds;


            const double galloping_ns =
                (
                    galloping_learned.seconds *
                    1.0e9
                )
                /
                total_searches;


            const double galloping_qps =
                total_searches /
                galloping_learned.seconds;


            const double galloping_speedup =
                binary_ns /
                galloping_ns;


            const double boundary_ns =
                (
                    boundary_learned.seconds *
                    1.0e9
                )
                /
                total_searches;


            const double boundary_qps =
                total_searches /
                boundary_learned.seconds;


            const double boundary_speedup =
                binary_ns /
                boundary_ns;


            const double speedup =
                binary_ns /
                learned_ns;


            const double fallback_rate =
                stats.lookups == 0
                ?
                0.0
                :
                static_cast<double>(
                    stats.fallback
                )
                /
                static_cast<double>(
                    stats.lookups
                );


            std::cout
                << "RMI"
                << '\t'
                << "leaves="
                << leaf_count
                << '\t'
                << "build_seconds="
                << rmi_build_seconds
                << '\t'
                << "seconds="
                << learned.seconds
                << '\t'
                << "ns_per_query="
                << learned_ns
                << '\t'
                << "queries_per_sec="
                << learned_qps
                << '\t'
                << "speedup_vs_binary="
                << speedup
                << '\t'
                << "lookups="
                << stats.lookups
                << '\t'
                << "local="
                << stats.local
                << '\t'
                << "fallback="
                << stats.fallback
                << '\t'
                << "fallback_rate="
                << fallback_rate
                << '\t'
                << "mean_leaf_error="
                << static_cast<double>(
                    rmi.mean_leaf_training_error()
                )
                << '\t'
                << "max_leaf_error="
                << rmi.maximum_leaf_training_error()
                << '\t'
                << "checksum="
                << learned.checksum
                << '\n';


            std::cout
                << "RMI_GALLOP"
                << '\t'
                << "leaves="
                << leaf_count
                << '\t'
                << "seconds="
                << galloping_learned.seconds
                << '\t'
                << "ns_per_query="
                << galloping_ns
                << '\t'
                << "queries_per_sec="
                << galloping_qps
                << '\t'
                << "speedup_vs_binary="
                << galloping_speedup
                << '\t'
                << "checksum="
                << galloping_learned.checksum
                << '\n';


            std::cout
                << "RMI_BOUNDARY"
                << '\t'
                << "leaves="
                << leaf_count
                << '\t'
                << "seconds="
                << boundary_learned.seconds
                << '\t'
                << "ns_per_query="
                << boundary_ns
                << '\t'
                << "queries_per_sec="
                << boundary_qps
                << '\t'
                << "speedup_vs_binary="
                << boundary_speedup
                << '\t'
                << "checksum="
                << boundary_learned.checksum
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
