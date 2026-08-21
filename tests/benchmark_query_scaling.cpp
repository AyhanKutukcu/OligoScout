#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "primerpair/fm_index.hpp"

namespace {

std::string generate_reference(
    const std::size_t length
) {
    static constexpr char alphabet[] = {
        'A', 'C', 'G', 'T'
    };

    std::string sequence(length, 'A');

    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;

    for (std::size_t i = 0;
         i < length;
         ++i) {

        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;

        sequence[i] =
            alphabet[state & 3ULL];
    }

    return sequence;
}

}  // namespace

int main(
    int argc,
    char* argv[]
) {
    try {
        if (argc < 2 || argc > 4) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <reference_length_bp>"
                << " [search_iterations]"
                << " [sa_sample_rate]\n";

            return 2;
        }

        const std::size_t reference_length =
            static_cast<std::size_t>(
                std::stoull(argv[1])
            );

        const std::size_t search_iterations =
            argc >= 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 1'000'000;

        const std::size_t sa_sample_rate =
            argc >= 4
                ? static_cast<std::size_t>(
                      std::stoull(argv[3])
                  )
                : primerpair::FMIndex::kDefaultSuffixArraySampleRate;

        constexpr std::size_t pattern_length = 12;
        constexpr std::size_t max_pattern_pool = 1024;
        constexpr std::size_t warmup_iterations = 50'000;

        if (reference_length < 100) {
            throw std::invalid_argument(
                "Reference must be at least 100 bp."
            );
        }

        if (search_iterations == 0) {
            throw std::invalid_argument(
                "Search iterations must be greater than zero."
            );
        }

        if (sa_sample_rate == 0) {
            throw std::invalid_argument(
                "SA sample rate must be greater than zero."
            );
        }

        std::string reference =
            generate_reference(
                reference_length
            );

        const std::size_t available_positions =
            reference_length -
            pattern_length +
            1;

        const std::size_t pattern_pool_size =
            std::min(
                max_pattern_pool,
                available_positions
            );

        std::vector<std::string> patterns;
        patterns.reserve(pattern_pool_size);

        constexpr std::uint64_t stride =
            104729ULL;

        for (std::size_t i = 0;
             i < pattern_pool_size;
             ++i) {

            const std::size_t position =
                static_cast<std::size_t>(
                    (
                        static_cast<std::uint64_t>(i) *
                        stride
                    ) %
                    static_cast<std::uint64_t>(
                        available_positions
                    )
                );

            patterns.push_back(
                reference.substr(
                    position,
                    pattern_length
                )
            );
        }

        /*
         * İndeks timer başlamadan önce oluşturulur.
         * Böylece yalnızca query latency ölçülür.
         */
        primerpair::FMIndex index(
            std::move(reference),
            sa_sample_rate
        );

        std::uint64_t checksum = 0;

        /*
         * Warm-up:
         * cache ve branch predictor etkilerinin
         * daha kararlı hale gelmesi için.
         */
        for (std::size_t i = 0;
             i < warmup_iterations;
             ++i) {

            const std::string& pattern =
                patterns[
                    i % pattern_pool_size
                ];

            const auto interval =
                index.backward_search(
                    pattern
                );

            checksum += interval.begin;
            checksum += interval.end;
        }

        /*
         * ----------------------------------------------------
         * COUNT-ONLY benchmark
         * ----------------------------------------------------
         */
        const auto search_start =
            std::chrono::steady_clock::now();

        for (std::size_t i = 0;
             i < search_iterations;
             ++i) {

            const std::string& pattern =
                patterns[
                    i % pattern_pool_size
                ];

            const auto interval =
                index.backward_search(
                    pattern
                );

            checksum += interval.begin;
            checksum += interval.size();
        }

        const auto search_stop =
            std::chrono::steady_clock::now();

        const auto search_ns =
            std::chrono::duration_cast<
                std::chrono::nanoseconds
            >(
                search_stop -
                search_start
            ).count();

        /*
         * ----------------------------------------------------
         * SEARCH + LOCATE benchmark
         * ----------------------------------------------------
         *
         * Locate daha pahalı olduğu için
         * search_iterations / 10 kadar çalıştırılır.
         */
        const std::size_t locate_iterations =
            std::max<std::size_t>(
                1,
                search_iterations / 10
            );

        std::uint64_t total_hits = 0;

        const auto locate_start =
            std::chrono::steady_clock::now();

        for (std::size_t i = 0;
             i < locate_iterations;
             ++i) {

            const std::string& pattern =
                patterns[
                    i % pattern_pool_size
                ];

            const auto interval =
                index.backward_search(
                    pattern
                );

            const auto positions =
                index.locate(
                    interval
                );

            total_hits +=
                static_cast<std::uint64_t>(
                    positions.size()
                );

            checksum +=
                static_cast<std::uint64_t>(
                    positions.size()
                );

            if (!positions.empty()) {
                checksum += positions.front();
            }
        }

        const auto locate_stop =
            std::chrono::steady_clock::now();

        const auto locate_ns =
            std::chrono::duration_cast<
                std::chrono::nanoseconds
            >(
                locate_stop -
                locate_start
            ).count();

        const double search_ns_per_query =
            static_cast<double>(
                search_ns
            ) /
            static_cast<double>(
                search_iterations
            );

        const double locate_ns_per_query =
            static_cast<double>(
                locate_ns
            ) /
            static_cast<double>(
                locate_iterations
            );

        const double search_qps =
            1'000'000'000.0 /
            search_ns_per_query;

        const double locate_qps =
            1'000'000'000.0 /
            locate_ns_per_query;

        const double mean_hits =
            static_cast<double>(
                total_hits
            ) /
            static_cast<double>(
                locate_iterations
            );

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "reference_length_bp\t"
            << reference_length
            << '\n';

        std::cout
            << "pattern_length\t"
            << pattern_length
            << '\n';

        std::cout
            << "pattern_pool_size\t"
            << pattern_pool_size
            << '\n';

        std::cout
            << "search_iterations\t"
            << search_iterations
            << '\n';

        std::cout
            << "search_ns_per_query\t"
            << search_ns_per_query
            << '\n';

        std::cout
            << "search_queries_per_second\t"
            << search_qps
            << '\n';

        std::cout
            << "search_plus_locate_iterations\t"
            << locate_iterations
            << '\n';

        std::cout
            << "search_plus_locate_ns_per_query\t"
            << locate_ns_per_query
            << '\n';

        std::cout
            << "search_plus_locate_queries_per_second\t"
            << locate_qps
            << '\n';

        std::cout
            << "mean_hits_per_query\t"
            << mean_hits
            << '\n';

        std::cout
            << "checksum\t"
            << checksum
            << '\n';

        return 0;

    } catch (const std::exception& exception) {
        std::cerr
            << "Fatal error: "
            << exception.what()
            << '\n';

        return 1;
    }
}
