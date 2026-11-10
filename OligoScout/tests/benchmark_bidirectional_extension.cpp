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

#include "primerpair/bidirectional_fm_index.hpp"

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

        sequence.at(i) =
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
                << " [iterations]"
                << " [sa_sample_rate]\n";

            return 2;
        }

        const std::size_t reference_length =
            static_cast<std::size_t>(
                std::stoull(argv[1])
            );

        const std::size_t iterations =
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
                : primerpair::FMIndex::
                      kDefaultSuffixArraySampleRate;

        constexpr std::size_t seed_length = 12;
        constexpr std::size_t full_length = 13;
        constexpr std::size_t max_pattern_pool = 1024;
        constexpr std::size_t warmup_iterations = 50'000;

        if (reference_length < 100) {
            throw std::invalid_argument(
                "Reference must be at least 100 bp."
            );
        }

        if (iterations == 0) {
            throw std::invalid_argument(
                "Iterations must be greater than zero."
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
            full_length +
            1;

        const std::size_t pattern_pool_size =
            std::min(
                max_pattern_pool,
                available_positions
            );

        std::vector<std::string> full_patterns;
        std::vector<std::string> left_seed_patterns;
        std::vector<std::string> right_seed_patterns;
        std::vector<char> left_bases;
        std::vector<char> right_bases;

        full_patterns.reserve(pattern_pool_size);
        left_seed_patterns.reserve(pattern_pool_size);
        right_seed_patterns.reserve(pattern_pool_size);
        left_bases.reserve(pattern_pool_size);
        right_bases.reserve(pattern_pool_size);

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

            full_patterns.push_back(
                reference.substr(
                    position,
                    full_length
                )
            );

            /*
             * left extension:
             *
             * [base][12-mer seed]
             */
            left_bases.push_back(
                reference.at(position)
            );

            left_seed_patterns.push_back(
                reference.substr(
                    position + 1,
                    seed_length
                )
            );

            /*
             * right extension:
             *
             * [12-mer seed][base]
             */
            right_seed_patterns.push_back(
                reference.substr(
                    position,
                    seed_length
                )
            );

            right_bases.push_back(
                reference.at(
                    position +
                    seed_length
                )
            );
        }

        primerpair::BidirectionalFMIndex index(
            std::move(reference),
            sa_sample_rate
        );

        /*
         * Extension benchmarklarında seed search
         * maliyetini benchmark dışında tutuyoruz.
         */
        std::vector<
            primerpair::BidirectionalInterval
        > left_states;

        std::vector<
            primerpair::BidirectionalInterval
        > right_states;

        left_states.reserve(pattern_pool_size);
        right_states.reserve(pattern_pool_size);

        for (std::size_t i = 0;
             i < pattern_pool_size;
             ++i) {

            left_states.push_back(
                index.search(
                    left_seed_patterns.at(i)
                )
            );

            right_states.push_back(
                index.search(
                    right_seed_patterns.at(i)
                )
            );
        }

        /*
         * Benchmark yalnız başarılı/non-empty
         * extension'lardan oluşmalı.
         */
        for (std::size_t i = 0;
             i < pattern_pool_size;
             ++i) {

            const auto left =
                index.extend_left(
                    left_states.at(i),
                    left_bases.at(i)
                );

            const auto right =
                index.extend_right(
                    right_states.at(i),
                    right_bases.at(i)
                );

            if (
                left.empty() ||
                right.empty()
            ) {
                throw std::logic_error(
                    "Generated benchmark extension "
                    "was unexpectedly empty."
                );
            }
        }

        std::uint64_t checksum = 0;

        /*
         * Warm-up
         */
        for (std::size_t i = 0;
             i < warmup_iterations;
             ++i) {

            const std::size_t slot =
                i % pattern_pool_size;

            const auto fm_interval =
                index.forward_index()
                    .backward_search(
                        full_patterns.at(slot)
                    );

            const auto bifm_state =
                index.search(
                    full_patterns.at(slot)
                );

            const auto left =
                index.extend_left(
                    left_states.at(slot),
                    left_bases.at(slot)
                );

            const auto right =
                index.extend_right(
                    right_states.at(slot),
                    right_bases.at(slot)
                );

            checksum += fm_interval.begin;
            checksum += bifm_state.forward.begin;
            checksum += left.forward.begin;
            checksum += right.forward.begin;
        }

        /*
         * --------------------------------------------------
         * Single FM 13-mer search
         * --------------------------------------------------
         */
        const auto fm_start =
            std::chrono::steady_clock::now();

        for (std::size_t i = 0;
             i < iterations;
             ++i) {

            const std::size_t slot =
                i % pattern_pool_size;

            const auto interval =
                index.forward_index()
                    .backward_search(
                        full_patterns.at(slot)
                    );

            checksum += interval.begin;
            checksum += interval.size();
        }

        const auto fm_stop =
            std::chrono::steady_clock::now();

        /*
         * --------------------------------------------------
         * Full BiFM 13-mer search
         * --------------------------------------------------
         */
        const auto bifm_start =
            std::chrono::steady_clock::now();

        for (std::size_t i = 0;
             i < iterations;
             ++i) {

            const std::size_t slot =
                i % pattern_pool_size;

            const auto state =
                index.search(
                    full_patterns.at(slot)
                );

            checksum += state.forward.begin;
            checksum += state.reverse.begin;
            checksum += state.size();
        }

        const auto bifm_stop =
            std::chrono::steady_clock::now();

        /*
         * --------------------------------------------------
         * Successful left extension
         * --------------------------------------------------
         */
        const auto left_start =
            std::chrono::steady_clock::now();

        for (std::size_t i = 0;
             i < iterations;
             ++i) {

            const std::size_t slot =
                i % pattern_pool_size;

            const auto state =
                index.extend_left(
                    left_states.at(slot),
                    left_bases.at(slot)
                );

            checksum += state.forward.begin;
            checksum += state.reverse.begin;
            checksum += state.size();
        }

        const auto left_stop =
            std::chrono::steady_clock::now();

        /*
         * --------------------------------------------------
         * Successful right extension
         * --------------------------------------------------
         */
        const auto right_start =
            std::chrono::steady_clock::now();

        for (std::size_t i = 0;
             i < iterations;
             ++i) {

            const std::size_t slot =
                i % pattern_pool_size;

            const auto state =
                index.extend_right(
                    right_states.at(slot),
                    right_bases.at(slot)
                );

            checksum += state.forward.begin;
            checksum += state.reverse.begin;
            checksum += state.size();
        }

        const auto right_stop =
            std::chrono::steady_clock::now();

        const auto elapsed_ns =
            [](
                const auto start,
                const auto stop
            ) -> double {

                return static_cast<double>(
                    std::chrono::duration_cast<
                        std::chrono::nanoseconds
                    >(
                        stop - start
                    ).count()
                );
            };

        const double denominator =
            static_cast<double>(
                iterations
            );

        const double fm_per_query =
            elapsed_ns(
                fm_start,
                fm_stop
            ) /
            denominator;

        const double bifm_per_query =
            elapsed_ns(
                bifm_start,
                bifm_stop
            ) /
            denominator;

        const double left_per_query =
            elapsed_ns(
                left_start,
                left_stop
            ) /
            denominator;

        const double right_per_query =
            elapsed_ns(
                right_start,
                right_stop
            ) /
            denominator;

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "reference_length_bp\t"
            << reference_length
            << '\n';

        std::cout
            << "sa_sample_rate\t"
            << sa_sample_rate
            << '\n';

        std::cout
            << "seed_length\t"
            << seed_length
            << '\n';

        std::cout
            << "full_pattern_length\t"
            << full_length
            << '\n';

        std::cout
            << "pattern_pool_size\t"
            << pattern_pool_size
            << '\n';

        std::cout
            << "iterations\t"
            << iterations
            << '\n';

        std::cout
            << "fm_search_13mer_ns_per_query\t"
            << fm_per_query
            << '\n';

        std::cout
            << "bifm_search_13mer_ns_per_query\t"
            << bifm_per_query
            << '\n';

        std::cout
            << "extend_left_ns_per_query\t"
            << left_per_query
            << '\n';

        std::cout
            << "extend_right_ns_per_query\t"
            << right_per_query
            << '\n';

        std::cout
            << "fm_search_queries_per_second\t"
            << 1'000'000'000.0 /
               fm_per_query
            << '\n';

        std::cout
            << "bifm_search_queries_per_second\t"
            << 1'000'000'000.0 /
               bifm_per_query
            << '\n';

        std::cout
            << "extend_left_queries_per_second\t"
            << 1'000'000'000.0 /
               left_per_query
            << '\n';

        std::cout
            << "extend_right_queries_per_second\t"
            << 1'000'000'000.0 /
               right_per_query
            << '\n';

        std::cout
            << "checksum\t"
            << checksum
            << '\n';

        return 0;

    } catch (const std::exception& error) {
        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
