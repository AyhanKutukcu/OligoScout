#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "primerpair/approximate_anchor_search.hpp"

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
        if (argc < 2 || argc > 5) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <reference_length_bp>"
                << " [iterations]"
                << " [primer_length]"
                << " [anchor_length]\n";

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
                : 10'000;

        const std::size_t primer_length =
            argc >= 4
                ? static_cast<std::size_t>(
                      std::stoull(argv[3])
                  )
                : 20;

        const std::size_t anchor_length =
            argc >= 5
                ? static_cast<std::size_t>(
                      std::stoull(argv[4])
                  )
                : 12;

        constexpr std::size_t max_pool = 512;
        constexpr std::size_t warmup = 200;

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

        if (primer_length == 0) {
            throw std::invalid_argument(
                "Primer length must be greater than zero."
            );
        }

        if (
            anchor_length == 0 ||
            anchor_length > primer_length
        ) {
            throw std::invalid_argument(
                "Invalid anchor length."
            );
        }

        std::string reference =
            generate_reference(
                reference_length
            );

        const std::size_t available =
            reference_length -
            primer_length +
            1;

        const std::size_t pool_size =
            std::min(
                max_pool,
                available
            );

        std::vector<std::string> primers;
        primers.reserve(pool_size);

        constexpr std::uint64_t stride =
            104729ULL;

        for (std::size_t i = 0;
             i < pool_size;
             ++i) {

            const std::size_t position =
                static_cast<std::size_t>(
                    (
                        static_cast<std::uint64_t>(i) *
                        stride
                    ) %
                    static_cast<std::uint64_t>(
                        available
                    )
                );

            primers.push_back(
                reference.substr(
                    position,
                    primer_length
                )
            );
        }

        primerpair::BidirectionalFMIndex index(
            std::move(reference)
        );

        const primerpair::ApproximateAnchorSearcher
            searcher(index);

        std::uint64_t checksum = 0;

        /*
         * Kısa warm-up.
         */
        for (std::size_t i = 0;
             i < warmup;
             ++i) {

            const auto result =
                searcher.search_5prime_mismatches(
                    primers.at(
                        i % pool_size
                    ),
                    anchor_length,
                    i % 4
                );

            checksum +=
                result.total_match_count();

            checksum +=
                static_cast<std::uint64_t>(
                    result.hits.size()
                );
        }

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "reference_length_bp\t"
            << reference_length
            << '\n';

        std::cout
            << "primer_length\t"
            << primer_length
            << '\n';

        std::cout
            << "anchor_length\t"
            << anchor_length
            << '\n';

        std::cout
            << "pattern_pool_size\t"
            << pool_size
            << '\n';

        std::cout
            << "iterations_per_budget\t"
            << iterations
            << '\n';

        std::cout
            << "mismatch_budget"
            << '\t'
            << "ns_per_query"
            << '\t'
            << "queries_per_second"
            << '\t'
            << "mean_final_branches"
            << '\t'
            << "mean_match_count"
            << '\n';

        for (std::size_t budget = 0;
             budget <= 3;
             ++budget) {

            std::uint64_t total_branches = 0;
            std::uint64_t total_matches = 0;

            const auto start =
                std::chrono::steady_clock::now();

            for (std::size_t i = 0;
                 i < iterations;
                 ++i) {

                const auto result =
                    searcher.search_5prime_mismatches(
                        primers.at(
                            i % pool_size
                        ),
                        anchor_length,
                        budget
                    );

                total_branches +=
                    static_cast<std::uint64_t>(
                        result.hits.size()
                    );

                total_matches +=
                    result.total_match_count();

                checksum +=
                    result.total_match_count();

                checksum +=
                    static_cast<std::uint64_t>(
                        result.hits.size()
                    );
            }

            const auto stop =
                std::chrono::steady_clock::now();

            const double elapsed_ns =
                static_cast<double>(
                    std::chrono::duration_cast<
                        std::chrono::nanoseconds
                    >(
                        stop - start
                    ).count()
                );

            const double denominator =
                static_cast<double>(
                    iterations
                );

            const double ns_per_query =
                elapsed_ns /
                denominator;

            const double qps =
                1'000'000'000.0 /
                ns_per_query;

            const double mean_branches =
                static_cast<double>(
                    total_branches
                ) /
                denominator;

            const double mean_matches =
                static_cast<double>(
                    total_matches
                ) /
                denominator;

            std::cout
                << budget
                << '\t'
                << ns_per_query
                << '\t'
                << qps
                << '\t'
                << mean_branches
                << '\t'
                << mean_matches
                << '\n';
        }

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
