#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

#include "primerpair/fm_index.hpp"

int main() {
    const std::string reference =
        "ACGTACGTACGTGATTACAGATTACACCCCC"
        "GGGGGTTTTTAAAAAACGTACGTACGT";

    const std::string pattern =
        "ACGTACGTACGT";

    /*
     * İndeks timer başlamadan önce bir kez oluşturulur.
     * Böylece indeksleme süresi sorgu süresine karışmaz.
     */
    const primerpair::FMIndex index(reference);

    constexpr std::size_t warmup_iterations = 100'000;
    constexpr std::size_t search_iterations = 5'000'000;
    constexpr std::size_t locate_iterations = 500'000;

    std::uint64_t checksum = 0;

    // CPU cache ve branch predictor ısınması.
    for (std::size_t i = 0;
         i < warmup_iterations;
         ++i) {

        const auto interval =
            index.backward_search(pattern);

        checksum += interval.begin;
        checksum += interval.end;
    }

    // --------------------------------------------------------
    // Yalnızca backward-search benchmarkı
    // --------------------------------------------------------

    const auto search_start =
        std::chrono::steady_clock::now();

    for (std::size_t i = 0;
         i < search_iterations;
         ++i) {

        const auto interval =
            index.backward_search(pattern);

        /*
         * Derleyicinin döngüyü kaldırmasını önler.
         */
        checksum += interval.begin;
        checksum += interval.end;
    }

    const auto search_stop =
        std::chrono::steady_clock::now();

    const auto search_total_ns =
        std::chrono::duration_cast<
            std::chrono::nanoseconds
        >(search_stop - search_start).count();

    const double search_ns_per_query =
        static_cast<double>(search_total_ns) /
        static_cast<double>(search_iterations);

    const double search_queries_per_second =
        1'000'000'000.0 / search_ns_per_query;

    // --------------------------------------------------------
    // Backward-search + locate benchmarkı
    // --------------------------------------------------------

    const auto locate_start =
        std::chrono::steady_clock::now();

    for (std::size_t i = 0;
         i < locate_iterations;
         ++i) {

        const auto interval =
            index.backward_search(pattern);

        const auto positions =
            index.locate(interval);

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

    const auto locate_total_ns =
        std::chrono::duration_cast<
            std::chrono::nanoseconds
        >(locate_stop - locate_start).count();

    const double locate_ns_per_query =
        static_cast<double>(locate_total_ns) /
        static_cast<double>(locate_iterations);

    const double locate_queries_per_second =
        1'000'000'000.0 / locate_ns_per_query;

    std::cout << std::fixed
              << std::setprecision(2);

    std::cout
        << "reference_length\t"
        << reference.size()
        << '\n';

    std::cout
        << "pattern\t"
        << pattern
        << '\n';

    std::cout
        << "pattern_length\t"
        << pattern.size()
        << '\n';

    std::cout
        << "search_iterations\t"
        << search_iterations
        << '\n';

    std::cout
        << "search_total_ms\t"
        << static_cast<double>(search_total_ns) /
               1'000'000.0
        << '\n';

    std::cout
        << "search_ns_per_query\t"
        << search_ns_per_query
        << '\n';

    std::cout
        << "search_queries_per_second\t"
        << search_queries_per_second
        << '\n';

    std::cout
        << "search_plus_locate_iterations\t"
        << locate_iterations
        << '\n';

    std::cout
        << "search_plus_locate_total_ms\t"
        << static_cast<double>(locate_total_ns) /
               1'000'000.0
        << '\n';

    std::cout
        << "search_plus_locate_ns_per_query\t"
        << locate_ns_per_query
        << '\n';

    std::cout
        << "search_plus_locate_queries_per_second\t"
        << locate_queries_per_second
        << '\n';

    std::cout
        << "checksum\t"
        << checksum
        << '\n';

    return 0;
}
