#include "primerpair/ppfm_io.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double elapsed_seconds(
    const std::chrono::steady_clock::time_point begin,
    const std::chrono::steady_clock::time_point end
) {
    return
        std::chrono::duration<double>(
            end - begin
        ).count();
}


std::uint64_t proc_status_kb(
    const std::string& wanted_key
) {
    std::ifstream input(
        "/proc/self/status"
    );

    if (!input) {
        return 0;
    }

    std::string line;

    while (
        std::getline(
            input,
            line
        )
    ) {
        if (
            line.rfind(
                wanted_key,
                0
            ) != 0
        ) {
            continue;
        }

        std::istringstream parser(
            line.substr(
                wanted_key.size()
            )
        );

        std::uint64_t value = 0;
        std::string unit;

        if (
            parser >>
            value >>
            unit
        ) {
            return value;
        }

        return 0;
    }

    return 0;
}

}  // namespace


int main(
    const int argc,
    char** argv
) {
    using namespace primerpair;

    try {
        if (argc != 4) {
            std::cerr
                << "Usage:\n"
                << "  benchmark_ppfm_load "
                << "<index.ppfm> "
                << "<pattern> "
                << "<expected_position>\n";

            return 2;
        }

        const std::filesystem::path
            path =
                argv[1];

        const std::string pattern =
            argv[2];

        const std::uint64_t
            expected_position =
                std::stoull(
                    argv[3]
                );

        const std::uint64_t
            rss_before_load_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        const auto load_begin =
            std::chrono::
                steady_clock::now();

        auto loaded =
            PpfmIO::load_shard(
                path
            );

        const auto load_end =
            std::chrono::
                steady_clock::now();

        const std::uint64_t
            rss_after_load_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        const std::uint64_t
            hwm_after_load_kb =
                proc_status_kb(
                    "VmHWM:"
                );

        std::cout
            << "chromosome\t"
            << loaded.chromosome
            << '\n';

        std::cout
            << "reference_length\t"
            << loaded.reference.size()
            << '\n';

        std::cout
            << "sa_rate\t"
            << loaded.index
                .forward_index()
                .suffix_array_sample_rate()
            << '\n';

        std::cout
            << "ppfm_bytes\t"
            << std::filesystem::
                file_size(
                    path
                )
            << '\n';

        std::cout
            << "load_seconds\t"
            << elapsed_seconds(
                load_begin,
                load_end
            )
            << '\n';

        std::cout
            << "rss_before_load_kb\t"
            << rss_before_load_kb
            << '\n';

        std::cout
            << "rss_after_load_kb\t"
            << rss_after_load_kb
            << '\n';

        std::cout
            << "hwm_after_load_kb\t"
            << hwm_after_load_kb
            << '\n';

        std::vector<double>
            search_times_us;

        std::vector<std::uint64_t>
            final_positions;

        constexpr std::size_t repeats = 5;

        for (
            std::size_t repeat = 0;
            repeat < repeats;
            ++repeat
        ) {
            const auto search_begin =
                std::chrono::
                    steady_clock::now();

            const auto interval =
                loaded.index
                    .forward_index()
                    .backward_search(
                        pattern
                    );

            final_positions =
                loaded.index
                    .forward_index()
                    .locate(
                        interval
                    );

            const auto search_end =
                std::chrono::
                    steady_clock::now();

            const double microseconds =
                std::chrono::duration<
                    double,
                    std::micro
                >(
                    search_end -
                    search_begin
                ).count();

            search_times_us.push_back(
                microseconds
            );

            std::cout
                << "SEARCH\t"
                << repeat
                << '\t'
                << microseconds
                << '\t'
                << final_positions.size()
                << '\n';
        }

        std::sort(
            search_times_us.begin(),
            search_times_us.end()
        );

        const double median =
            search_times_us.at(
                search_times_us.size() / 2
            );

        const bool expected_found =
            std::find(
                final_positions.begin(),
                final_positions.end(),
                expected_position
            )
            !=
            final_positions.end();

        const std::size_t
            exact_distance =
                loaded.reference
                    .bounded_hamming_distance(
                        expected_position,
                        pattern,
                        0
                    );

        std::cout
            << "search_median_us\t"
            << median
            << '\n';

        std::cout
            << "final_hit_count\t"
            << final_positions.size()
            << '\n';

        for (
            const auto position :
            final_positions
        ) {
            std::cout
                << "HIT\t"
                << position
                << '\n';
        }

        std::cout
            << "expected_position\t"
            << expected_position
            << '\n';

        std::cout
            << "expected_found\t"
            << (
                expected_found
                    ? "YES"
                    : "NO"
            )
            << '\n';

        std::cout
            << "reference_exact_at_expected\t"
            << (
                exact_distance == 0
                    ? "YES"
                    : "NO"
            )
            << '\n';

        const std::uint64_t
            rss_after_search_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        std::cout
            << "rss_after_search_kb\t"
            << rss_after_search_kb
            << '\n';

        if (
            !expected_found
            ||
            exact_distance != 0
        ) {
            throw std::runtime_error(
                "Real PPFM exact-position "
                "verification failed."
            );
        }

        std::cout
            << "ALL_CHECKS\tYES\n";

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
