#include "primerpair/ppfm_io.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/sensitive_pair_constrained_search.hpp"

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


double seconds_between(
    const std::chrono::steady_clock::time_point begin,
    const std::chrono::steady_clock::time_point end
) {
    return
        std::chrono::duration<double>(
            end - begin
        ).count();
}


double microseconds_between(
    const std::chrono::steady_clock::time_point begin,
    const std::chrono::steady_clock::time_point end
) {
    return
        std::chrono::duration<
            double,
            std::micro
        >(
            end - begin
        ).count();
}


void print_hit(
    const primerpair::PpfmShardData& shard,
    const std::size_t index,
    const primerpair::PrimerPairHit& hit
) {
    std::cout
        << "AMPLICON\t"
        << shard.chromosome
        << '\t'
        << index
        << '\t'
        << primerpair::to_string(
            hit.left_primer
        )
        << '\t'
        << primerpair::to_string(
            hit.right_primer
        )
        << '\t'
        << hit.left_position
        << '\t'
        << hit.right_position
        << '\t'
        << hit.amplicon_start
        << '\t'
        << hit.amplicon_end_exclusive
        << '\t'
        << hit.amplicon_length
        << '\t'
        << hit.left_mismatches
        << '\t'
        << hit.right_mismatches
        << '\t'
        << hit.total_mismatches()
        << '\t'
        << hit.left_mismatch_mask
        << '\t'
        << hit.right_mismatch_mask
        << '\n';
}


void validate_expected_hit(
    const primerpair::PpfmShardData& shard,
    const primerpair::SensitivePairConstrainedSearchResult& result,
    const std::uint64_t expected_start,
    const std::uint64_t expected_length,
    const std::size_t primer2_length
) {
    const auto& hits =
        result.pair_result.amplicons;

    if (hits.size() != 1) {
        throw std::runtime_error(
            shard.chromosome +
            ": expected exactly one amplicon."
        );
    }

    const auto& hit =
        hits.front();

    const std::uint64_t expected_end =
        expected_start +
        expected_length;

    const std::uint64_t
        expected_right_position =
            expected_end -
            static_cast<std::uint64_t>(
                primer2_length
            );

    if (
        hit.left_primer !=
            primerpair::PrimerIdentity::Primer1
        ||
        hit.right_primer !=
            primerpair::PrimerIdentity::Primer2
        ||
        hit.left_position !=
            expected_start
        ||
        hit.right_position !=
            expected_right_position
        ||
        hit.amplicon_start !=
            expected_start
        ||
        hit.amplicon_end_exclusive !=
            expected_end
        ||
        hit.amplicon_length !=
            expected_length
        ||
        hit.left_mismatches != 0
        ||
        hit.right_mismatches != 0
        ||
        hit.left_mismatch_mask != 0
        ||
        hit.right_mismatch_mask != 0
    ) {
        throw std::runtime_error(
            shard.chromosome +
            ": persistent pair-search hit does "
            "not match expected exact product."
        );
    }
}

}  // namespace


int main(
    const int argc,
    char** argv
) {
    using namespace primerpair;

    try {
        if (argc != 8) {
            std::cerr
                << "Usage:\n"
                << "  benchmark_ppfm_two_shards_pair "
                << "<ppfm1> "
                << "<expected_start1> "
                << "<ppfm2> "
                << "<expected_start2> "
                << "<primer1> "
                << "<primer2> "
                << "<expected_amplicon_length>\n";

            return 2;
        }

        const std::filesystem::path path1 =
            argv[1];

        const std::uint64_t expected_start1 =
            std::stoull(
                argv[2]
            );

        const std::filesystem::path path2 =
            argv[3];

        const std::uint64_t expected_start2 =
            std::stoull(
                argv[4]
            );

        const std::string primer1 =
            argv[5];

        const std::string primer2 =
            argv[6];

        const std::uint64_t
            expected_amplicon_length =
                std::stoull(
                    argv[7]
                );

        const std::uint64_t
            rss_before_load_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        const auto load1_begin =
            std::chrono::
                steady_clock::now();

        auto shard1 =
            PpfmIO::load_shard(
                path1
            );

        const auto load1_end =
            std::chrono::
                steady_clock::now();

        const std::uint64_t
            rss_after_first_load_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        const auto load2_begin =
            std::chrono::
                steady_clock::now();

        auto shard2 =
            PpfmIO::load_shard(
                path2
            );

        const auto load2_end =
            std::chrono::
                steady_clock::now();

        const std::uint64_t
            rss_after_both_loads_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        const std::uint64_t
            hwm_after_both_loads_kb =
                proc_status_kb(
                    "VmHWM:"
                );

        const double load1_seconds =
            seconds_between(
                load1_begin,
                load1_end
            );

        const double load2_seconds =
            seconds_between(
                load2_begin,
                load2_end
            );

        std::cout
            << "SHARD_LOAD\t"
            << shard1.chromosome
            << '\t'
            << shard1.reference.size()
            << '\t'
            << load1_seconds
            << '\t'
            << std::filesystem::file_size(
                path1
            )
            << '\n';

        std::cout
            << "SHARD_LOAD\t"
            << shard2.chromosome
            << '\t'
            << shard2.reference.size()
            << '\t'
            << load2_seconds
            << '\t'
            << std::filesystem::file_size(
                path2
            )
            << '\n';

        std::cout
            << "combined_load_seconds\t"
            << (
                load1_seconds +
                load2_seconds
            )
            << '\n';

        std::cout
            << "rss_before_load_kb\t"
            << rss_before_load_kb
            << '\n';

        std::cout
            << "rss_after_first_load_kb\t"
            << rss_after_first_load_kb
            << '\n';

        std::cout
            << "rss_after_both_loads_kb\t"
            << rss_after_both_loads_kb
            << '\n';

        std::cout
            << "hwm_after_both_loads_kb\t"
            << hwm_after_both_loads_kb
            << '\n';

        SensitivePairConstrainedSearchEngine
            searcher1(
                shard1.index,
                shard1.reference
            );

        SensitivePairConstrainedSearchEngine
            searcher2(
                shard2.index,
                shard2.reference
            );

        constexpr std::size_t repeats = 5;

        std::vector<double>
            search_times_us;

        SensitivePairConstrainedSearchResult
            final1;

        SensitivePairConstrainedSearchResult
            final2;

        for (
            std::size_t repeat = 0;
            repeat < repeats;
            ++repeat
        ) {
            const auto begin =
                std::chrono::
                    steady_clock::now();

            final1 =
                searcher1.search(
                    primer1,
                    primer2,
                    3,
                    50,
                    3000
                );

            final2 =
                searcher2.search(
                    primer1,
                    primer2,
                    3,
                    50,
                    3000
                );

            const auto end =
                std::chrono::
                    steady_clock::now();

            const double us =
                microseconds_between(
                    begin,
                    end
                );

            search_times_us.push_back(
                us
            );

            const std::size_t total =
                final1
                    .pair_result
                    .amplicons
                    .size()
                +
                final2
                    .pair_result
                    .amplicons
                    .size();

            std::cout
                << "SEARCH\t"
                << repeat
                << '\t'
                << us
                << '\t'
                << total
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

        std::cout
            << "SEARCH_SHARD\t"
            << shard1.chromosome
            << '\t'
            << final1
                .pair_result
                .amplicons
                .size()
            << '\n';

        for (
            std::size_t i = 0;
            i <
                final1
                    .pair_result
                    .amplicons
                    .size();
            ++i
        ) {
            print_hit(
                shard1,
                i,
                final1
                    .pair_result
                    .amplicons
                    .at(i)
            );
        }

        std::cout
            << "SEARCH_SHARD\t"
            << shard2.chromosome
            << '\t'
            << final2
                .pair_result
                .amplicons
                .size()
            << '\n';

        for (
            std::size_t i = 0;
            i <
                final2
                    .pair_result
                    .amplicons
                    .size();
            ++i
        ) {
            print_hit(
                shard2,
                i,
                final2
                    .pair_result
                    .amplicons
                    .at(i)
            );
        }

        validate_expected_hit(
            shard1,
            final1,
            expected_start1,
            expected_amplicon_length,
            primer2.size()
        );

        validate_expected_hit(
            shard2,
            final2,
            expected_start2,
            expected_amplicon_length,
            primer2.size()
        );

        const std::size_t
            final_amplicons =
                final1
                    .pair_result
                    .amplicons
                    .size()
                +
                final2
                    .pair_result
                    .amplicons
                    .size();

        if (final_amplicons != 2) {
            throw std::runtime_error(
                "Expected exactly two total "
                "persistent amplicons."
            );
        }

        const std::uint64_t
            rss_after_search_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        std::cout
            << "search_median_us\t"
            << median
            << '\n';

        std::cout
            << "final_amplicons\t"
            << final_amplicons
            << '\n';

        std::cout
            << "rss_after_search_kb\t"
            << rss_after_search_kb
            << '\n';

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
