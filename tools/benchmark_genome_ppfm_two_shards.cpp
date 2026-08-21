#include "primerpair/genome_search.hpp"
#include "primerpair/primer_pair_search.hpp"

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


void validate_hit(
    const primerpair::GenomeShardPairSearchResult& shard,
    const std::uint64_t expected_start,
    const std::uint64_t expected_length,
    const std::size_t primer2_length
) {
    const auto& hits =
        shard
            .search_result
            .pair_result
            .amplicons;

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

    const std::uint64_t expected_right =
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
            expected_right
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
            ": persistent production hit mismatch."
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
                << "  benchmark_genome_ppfm_two_shards "
                << "<chr21.ppfm> "
                << "<chr21_start> "
                << "<chr22.ppfm> "
                << "<chr22_start> "
                << "<primer1> "
                << "<primer2> "
                << "<amplicon_length>\n";

            return 2;
        }

        const std::filesystem::path
            path1 =
                argv[1];

        const std::uint64_t start1 =
            std::stoull(
                argv[2]
            );

        const std::filesystem::path
            path2 =
                argv[3];

        const std::uint64_t start2 =
            std::stoull(
                argv[4]
            );

        const std::string primer1 =
            argv[5];

        const std::string primer2 =
            argv[6];

        const std::uint64_t
            expected_length =
                std::stoull(
                    argv[7]
                );

        GenomeSearchEngine genome(
            8
        );

        const std::uint64_t
            rss_before_load_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        const auto load1_begin =
            std::chrono::
                steady_clock::now();

        const std::size_t id1 =
            genome.load_ppfm_shard(
                path1
            );

        const auto load1_end =
            std::chrono::
                steady_clock::now();

        const auto load2_begin =
            std::chrono::
                steady_clock::now();

        const std::size_t id2 =
            genome.load_ppfm_shard(
                path2
            );

        const auto load2_end =
            std::chrono::
                steady_clock::now();

        if (
            id1 != 0
            ||
            id2 != 1
        ) {
            throw std::runtime_error(
                "Unexpected persistent shard IDs."
            );
        }

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
            << genome.shard(0).chromosome()
            << '\t'
            << genome.shard(0).sequence_length()
            << '\t'
            << load1_seconds
            << '\n';

        std::cout
            << "SHARD_LOAD\t"
            << genome.shard(1).chromosome()
            << '\t'
            << genome.shard(1).sequence_length()
            << '\t'
            << load2_seconds
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
            << "rss_after_load_kb\t"
            << proc_status_kb(
                "VmRSS:"
            )
            << '\n';

        std::cout
            << "hwm_after_load_kb\t"
            << proc_status_kb(
                "VmHWM:"
            )
            << '\n';

        constexpr std::size_t repeats = 5;

        std::vector<double>
            search_times_us;

        GenomePairSearchResult
            final_result;

        for (
            std::size_t repeat = 0;
            repeat < repeats;
            ++repeat
        ) {
            const auto begin =
                std::chrono::
                    steady_clock::now();

            final_result =
                genome.search_pair(
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

            std::cout
                << "SEARCH\t"
                << repeat
                << '\t'
                << us
                << '\t'
                << final_result
                    .total_amplicon_count()
                << '\n';
        }

        std::sort(
            search_times_us.begin(),
            search_times_us.end()
        );

        if (
            final_result.shards.size() !=
            2
        ) {
            throw std::runtime_error(
                "Expected exactly two searched shards."
            );
        }

        for (
            const auto& shard :
            final_result.shards
        ) {
            std::cout
                << "SEARCH_SHARD\t"
                << shard.shard_id
                << '\t'
                << shard.chromosome
                << '\t'
                << shard
                    .search_result
                    .pair_result
                    .amplicons
                    .size()
                << '\n';

            for (
                std::size_t i = 0;
                i <
                    shard
                        .search_result
                        .pair_result
                        .amplicons
                        .size();
                ++i
            ) {
                const auto& hit =
                    shard
                        .search_result
                        .pair_result
                        .amplicons
                        .at(i);

                std::cout
                    << "AMPLICON\t"
                    << shard.shard_id
                    << '\t'
                    << shard.chromosome
                    << '\t'
                    << i
                    << '\t'
                    << to_string(
                        hit.left_primer
                    )
                    << '\t'
                    << to_string(
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
        }

        validate_hit(
            final_result.shards.at(0),
            start1,
            expected_length,
            primer2.size()
        );

        validate_hit(
            final_result.shards.at(1),
            start2,
            expected_length,
            primer2.size()
        );

        if (
            final_result
                .total_amplicon_count()
            !=
            2
        ) {
            throw std::runtime_error(
                "Expected exactly two total amplicons."
            );
        }

        const double median =
            search_times_us.at(
                search_times_us.size() / 2
            );

        std::cout
            << "search_median_us\t"
            << median
            << '\n';

        std::cout
            << "final_amplicons\t"
            << final_result
                .total_amplicon_count()
            << '\n';

        std::cout
            << "rss_after_search_kb\t"
            << proc_status_kb(
                "VmRSS:"
            )
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
