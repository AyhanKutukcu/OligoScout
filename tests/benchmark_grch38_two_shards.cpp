#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "primerpair/genome_search.hpp"
#include "primerpair/grch38_genome_loader.hpp"

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

    std::string key;

    while (input >> key) {

        if (key == wanted_key) {
            std::uint64_t value = 0;
            std::string unit;

            input
                >> value
                >> unit;

            return value;
        }

        std::string rest;
        std::getline(
            input,
            rest
        );
    }

    return 0;
}


double median(
    std::vector<double> values
) {
    if (values.empty()) {
        return 0.0;
    }

    std::sort(
        values.begin(),
        values.end()
    );

    const std::size_t n =
        values.size();

    if (n % 2 == 1) {
        return values.at(
            n / 2
        );
    }

    return (
        values.at(n / 2 - 1)
        +
        values.at(n / 2)
    ) /
    2.0;
}


}  // namespace


int main(
    int argc,
    char* argv[]
) {
    try {
        if (
            argc < 4
            ||
            argc > 6
        ) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <GRCh38.fa>"
                << " <primer1>"
                << " <primer2>"
                << " [repeats=5]"
                << " [sa_rate=8]\n";

            return 2;
        }


        const std::string fasta =
            argv[1];

        const std::string primer1 =
            argv[2];

        const std::string primer2 =
            argv[3];

        const std::size_t repeats =
            argc >= 5
                ?
                static_cast<std::size_t>(
                    std::stoull(
                        argv[4]
                    )
                )
                :
                5;

        const std::size_t sa_rate =
            argc >= 6
                ?
                static_cast<std::size_t>(
                    std::stoull(
                        argv[5]
                    )
                )
                :
                8;


        if (
            repeats == 0
            ||
            sa_rate == 0
        ) {
            throw std::invalid_argument(
                "repeats and sa_rate must be > 0."
            );
        }


        std::cout
            << std::fixed
            << std::setprecision(
                3
            );


        const std::uint64_t
            rss_before_build_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        const auto build_start =
            std::chrono::
                steady_clock::now();


        primerpair::GenomeSearchEngine
            genome(
                sa_rate
            );


        const auto load_summary =
            primerpair::
                load_grch38_primary_shards(
                    genome,
                    fasta,
                    {
                        "chr21",
                        "chr22"
                    }
                );


        const auto build_stop =
            std::chrono::
                steady_clock::now();


        const std::uint64_t
            rss_after_build_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        const std::uint64_t
            hwm_after_build_kb =
                proc_status_kb(
                    "VmHWM:"
                );


        const double build_seconds =
            static_cast<double>(
                std::chrono::
                    duration_cast<
                        std::chrono::
                            milliseconds
                    >(
                        build_stop -
                        build_start
                    ).count()
            )
            /
            1000.0;


        std::cout
            << "loaded_shards\t"
            << load_summary.records_loaded
            << '\n';

        std::cout
            << "loaded_bases\t"
            << load_summary.bases_loaded
            << '\n';

        std::cout
            << "sa_rate\t"
            << sa_rate
            << '\n';

        std::cout
            << "build_seconds\t"
            << build_seconds
            << '\n';

        std::cout
            << "rss_before_build_kb\t"
            << rss_before_build_kb
            << '\n';

        std::cout
            << "rss_after_build_kb\t"
            << rss_after_build_kb
            << '\n';

        std::cout
            << "hwm_after_build_kb\t"
            << hwm_after_build_kb
            << '\n';


        for (
            std::size_t i = 0;
            i < genome.shard_count();
            ++i
        ) {
            const auto& shard =
                genome.shard(
                    i
                );

            const auto& forward =
                shard
                    .index()
                    .forward_index();

            const auto& reverse =
                shard
                    .index()
                    .reverse_index();


            std::cout
                << "SHARD\t"
                << i
                << '\t'
                << shard.chromosome()
                << '\t'
                << shard.sequence_length()
                << '\t'
                << forward.sampled_sa_memory_bytes()
                << '\t'
                << reverse.sampled_sa_memory_bytes()
                << '\n';
        }


        std::vector<double>
            search_times;

        search_times.reserve(
            repeats
        );


        std::uint64_t
            final_amplicons = 0;


        for (
            std::size_t repeat = 0;
            repeat < repeats;
            ++repeat
        ) {
            const auto start =
                std::chrono::
                    steady_clock::now();


            const auto result =
                genome.search_pair(
                    primer1,
                    primer2,
                    3,
                    50,
                    3000
                );


            const auto stop =
                std::chrono::
                    steady_clock::now();


            final_amplicons =
                result.total_amplicon_count();


            search_times.push_back(
                static_cast<double>(
                    std::chrono::
                        duration_cast<
                            std::chrono::
                                nanoseconds
                        >(
                            stop -
                            start
                        ).count()
                )
                /
                1000.0
            );


            std::cout
                << "SEARCH\t"
                << repeat
                << '\t'
                << search_times.back()
                << '\t'
                << final_amplicons
                << '\n';


            if (
                repeat + 1 ==
                repeats
            ) {
                for (
                    const auto& shard_result :
                    result.shards
                ) {
                    std::cout
                        << "SEARCH_SHARD\t"
                        << shard_result.shard_id
                        << '\t'
                        << shard_result.chromosome
                        << '\t'
                        << shard_result
                            .search_result
                            .pair_result
                            .amplicons
                            .size()
                        << '\n';


                    const auto& amplicons =
                        shard_result
                            .search_result
                            .pair_result
                            .amplicons;

                    for (
                        std::size_t amplicon_index = 0;
                        amplicon_index <
                            amplicons.size();
                        ++amplicon_index
                    ) {
                        const auto& hit =
                            amplicons.at(
                                amplicon_index
                            );

                        std::cout
                            << "AMPLICON\t"
                            << shard_result.shard_id
                            << '\t'
                            << shard_result.chromosome
                            << '\t'
                            << amplicon_index
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
                }
            }
        }


        const std::uint64_t
            rss_after_search_kb =
                proc_status_kb(
                    "VmRSS:"
                );

        std::cout
            << "rss_after_search_kb\t"
            << rss_after_search_kb
            << '\n';


        std::cout
            << "search_median_us\t"
            << median(
                search_times
            )
            << '\n';

        std::cout
            << "final_amplicons\t"
            << final_amplicons
            << '\n';


        if (
            final_amplicons == 0
        ) {
            throw std::runtime_error(
                "Expected real chr22 amplicon "
                "was not recovered."
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
