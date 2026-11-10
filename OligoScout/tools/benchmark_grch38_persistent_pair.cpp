#include "primerpair/persistent_genome_search.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/mismatch_features.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

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


bool contains_expected_exact(
    const primerpair::GenomePairSearchResult& result,
    const std::string& chromosome,
    const std::uint64_t expected_start,
    const std::uint64_t expected_length,
    const std::size_t primer2_length
) {
    const std::uint64_t expected_end =
        expected_start +
        expected_length;

    const std::uint64_t expected_right =
        expected_end -
        static_cast<std::uint64_t>(
            primer2_length
        );

    for (
        const auto& shard :
        result.shards
    ) {
        if (
            shard.chromosome !=
            chromosome
        ) {
            continue;
        }

        for (
            const auto& hit :
            shard
                .search_result
                .pair_result
                .amplicons
        ) {
            if (
                hit.left_primer ==
                    primerpair::PrimerIdentity::Primer1
                &&
                hit.right_primer ==
                    primerpair::PrimerIdentity::Primer2
                &&
                hit.left_position ==
                    expected_start
                &&
                hit.right_position ==
                    expected_right
                &&
                hit.amplicon_start ==
                    expected_start
                &&
                hit.amplicon_end_exclusive ==
                    expected_end
                &&
                hit.amplicon_length ==
                    expected_length
                &&
                hit.left_mismatches == 0
                &&
                hit.right_mismatches == 0
                &&
                hit.left_mismatch_mask == 0
                &&
                hit.right_mismatch_mask == 0
            ) {
                return true;
            }
        }
    }

    return false;
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
                << "  benchmark_grch38_persistent_pair "
                << "<manifest> "
                << "<index_dir> "
                << "<primer1> "
                << "<primer2> "
                << "<chr21_start> "
                << "<chr22_start> "
                << "<expected_length>\n";

            return 2;
        }

        const std::string manifest_path =
            argv[1];

        const std::string index_dir =
            argv[2];

        const std::string primer1 =
            argv[3];

        const std::string primer2 =
            argv[4];

        const std::uint64_t chr21_start =
            std::stoull(
                argv[5]
            );

        const std::uint64_t chr22_start =
            std::stoull(
                argv[6]
            );

        const std::uint64_t expected_length =
            std::stoull(
                argv[7]
            );


        const auto manifest =
            PpfmManifest::load(
                manifest_path,
                index_dir
            );

        if (manifest.size() != 24) {
            throw std::runtime_error(
                "Expected exactly 24 canonical "
                "GRCh38 shards."
            );
        }


        PersistentGenomeSearchEngine genome(
            manifest,
            2,
            8
        );


        std::cout
            << "manifest_shards\t"
            << manifest.size()
            << '\n';

        std::cout
            << "cache_capacity\t"
            << genome.cache().capacity()
            << '\n';

        std::cout
            << "rss_before_search_kb\t"
            << proc_status_kb(
                "VmRSS:"
            )
            << '\n';


        const auto begin =
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

        const auto end =
            std::chrono::
                steady_clock::now();


        const double seconds =
            std::chrono::duration<double>(
                end - begin
            ).count();


        std::cout
            << "search_seconds\t"
            << seconds
            << '\n';

        std::cout
            << "shards_searched\t"
            << result.shards.size()
            << '\n';

        std::cout
            << "total_amplicons\t"
            << result.total_amplicon_count()
            << '\n';


        for (
            const auto& shard :
            result.shards
        ) {
            const auto& hits =
                shard
                    .search_result
                    .pair_result
                    .amplicons;

            std::cout
                << "SHARD\t"
                << shard.shard_id
                << '\t'
                << shard.chromosome
                << '\t'
                << hits.size()
                << '\n';

            for (
                std::size_t i = 0;
                i < hits.size();
                ++i
            ) {
                const auto& hit =
                    hits.at(i);

                const std::size_t left_primer_length =
                    (
                        hit.left_primer ==
                        PrimerIdentity::Primer1
                    )
                    ? primer1.size()
                    : primer2.size();

                const std::size_t right_primer_length =
                    (
                        hit.right_primer ==
                        PrimerIdentity::Primer1
                    )
                    ? primer1.size()
                    : primer2.size();

                const auto mismatch_features =
                    extract_pair_mismatch_features(
                        hit,
                        left_primer_length,
                        right_primer_length
                    );

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

                std::cout
                    << "MISMATCH_FEATURES\t"
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
                    << mismatch_features.left.primer_length
                    << '\t'
                    << mismatch_features.left.mismatch_count
                    << '\t'
                    << mismatch_features.left.first_5_count
                    << '\t'
                    << mismatch_features.left.last_8_count
                    << '\t'
                    << mismatch_features.left.last_5_count
                    << '\t'
                    << mismatch_features.left.last_3_count
                    << '\t'
                    << (
                        mismatch_features
                            .left
                            .terminal_3prime_mismatch
                        ? "YES"
                        : "NO"
                    )
                    << '\t'
                    << (
                        mismatch_features
                            .left
                            .nearest_mismatch_to_3prime
                            .has_value()
                        ? std::to_string(
                              *mismatch_features
                                  .left
                                  .nearest_mismatch_to_3prime
                          )
                        : "NA"
                    )
                    << '\t'
                    << to_string(
                        hit.right_primer
                    )
                    << '\t'
                    << mismatch_features.right.primer_length
                    << '\t'
                    << mismatch_features.right.mismatch_count
                    << '\t'
                    << mismatch_features.right.first_5_count
                    << '\t'
                    << mismatch_features.right.last_8_count
                    << '\t'
                    << mismatch_features.right.last_5_count
                    << '\t'
                    << mismatch_features.right.last_3_count
                    << '\t'
                    << (
                        mismatch_features
                            .right
                            .terminal_3prime_mismatch
                        ? "YES"
                        : "NO"
                    )
                    << '\t'
                    << (
                        mismatch_features
                            .right
                            .nearest_mismatch_to_3prime
                            .has_value()
                        ? std::to_string(
                              *mismatch_features
                                  .right
                                  .nearest_mismatch_to_3prime
                          )
                        : "NA"
                    )
                    << '\t'
                    << mismatch_features.total_mismatches
                    << '\n';
            }
        }


        const bool chr21_expected =
            contains_expected_exact(
                result,
                "chr21",
                chr21_start,
                expected_length,
                primer2.size()
            );

        const bool chr22_expected =
            contains_expected_exact(
                result,
                "chr22",
                chr22_start,
                expected_length,
                primer2.size()
            );


        std::cout
            << "chr21_expected_exact\t"
            << (
                chr21_expected
                ? "YES"
                : "NO"
            )
            << '\n';

        std::cout
            << "chr22_expected_exact\t"
            << (
                chr22_expected
                ? "YES"
                : "NO"
            )
            << '\n';


        std::cout
            << "cache_resident_final\t"
            << genome.cache().size()
            << '\n';

        std::cout
            << "cache_loads\t"
            << genome.cache().load_count()
            << '\n';

        std::cout
            << "cache_hits\t"
            << genome.cache().hit_count()
            << '\n';

        std::cout
            << "cache_evictions\t"
            << genome.cache().eviction_count()
            << '\n';

        std::cout
            << "rss_after_search_kb\t"
            << proc_status_kb(
                "VmRSS:"
            )
            << '\n';

        std::cout
            << "hwm_after_search_kb\t"
            << proc_status_kb(
                "VmHWM:"
            )
            << '\n';


        if (
            result.shards.size() != 24
        ) {
            throw std::runtime_error(
                "Not all 24 GRCh38 shards were searched."
            );
        }

        if (
            genome.cache().size() > 2
        ) {
            throw std::runtime_error(
                "Cache exceeded configured capacity."
            );
        }

        if (
            genome.cache().load_count() != 24
        ) {
            throw std::runtime_error(
                "Expected exactly 24 initial shard loads."
            );
        }

        if (
            !chr21_expected
            ||
            !chr22_expected
        ) {
            throw std::runtime_error(
                "Known chr21/chr22 exact regression "
                "product was lost."
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
