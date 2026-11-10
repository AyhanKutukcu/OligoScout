#include "primerpair/persistent_multiplex_primer_search_v2.hpp"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct ExpectedShard {
    std::size_t total_hits{0};
    std::size_t window_candidates{0};
    std::size_t cross_products{0};
};


struct OwnedPrimerPair {
    std::string primer1;
    std::string primer2;
};


std::vector<OwnedPrimerPair>
load_panel(
    const std::string& path
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open panel: " + path
        );
    }


    std::string line;

    if (!std::getline(input, line)) {
        throw std::runtime_error(
            "Panel is empty."
        );
    }


    std::vector<OwnedPrimerPair>
        panel;


    while (
        std::getline(
            input,
            line
        )
    ) {
        if (line.empty()) {
            continue;
        }


        std::istringstream row(line);

        std::size_t pair_id = 0;

        std::string primer1;
        std::string primer2;

        std::uint64_t left = 0;
        std::uint64_t right = 0;
        std::uint64_t expected_length = 0;


        if (
            !(
                row
                >> pair_id
                >> primer1
                >> primer2
                >> left
                >> right
                >> expected_length
            )
        ) {
            throw std::runtime_error(
                "Malformed panel row."
            );
        }


        const auto valid_primer =
            [](
                const std::string& sequence
            ) {
                if (sequence.empty()) {
                    return false;
                }

                for (
                    const char base :
                    sequence
                ) {
                    if (
                        base != 'A' &&
                        base != 'C' &&
                        base != 'G' &&
                        base != 'T'
                    ) {
                        return false;
                    }
                }

                return true;
            };


        if (
            !valid_primer(primer1) ||
            !valid_primer(primer2)
        ) {
            throw std::runtime_error(
                "Panel contains non-ACGT "
                "primer at pair " +
                std::to_string(
                    pair_id
                )
            );
        }


        panel.push_back(
            OwnedPrimerPair{
                std::move(primer1),
                std::move(primer2)
            }
        );
    }


    if (panel.empty()) {
        throw std::runtime_error(
            "Panel contains no primer pairs."
        );
    }


    return panel;
}


std::unordered_map<
    std::string,
    ExpectedShard
>
load_expected(
    const std::string& path
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open expected summary: " +
            path
        );
    }

    std::string header;

    if (!std::getline(input, header)) {
        throw std::runtime_error(
            "Expected summary is empty."
        );
    }

    std::unordered_map<
        std::string,
        ExpectedShard
    > expected;

    std::string line;

    while (
        std::getline(
            input,
            line
        )
    ) {
        if (line.empty()) {
            continue;
        }

        std::istringstream row(line);

        std::string chromosome;

        std::uint64_t sequence_length = 0;

        std::size_t total_hits = 0;
        std::size_t forward_hits = 0;
        std::size_t reverse_hits = 0;

        std::size_t logical_requests = 0;
        std::size_t window_candidates = 0;
        std::size_t cross_products = 0;

        double load_ms = 0.0;
        double search_ms = 0.0;
        double v1_join_ms = 0.0;
        double v2_join_ms = 0.0;
        double speedup = 0.0;

        std::size_t rss_kb = 0;

        if (
            !(
                row
                >> chromosome
                >> sequence_length
                >> total_hits
                >> forward_hits
                >> reverse_hits
                >> logical_requests
                >> window_candidates
                >> cross_products
                >> load_ms
                >> search_ms
                >> v1_join_ms
                >> v2_join_ms
                >> speedup
                >> rss_kb
            )
        ) {
            throw std::runtime_error(
                "Malformed expected summary row."
            );
        }

        expected.emplace(
            chromosome,
            ExpectedShard{
                total_hits,
                window_candidates,
                cross_products
            }
        );
    }

    return expected;
}

}  // namespace


int main(
    const int argc,
    char** argv
) {
    try {
        using namespace primerpair;

        if (argc != 5) {
            std::cerr
                << "Usage:\n"
                << "validate_persistent_multiplex_v2_grch38 "
                << "<manifest.tsv> "
                << "<index_dir> "
                << "<panel.tsv> "
                << "<expected_summary.tsv>\n";

            return 2;
        }

        const auto manifest =
            PpfmManifest::load(
                argv[1],
                argv[2]
            );

        /*
         * MultiplexPrimerPairRequest contains string_views.
         *
         * Keep the owning primer strings alive until after
         * engine.search() has completed.
         */
        const auto owned_panel =
            load_panel(
                argv[3]
            );


        std::vector<
            MultiplexPrimerPairRequest
        > requests;


        requests.reserve(
            owned_panel.size()
        );


        for (
            const auto& pair :
            owned_panel
        ) {
            requests.push_back(
                MultiplexPrimerPairRequest{
                    pair.primer1,
                    pair.primer2,
                    3,
                    50,
                    3000
                }
            );
        }

        const auto expected =
            load_expected(
                argv[4]
            );

        PersistentMultiplexPrimerSearchEngineV2
            engine(
                manifest,
                1,
                8
            );

        const auto result =
            engine.search(
                requests,
                12,
                true,
                50,
                3000
            );

        if (
            result.shards.size() !=
            expected.size()
        ) {
            throw std::runtime_error(
                "Shard-count mismatch."
            );
        }

        std::size_t total_hits = 0;
        std::size_t total_candidates = 0;
        std::size_t total_cross = 0;

        for (
            const auto& shard :
            result.shards
        ) {
            const auto found =
                expected.find(
                    shard.chromosome
                );

            if (
                found ==
                expected.end()
            ) {
                throw std::runtime_error(
                    "Unexpected chromosome: " +
                    shard.chromosome
                );
            }

            const auto& exp =
                found->second;

            if (
                shard.stats.total_primer_hits !=
                exp.total_hits
            ) {
                throw std::runtime_error(
                    shard.chromosome +
                    ": total-hit mismatch."
                );
            }

            if (
                shard.global_cross_stats
                    .window_candidates
                !=
                exp.window_candidates
            ) {
                throw std::runtime_error(
                    shard.chromosome +
                    ": window-candidate mismatch."
                );
            }

            if (
                shard.cross_amplicons.size() !=
                exp.cross_products
            ) {
                throw std::runtime_error(
                    shard.chromosome +
                    ": cross-product mismatch."
                );
            }

            if (
                shard.intended_pairs.size() !=
                requests.size()
            ) {
                throw std::runtime_error(
                    shard.chromosome +
                    ": intended-pair count mismatch."
                );
            }

            total_hits +=
                shard.stats.total_primer_hits;

            total_candidates +=
                shard.global_cross_stats
                    .window_candidates;

            total_cross +=
                shard.cross_amplicons.size();

            std::cout
                << shard.chromosome
                << '\t'
                << shard.stats.total_primer_hits
                << '\t'
                << shard.global_cross_stats
                    .window_candidates
                << '\t'
                << shard.cross_amplicons.size()
                << '\t'
                << "YES"
                << '\n';
        }

        if (
            total_hits !=
            result.total_primer_hits()
        ) {
            throw std::runtime_error(
                "Whole-genome total-hit "
                "aggregation mismatch."
            );
        }

        if (
            total_candidates !=
            result.total_window_candidates()
        ) {
            throw std::runtime_error(
                "Whole-genome candidate "
                "aggregation mismatch."
            );
        }

        if (
            total_cross !=
            result.total_cross_amplicons()
        ) {
            throw std::runtime_error(
                "Whole-genome cross-product "
                "aggregation mismatch."
            );
        }

        std::cout
            << "validated_shards\t"
            << result.shards.size()
            << '\n';

        std::cout
            << "total_hits\t"
            << total_hits
            << '\n';

        std::cout
            << "total_window_candidates\t"
            << total_candidates
            << '\n';

        std::cout
            << "total_cross_amplicons\t"
            << total_cross
            << '\n';

        std::cout
            << "cache_loads\t"
            << engine.cache().load_count()
            << '\n';

        std::cout
            << "cache_evictions\t"
            << engine.cache().eviction_count()
            << '\n';

        std::cout
            << "VERIFY_PERSISTENT_MULTIPLEX_V2"
            << "\tYES\n";

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
