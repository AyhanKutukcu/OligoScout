#include "primerpair/persistent_multiplex_primer_search_v2.hpp"
#include "primerpair/ppfm_manifest.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace {


struct OwnedPrimerPair {
    std::string primer1;
    std::string primer2;
};


struct MemoryStatus {
    std::uint64_t rss_kb{0};
    std::uint64_t hwm_kb{0};
};


MemoryStatus memory_status() {
    std::ifstream input(
        "/proc/self/status"
    );

    if (!input) {
        throw std::runtime_error(
            "Cannot read /proc/self/status."
        );
    }


    MemoryStatus result;

    std::string line;


    while (
        std::getline(
            input,
            line
        )
    ) {
        std::istringstream row(
            line
        );

        std::string key;

        row >> key;


        if (key == "VmRSS:") {
            row >> result.rss_kb;
        }

        if (key == "VmHWM:") {
            row >> result.hwm_kb;
        }
    }


    return result;
}


std::vector<OwnedPrimerPair>
load_panel(
    const std::string& path
) {
    std::ifstream input(
        path
    );

    if (!input) {
        throw std::runtime_error(
            "Cannot open panel."
        );
    }


    std::string line;

    if (
        !std::getline(
            input,
            line
        )
    ) {
        throw std::runtime_error(
            "Panel is empty."
        );
    }


    std::vector<
        OwnedPrimerPair
    > result;


    while (
        std::getline(
            input,
            line
        )
    ) {
        if (line.empty()) {
            continue;
        }


        std::istringstream row(
            line
        );

        std::size_t pair_id = 0;

        std::string primer1;
        std::string primer2;

        std::uint64_t left = 0;
        std::uint64_t right = 0;
        std::uint64_t length = 0;


        if (
            !(
                row
                >> pair_id
                >> primer1
                >> primer2
                >> left
                >> right
                >> length
            )
        ) {
            throw std::runtime_error(
                "Malformed panel."
            );
        }


        result.push_back(
            OwnedPrimerPair{
                std::move(
                    primer1
                ),
                std::move(
                    primer2
                )
            }
        );
    }


    return result;
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
                << "Usage: "
                << "profile_persistent_multiplex_v2_search_many "
                << "<manifest.tsv> "
                << "<index_dir> "
                << "<panel.tsv> "
                << "<query_count>\n";

            return 2;
        }


        const std::size_t query_count =
            static_cast<std::size_t>(
                std::stoull(
                    argv[4]
                )
            );


        if (query_count == 0) {
            throw std::invalid_argument(
                "Query count must be > 0."
            );
        }


        const auto manifest =
            PpfmManifest::load(
                argv[1],
                argv[2]
            );


        const auto owned =
            load_panel(
                argv[3]
            );


        std::vector<
            MultiplexPrimerPairRequest
        > panel;


        panel.reserve(
            owned.size()
        );


        for (
            const auto& pair :
            owned
        ) {
            panel.push_back(
                MultiplexPrimerPairRequest{
                    pair.primer1,
                    pair.primer2,
                    3,
                    50,
                    3000
                }
            );
        }


        std::vector<
            std::vector<
                MultiplexPrimerPairRequest
            >
        > panels;


        panels.reserve(
            query_count
        );


        for (
            std::size_t i = 0;
            i < query_count;
            ++i
        ) {
            panels.push_back(
                panel
            );
        }


        PersistentMultiplexPrimerSearchEngineV2
            engine(
                manifest,
                1,
                8
            );


        const auto before =
            memory_status();


        std::cout
            << "MEMORY_BEFORE"
            << '\t'
            << "rss_kb="
            << before.rss_kb
            << '\t'
            << "hwm_kb="
            << before.hwm_kb
            << '\n';


        const auto start =
            std::chrono::
                steady_clock::
                now();


        const auto results =
            engine.search_many(
                panels,
                12,
                true,
                50,
                3000
            );


        const auto end =
            std::chrono::
                steady_clock::
                now();


        const double elapsed_ms =
            std::chrono::duration<
                double,
                std::milli
            >(
                end - start
            ).count();


        if (
            results.size() !=
            query_count
        ) {
            throw std::runtime_error(
                "search_many result count mismatch."
            );
        }


        for (
            std::size_t i = 0;
            i < results.size();
            ++i
        ) {
            const auto& result =
                results.at(i);


            if (
                result.total_primer_hits() !=
                    675059
                ||
                result.total_window_candidates() !=
                    289686
                ||
                result.total_cross_amplicons() !=
                    123565
            ) {
                throw std::runtime_error(
                    "GRCh38 search_many totals changed."
                );
            }


            std::cout
                << "QUERY_RESULT"
                << '\t'
                << (i + 1)
                << '\t'
                << "hits="
                << result.total_primer_hits()
                << '\t'
                << "candidates="
                << result
                    .total_window_candidates()
                << '\t'
                << "cross="
                << result
                    .total_cross_amplicons()
                << '\n';
        }


        const std::uint64_t expected_loads =
            static_cast<std::uint64_t>(
                manifest.size()
            );


        const std::uint64_t expected_evictions =
            expected_loads > 0
            ?
            expected_loads - 1
            :
            0;


        if (
            engine.cache().load_count() !=
            expected_loads
        ) {
            throw std::runtime_error(
                "Shard-major load count mismatch."
            );
        }


        if (
            engine.cache().eviction_count() !=
            expected_evictions
        ) {
            throw std::runtime_error(
                "Shard-major eviction count mismatch."
            );
        }


        const auto after =
            memory_status();


        std::cout
            << "query_count\t"
            << query_count
            << '\n';

        std::cout
            << "shard_count\t"
            << manifest.size()
            << '\n';

        std::cout
            << "elapsed_ms\t"
            << elapsed_ms
            << '\n';

        std::cout
            << "expected_loads\t"
            << expected_loads
            << '\n';

        std::cout
            << "actual_loads\t"
            << engine
                .cache()
                .load_count()
            << '\n';

        std::cout
            << "cache_hits\t"
            << engine
                .cache()
                .hit_count()
            << '\n';

        std::cout
            << "expected_evictions\t"
            << expected_evictions
            << '\n';

        std::cout
            << "actual_evictions\t"
            << engine
                .cache()
                .eviction_count()
            << '\n';

        std::cout
            << "MEMORY_AFTER"
            << '\t'
            << "rss_kb="
            << after.rss_kb
            << '\t'
            << "hwm_kb="
            << after.hwm_kb
            << '\n';

        std::cout
            << "SEARCH_MANY_GRCH38_STABLE\tYES\n";

        std::cout
            << "SHARD_MAJOR_CACHE_MODEL\tYES\n";

        std::cout
            << "SEARCH_MANY_PROFILE_COMPLETE\tYES\n";


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
