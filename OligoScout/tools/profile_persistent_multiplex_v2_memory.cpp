#include "primerpair/persistent_multiplex_primer_search_v2.hpp"
#include "primerpair/ppfm_manifest.hpp"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>


namespace {


struct OwnedPrimerPair {
    std::string primer1;
    std::string primer2;
};


struct MemoryStatus {
    std::uint64_t vm_rss_kb{0};
    std::uint64_t vm_hwm_kb{0};
};


MemoryStatus read_memory_status() {
    std::ifstream input(
        "/proc/self/status"
    );

    if (!input) {
        throw std::runtime_error(
            "Cannot read /proc/self/status."
        );
    }


    MemoryStatus status;

    std::string line;


    while (
        std::getline(
            input,
            line
        )
    ) {
        std::istringstream row(line);

        std::string key;

        row >> key;


        if (key == "VmRSS:") {
            row >> status.vm_rss_kb;
        }

        if (key == "VmHWM:") {
            row >> status.vm_hwm_kb;
        }
    }


    return status;
}


void report_memory(
    const std::string_view stage
) {
    const auto memory =
        read_memory_status();

    std::cout
        << "MEMORY"
        << '\t'
        << stage
        << '\t'
        << "rss_kb="
        << memory.vm_rss_kb
        << '\t'
        << "hwm_kb="
        << memory.vm_hwm_kb
        << '\n';
}


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


        panel.push_back(
            OwnedPrimerPair{
                std::move(primer1),
                std::move(primer2)
            }
        );
    }


    return panel;
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
                << "profile_persistent_multiplex_v2_memory "
                << "<manifest.tsv> "
                << "<index_dir> "
                << "<panel.tsv> "
                << "<full|intended-only>\n";

            return 2;
        }


        const std::string mode =
            argv[4];


        const bool include_cross_pairs =
            mode == "full";


        if (
            mode != "full" &&
            mode != "intended-only"
        ) {
            throw std::invalid_argument(
                "Mode must be full or intended-only."
            );
        }


        const auto manifest =
            PpfmManifest::load(
                argv[1],
                argv[2]
            );


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


        std::cout
            << "mode\t"
            << mode
            << '\n';

        std::cout
            << "pair_count\t"
            << requests.size()
            << '\n';

        std::cout
            << "manifest_shards\t"
            << manifest.size()
            << '\n';


        report_memory(
            "process_start"
        );


        {
            PersistentMultiplexPrimerSearchEngineV2
                engine(
                    manifest,
                    1,
                    8
                );


            report_memory(
                "engine_constructed"
            );


            auto result =
                engine.search(
                    requests,
                    12,
                    include_cross_pairs,
                    50,
                    3000
                );


            report_memory(
                "after_search"
            );


            std::size_t intended_products = 0;

            for (
                const auto& shard :
                result.shards
            ) {
                for (
                    const auto& intended :
                    shard.intended_pairs
                ) {
                    intended_products +=
                        intended
                            .amplicons
                            .size();
                }
            }


            std::cout
                << "result_shards\t"
                << result.shards.size()
                << '\n';

            std::cout
                << "total_primer_hits\t"
                << result.total_primer_hits()
                << '\n';

            std::cout
                << "intended_amplicons\t"
                << intended_products
                << '\n';

            std::cout
                << "cross_amplicons\t"
                << result.total_cross_amplicons()
                << '\n';

            std::cout
                << "window_candidates\t"
                << result.total_window_candidates()
                << '\n';

            std::cout
                << "cache_size\t"
                << engine.cache().size()
                << '\n';

            std::cout
                << "cache_loads\t"
                << engine.cache().load_count()
                << '\n';

            std::cout
                << "cache_evictions\t"
                << engine.cache().eviction_count()
                << '\n';


            /*
             * Release all retained result vectors.
             */
            result.shards.clear();
            result.shards.shrink_to_fit();


            report_memory(
                "after_result_clear"
            );


            const int trim_result =
                malloc_trim(0);


            std::cout
                << "malloc_trim_after_result\t"
                << trim_result
                << '\n';


            report_memory(
                "after_result_trim"
            );
        }


        /*
         * Engine destruction also releases the final
         * resident PPFM shard from capacity-one cache.
         */
        report_memory(
            "after_engine_destroy"
        );


        const int final_trim =
            malloc_trim(0);


        std::cout
            << "malloc_trim_after_engine\t"
            << final_trim
            << '\n';


        report_memory(
            "after_engine_trim"
        );


        std::cout
            << "MEMORY_PROFILE_COMPLETE\tYES\n";


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
