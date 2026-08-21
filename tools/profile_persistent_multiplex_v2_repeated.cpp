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
#include <string_view>
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

    while (std::getline(input, line)) {
        std::istringstream row(line);

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
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open panel."
        );
    }

    std::string line;

    if (!std::getline(input, line)) {
        throw std::runtime_error(
            "Panel empty."
        );
    }

    std::vector<OwnedPrimerPair> result;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream row(line);

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
                std::move(primer1),
                std::move(primer2)
            }
        );
    }

    return result;
}


void report(
    const std::size_t iteration,
    const char* stage
) {
    const auto memory =
        memory_status();

    std::cout
        << "MEMORY"
        << '\t'
        << iteration
        << '\t'
        << stage
        << '\t'
        << "rss_kb="
        << memory.rss_kb
        << '\t'
        << "hwm_kb="
        << memory.hwm_kb
        << '\n';
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
                << "profile_persistent_multiplex_v2_repeated "
                << "<manifest.tsv> "
                << "<index_dir> "
                << "<panel.tsv> "
                << "<repetitions>\n";

            return 2;
        }


        const std::size_t repetitions =
            static_cast<std::size_t>(
                std::stoull(
                    argv[4]
                )
            );

        if (repetitions == 0) {
            throw std::invalid_argument(
                "Repetitions must be > 0."
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
        > requests;

        requests.reserve(
            owned.size()
        );


        for (const auto& pair : owned) {
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


        PersistentMultiplexPrimerSearchEngineV2
            engine(
                manifest,
                1,
                8
            );


        report(
            0,
            "before_queries"
        );


        std::size_t baseline_hits = 0;
        std::size_t baseline_candidates = 0;
        std::size_t baseline_cross = 0;


        for (
            std::size_t iteration = 1;
            iteration <= repetitions;
            ++iteration
        ) {
            const auto start =
                std::chrono::
                    steady_clock::
                    now();


            {
                const auto result =
                    engine.search(
                        requests,
                        12,
                        true,
                        50,
                        3000
                    );


                const std::size_t hits =
                    result.total_primer_hits();

                const std::size_t candidates =
                    result.total_window_candidates();

                const std::size_t cross =
                    result.total_cross_amplicons();


                if (iteration == 1) {
                    baseline_hits =
                        hits;

                    baseline_candidates =
                        candidates;

                    baseline_cross =
                        cross;

                } else {
                    if (
                        hits != baseline_hits ||
                        candidates != baseline_candidates ||
                        cross != baseline_cross
                    ) {
                        throw std::runtime_error(
                            "Repeated-query result mismatch."
                        );
                    }
                }


                report(
                    iteration,
                    "result_live"
                );


                std::cout
                    << "QUERY_RESULT"
                    << '\t'
                    << iteration
                    << '\t'
                    << "hits="
                    << hits
                    << '\t'
                    << "candidates="
                    << candidates
                    << '\t'
                    << "cross="
                    << cross
                    << '\n';
            }


            const auto end =
                std::chrono::
                    steady_clock::
                    now();


            const double milliseconds =
                std::chrono::duration<
                    double,
                    std::milli
                >(
                    end - start
                ).count();


            report(
                iteration,
                "result_destroyed"
            );


            std::cout
                << "QUERY_TIME_MS"
                << '\t'
                << iteration
                << '\t'
                << milliseconds
                << '\n';


            std::cout
                << "CACHE"
                << '\t'
                << iteration
                << '\t'
                << "loads="
                << engine.cache().load_count()
                << '\t'
                << "hits="
                << engine.cache().hit_count()
                << '\t'
                << "evictions="
                << engine.cache().eviction_count()
                << '\t'
                << "resident="
                << engine.cache().size()
                << '\n';
        }


        const std::uint64_t expected_loads =
            static_cast<std::uint64_t>(
                repetitions
            )
            *
            static_cast<std::uint64_t>(
                manifest.size()
            );


        const std::uint64_t expected_evictions =
            expected_loads > 0
            ?
            expected_loads - 1
            :
            0;


        std::cout
            << "expected_loads\t"
            << expected_loads
            << '\n';

        std::cout
            << "actual_loads\t"
            << engine.cache().load_count()
            << '\n';

        std::cout
            << "expected_evictions\t"
            << expected_evictions
            << '\n';

        std::cout
            << "actual_evictions\t"
            << engine.cache().eviction_count()
            << '\n';


        if (
            engine.cache().load_count() !=
            expected_loads
        ) {
            throw std::runtime_error(
                "Unexpected repeated-query load count."
            );
        }


        if (
            engine.cache().eviction_count() !=
            expected_evictions
        ) {
            throw std::runtime_error(
                "Unexpected repeated-query eviction count."
            );
        }


        if (
            baseline_hits != 675059 ||
            baseline_candidates != 289686 ||
            baseline_cross != 123565
        ) {
            throw std::runtime_error(
                "GRCh38 baseline totals changed."
            );
        }


        report(
            repetitions,
            "final"
        );


        std::cout
            << "REPEATED_RESULTS_STABLE\tYES\n";

        std::cout
            << "CACHE_COUNT_MODEL_MATCH\tYES\n";

        std::cout
            << "REPEATED_PROFILE_COMPLETE\tYES\n";


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
