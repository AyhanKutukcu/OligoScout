#include "primerpair/ppfm_io.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void expect(
    const bool condition,
    const char* message
) {
    if (!condition) {
        throw std::runtime_error(
            message
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}


std::string reverse_complement(
    const std::string& sequence
) {
    std::string result;

    result.reserve(
        sequence.size()
    );

    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {
        switch (*it) {
            case 'A':
                result.push_back('T');
                break;

            case 'C':
                result.push_back('G');
                break;

            case 'G':
                result.push_back('C');
                break;

            case 'T':
                result.push_back('A');
                break;

            default:
                throw std::runtime_error(
                    "Unexpected nucleotide."
                );
        }
    }

    return result;
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string primer1 =
            "ACGTACGTTGCAACGTACGT";

        const std::string primer2 =
            "TTGGAACCTTGGAACCTTGG";

        const std::string right_site =
            reverse_complement(
                primer2
            );

        const std::string sequence =
            std::string(
                50,
                'N'
            )
            +
            primer1
            +
            std::string(
                60,
                'N'
            )
            +
            right_site
            +
            std::string(
                50,
                'N'
            );


        const auto suffix =
            std::chrono::
                high_resolution_clock::
                now().
                time_since_epoch().
                count();

        const auto temp =
            std::filesystem::
                temp_directory_path()
            /
            (
                "primerpair_cache_test_" +
                std::to_string(
                    suffix
                )
            );

        std::filesystem::
            create_directories(
                temp
            );


        const std::string chromosomes[] = {
            "chrA",
            "chrB",
            "chrC"
        };


        for (
            const auto& chromosome :
            chromosomes
        ) {
            PackedReference reference(
                sequence
            );

            BidirectionalFMIndex index(
                sequence,
                8
            );

            PpfmIO::save_shard(
                temp /
                    (
                        chromosome +
                        ".sa8.ppfm"
                    ),
                chromosome,
                reference,
                index
            );
        }


        const auto manifest_path =
            temp /
            "manifest.sa8.tsv";

        {
            std::ofstream output(
                manifest_path
            );

            output
                << "chromosome\t"
                << "file_bytes\t"
                << "sha256\t"
                << "status\n";

            for (
                const auto& chromosome :
                chromosomes
            ) {
                const auto path =
                    temp /
                    (
                        chromosome +
                        ".sa8.ppfm"
                    );

                output
                    << chromosome
                    << '\t'
                    << std::filesystem::
                        file_size(
                            path
                        )
                    << '\t'
                    << std::string(
                        64,
                        '0'
                    )
                    << '\t'
                    << "BUILT\n";
            }
        }


        auto manifest =
            PpfmManifest::load(
                manifest_path,
                temp
            );

        PpfmShardCache cache(
            std::move(
                manifest
            ),
            2,
            8
        );


        expect(
            cache.capacity() == 2,
            "Cache capacity is two"
        );

        expect(
            cache.size() == 0,
            "Cache starts empty"
        );


        const auto& shard_a =
            cache.get(
                "chrA"
            );

        expect(
            shard_a.chromosome() ==
                "chrA",
            "chrA loaded lazily"
        );

        expect(
            cache.size() == 1,
            "One resident shard after chrA"
        );

        expect(
            cache.load_count() == 1,
            "First access performs one load"
        );


        (void)cache.get(
            "chrB"
        );

        expect(
            cache.size() == 2,
            "Two resident shards after chrB"
        );

        expect(
            cache.load_count() == 2,
            "Second chromosome performs second load"
        );


        (void)cache.get(
            "chrA"
        );

        expect(
            cache.hit_count() == 1,
            "Repeated chrA access is cache hit"
        );


        (void)cache.get(
            "chrC"
        );

        expect(
            cache.size() == 2,
            "Cache remains bounded at two"
        );

        expect(
            cache.contains(
                "chrA"
            ),
            "Recently touched chrA retained"
        );

        expect(
            !cache.contains(
                "chrB"
            ),
            "Least-recently-used chrB evicted"
        );

        expect(
            cache.contains(
                "chrC"
            ),
            "chrC resident after load"
        );

        expect(
            cache.eviction_count() == 1,
            "One LRU eviction recorded"
        );


        /*
         * Verify a cached persistent GenomeShard
         * is immediately usable for PCR pair search.
         */
        const auto& searchable_a =
            cache.get(
                "chrA"
            );

        const auto search =
            searchable_a.search_pair(
                primer1,
                primer2,
                3,
                100,
                100
            );

        expect(
            search
                .pair_result
                .amplicons
                .size()
            ==
            1,
            "Cached shard supports pair search"
        );

        const auto& hit =
            search
                .pair_result
                .amplicons
                .front();

        expect(
            hit.amplicon_start == 50,
            "Cached product start correct"
        );

        expect(
            hit.amplicon_end_exclusive ==
                150,
            "Cached product end correct"
        );

        expect(
            hit.amplicon_length == 100,
            "Cached product length correct"
        );

        expect(
            hit.total_mismatches() == 0,
            "Cached product is exact"
        );


        bool missing_rejected = false;

        try {
            (void)cache.get(
                "chrMissing"
            );

        } catch (
            const std::out_of_range&
        ) {
            missing_rejected = true;
        }

        expect(
            missing_rejected,
            "Missing chromosome rejected"
        );


        bool zero_capacity_rejected = false;

        try {
            auto second_manifest =
                PpfmManifest::load(
                    manifest_path,
                    temp
                );

            PpfmShardCache invalid(
                std::move(
                    second_manifest
                ),
                0,
                8
            );

        } catch (
            const std::invalid_argument&
        ) {
            zero_capacity_rejected = true;
        }

        expect(
            zero_capacity_rejected,
            "Zero cache capacity rejected"
        );


        std::filesystem::remove_all(
            temp
        );


        std::cout
            << "cache_loads\t"
            << cache.load_count()
            << '\n';

        std::cout
            << "cache_hits\t"
            << cache.hit_count()
            << '\n';

        std::cout
            << "cache_evictions\t"
            << cache.eviction_count()
            << '\n';

        std::cout
            << "ALL_CHECKS\tYES\n";

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';

        return 1;
    }
}
