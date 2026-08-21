#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/ppfm_io.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>


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


struct TemporaryDirectory {
    std::filesystem::path path;

    ~TemporaryDirectory() {
        std::error_code error;

        std::filesystem::remove_all(
            path,
            error
        );
    }
};


}  // namespace


int main() {
    using namespace primerpair;

    try {
        const auto suffix =
            std::chrono::
                high_resolution_clock::
                now().
                time_since_epoch().
                count();


        TemporaryDirectory temporary{
            std::filesystem::
                temp_directory_path()
            /
            (
                "primerpair_evict_before_load_" +
                std::to_string(
                    suffix
                )
            )
        };


        std::filesystem::
            create_directories(
                temporary.path
            );


        /*
         * One valid PPFM shard.
         */
        const std::string sequence =
            std::string(
                80,
                'N'
            )
            +
            "ACGTACGTTGCAACGTACGT"
            +
            std::string(
                80,
                'N'
            );


        {
            PackedReference reference(
                sequence
            );

            BidirectionalFMIndex index(
                sequence,
                8
            );

            PpfmIO::save_shard(
                temporary.path /
                    "chrA.sa8.ppfm",
                "chrA",
                reference,
                index
            );
        }


        /*
         * One file that passes manifest-level file
         * existence/size handling but is not a valid
         * PPFM binary.
         */
        const auto broken_path =
            temporary.path /
            "chrBad.sa8.ppfm";


        {
            std::ofstream output(
                broken_path,
                std::ios::binary
            );

            if (!output) {
                throw std::runtime_error(
                    "Cannot create broken PPFM."
                );
            }

            output
                << "THIS_IS_NOT_A_PPFM_SHARD";
        }


        const auto manifest_path =
            temporary.path /
            "manifest.sa8.tsv";


        {
            std::ofstream output(
                manifest_path
            );

            if (!output) {
                throw std::runtime_error(
                    "Cannot create manifest."
                );
            }


            output
                << "chromosome\t"
                << "file_bytes\t"
                << "sha256\t"
                << "status\n";


            const auto valid_path =
                temporary.path /
                "chrA.sa8.ppfm";


            output
                << "chrA"
                << '\t'
                << std::filesystem::
                    file_size(
                        valid_path
                    )
                << '\t'
                << std::string(
                    64,
                    '0'
                )
                << '\t'
                << "BUILT\n";


            output
                << "chrBad"
                << '\t'
                << std::filesystem::
                    file_size(
                        broken_path
                    )
                << '\t'
                << std::string(
                    64,
                    '0'
                )
                << '\t'
                << "BUILT\n";
        }


        auto manifest =
            PpfmManifest::load(
                manifest_path,
                temporary.path
            );


        expect(
            manifest.size() == 2,
            "Manifest contains valid and broken shards"
        );


        PpfmShardCache cache(
            std::move(
                manifest
            ),
            1,
            8
        );


        /*
         * Fill the capacity-one cache.
         */
        const auto& shard_a =
            cache.get(
                "chrA"
            );


        expect(
            shard_a.chromosome() ==
                "chrA",
            "Valid shard loads"
        );


        expect(
            cache.size() == 1,
            "Cache is full before replacement"
        );


        expect(
            cache.contains(
                "chrA"
            ),
            "chrA resident before failed replacement"
        );


        expect(
            cache.load_count() == 1,
            "One successful load recorded"
        );


        /*
         * New production contract:
         *
         * chrA is evicted BEFORE chrBad loading begins.
         * chrBad then fails PPFM deserialization.
         *
         * chrA must therefore remain evicted.
         */
        bool rejected = false;


        try {
            (void)cache.get(
                "chrBad"
            );

        } catch (
            const std::exception&
        ) {
            rejected = true;
        }


        expect(
            rejected,
            "Broken replacement shard rejected"
        );


        expect(
            cache.eviction_count() == 1,
            "LRU eviction occurred before failed load"
        );


        expect(
            cache.size() == 0,
            "Cache empty after failed replacement load"
        );


        expect(
            !cache.contains(
                "chrA"
            ),
            "Previous LRU shard is not restored"
        );


        expect(
            !cache.contains(
                "chrBad"
            ),
            "Broken shard is never inserted"
        );


        /*
         * load_count counts successful loads only.
         */
        expect(
            cache.load_count() == 1,
            "Failed replacement does not increment "
            "successful load count"
        );


        /*
         * Cache remains usable after the failure.
         */
        const auto& reloaded =
            cache.get(
                "chrA"
            );


        expect(
            reloaded.chromosome() ==
                "chrA",
            "Valid shard can be reloaded after failure"
        );


        expect(
            cache.size() == 1,
            "Cache operational after failed replacement"
        );


        expect(
            cache.load_count() == 2,
            "Successful reload increments load count"
        );


        expect(
            cache.eviction_count() == 1,
            "Reload from empty cache adds no eviction"
        );


        std::cout
            << "cache_loads\t"
            << cache.load_count()
            << '\n';

        std::cout
            << "cache_evictions\t"
            << cache.eviction_count()
            << '\n';

        std::cout
            << "EVICT_BEFORE_LOAD_CONTRACT\tYES\n";

        std::cout
            << "FAILED_LOAD_LEAVES_CACHE_VALID\tYES\n";

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
