#include "primerpair/persistent_genome_search.hpp"
#include "primerpair/ppfm_io.hpp"

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

        /*
         * One exact 100-bp product per chromosome.
         */
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
                "primerpair_persistent_genome_" +
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

        PersistentGenomeSearchEngine genome(
            std::move(
                manifest
            ),
            2,
            8
        );


        expect(
            genome.manifest().size() == 3,
            "Three persistent shards registered"
        );

        expect(
            genome.cache().capacity() == 2,
            "Whole-genome cache capacity is two"
        );

        expect(
            genome.cache().size() == 0,
            "Whole-genome cache starts empty"
        );


        const auto result =
            genome.search_pair(
                primer1,
                primer2,
                3,
                100,
                100
            );


        expect(
            result.shards.size() == 3,
            "All three shards searched"
        );

        expect(
            result.total_amplicon_count() == 3,
            "Three chromosome-local products found"
        );

        expect(
            genome.cache().size() == 2,
            "Cache remains bounded after genome scan"
        );

        expect(
            genome.cache().load_count() == 3,
            "Each chromosome loaded exactly once"
        );

        expect(
            genome.cache().eviction_count() == 1,
            "One eviction required for three-shard scan"
        );


        for (
            std::size_t i = 0;
            i < result.shards.size();
            ++i
        ) {
            const auto& shard =
                result.shards.at(i);

            expect(
                shard.shard_id == i,
                "Manifest-order shard id preserved"
            );

            expect(
                shard.chromosome ==
                    chromosomes[i],
                "Manifest chromosome order preserved"
            );

            const auto& hits =
                shard
                    .search_result
                    .pair_result
                    .amplicons;

            expect(
                hits.size() == 1,
                "Each chromosome has one product"
            );

            const auto& hit =
                hits.front();

            expect(
                hit.left_position == 50,
                "Whole-genome left position correct"
            );

            expect(
                hit.right_position == 130,
                "Whole-genome right position correct"
            );

            expect(
                hit.amplicon_start == 50,
                "Whole-genome amplicon start correct"
            );

            expect(
                hit.amplicon_end_exclusive == 150,
                "Whole-genome amplicon end correct"
            );

            expect(
                hit.amplicon_length == 100,
                "Whole-genome amplicon length correct"
            );

            expect(
                hit.total_mismatches() == 0,
                "Whole-genome product is exact"
            );
        }


        std::filesystem::remove_all(
            temp
        );


        std::cout
            << "shards_searched\t"
            << result.shards.size()
            << '\n';

        std::cout
            << "total_amplicons\t"
            << result.total_amplicon_count()
            << '\n';

        std::cout
            << "cache_loads\t"
            << genome.cache().load_count()
            << '\n';

        std::cout
            << "cache_evictions\t"
            << genome.cache().eviction_count()
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
