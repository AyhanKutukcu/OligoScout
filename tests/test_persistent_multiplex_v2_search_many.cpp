#include "primerpair/persistent_multiplex_primer_search_v2.hpp"

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/ppfm_io.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <array>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>


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


void expect_equivalent(
    const primerpair::
        PersistentMultiplexSearchResultV2& batch,

    const primerpair::
        PersistentMultiplexSearchResultV2& scalar
) {
    expect(
        batch.shards.size() ==
            scalar.shards.size(),
        "Batch/scalar shard count equal"
    );

    expect(
        batch.total_primer_hits() ==
            scalar.total_primer_hits(),
        "Batch/scalar total hits equal"
    );

    expect(
        batch.total_window_candidates() ==
            scalar.total_window_candidates(),
        "Batch/scalar window candidates equal"
    );

    expect(
        batch.total_cross_amplicons() ==
            scalar.total_cross_amplicons(),
        "Batch/scalar cross-amplicons equal"
    );


    for (
        std::size_t i = 0;
        i < batch.shards.size();
        ++i
    ) {
        const auto& a =
            batch.shards.at(i);

        const auto& b =
            scalar.shards.at(i);


        expect(
            a.chromosome ==
                b.chromosome,
            "Chromosome ordering equal"
        );

        expect(
            a.sequence_length ==
                b.sequence_length,
            "Chromosome lengths equal"
        );

        expect(
            a.stats.total_primer_hits ==
                b.stats.total_primer_hits,
            "Per-shard primer hits equal"
        );

        expect(
            a.stats.forward_primer_hits ==
                b.stats.forward_primer_hits,
            "Per-shard forward hits equal"
        );

        expect(
            a.stats.reverse_primer_hits ==
                b.stats.reverse_primer_hits,
            "Per-shard reverse hits equal"
        );

        expect(
            a.global_cross_stats
                .window_candidates
            ==
            b.global_cross_stats
                .window_candidates,
            "Per-shard window candidates equal"
        );

        expect(
            a.cross_amplicons ==
                b.cross_amplicons,
            "Per-shard cross-products exact"
        );

        expect(
            a.intended_pairs.size() ==
                b.intended_pairs.size(),
            "Intended result-slot counts equal"
        );


        for (
            std::size_t p = 0;
            p < a.intended_pairs.size();
            ++p
        ) {
            expect(
                a.intended_pairs
                    .at(p)
                    .amplicons
                    .size()
                ==
                b.intended_pairs
                    .at(p)
                    .amplicons
                    .size(),
                "Per-pair intended counts equal"
            );
        }
    }
}


}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string primer_a =
            "ACGTTGCAACGATCGTACGA";

        const std::string primer_b =
            "TGCACCTAGGCTAACGTGCA";

        const std::string primer_c =
            "GATCCGATGCTAGTCAGTCA";

        const std::string primer_d =
            "CTAGGTCAGCATCGTACGTT";


        const std::string sequence_a =
            std::string(50, 'N')
            +
            primer_a
            +
            std::string(60, 'N')
            +
            reverse_complement(
                primer_d
            )
            +
            std::string(50, 'N');


        const std::string sequence_b =
            std::string(50, 'N')
            +
            primer_c
            +
            std::string(70, 'N')
            +
            reverse_complement(
                primer_b
            )
            +
            std::string(50, 'N');


        expect(
            sequence_a.size() == 200,
            "chrA length correct"
        );

        expect(
            sequence_b.size() == 210,
            "chrB length correct"
        );


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
                "primerpair_search_many_" +
                std::to_string(
                    suffix
                )
            )
        };


        std::filesystem::
            create_directories(
                temporary.path
            );


        const std::array<
            std::pair<
                std::string,
                std::string
            >,
            2
        > chromosomes{
            std::pair{
                std::string("chrA"),
                sequence_a
            },

            std::pair{
                std::string("chrB"),
                sequence_b
            }
        };


        for (
            const auto& [
                chromosome,
                sequence
            ] :
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
                temporary.path
                    /
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
            temporary.path
            /
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


            for (
                const auto& [
                    chromosome,
                    sequence
                ] :
                chromosomes
            ) {
                (void)sequence;

                const auto shard_path =
                    temporary.path
                    /
                    (
                        chromosome +
                        ".sa8.ppfm"
                    );


                output
                    << chromosome
                    << '\t'
                    << std::filesystem::
                        file_size(
                            shard_path
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


        const auto manifest =
            PpfmManifest::load(
                manifest_path,
                temporary.path
            );


        const std::vector<
            MultiplexPrimerPairRequest
        > panel{
            MultiplexPrimerPairRequest{
                primer_a,
                primer_b,
                0,
                80,
                150
            },

            MultiplexPrimerPairRequest{
                primer_c,
                primer_d,
                0,
                80,
                150
            }
        };


        const std::vector<
            std::vector<
                MultiplexPrimerPairRequest
            >
        > panels{
            panel,
            panel,
            panel
        };


        /*
         * Shard-major:
         *
         * 2 chromosomes total
         * => 2 loads, 1 eviction.
         */
        PersistentMultiplexPrimerSearchEngineV2
            batch_engine(
                manifest,
                1,
                8
            );


        const auto batch_results =
            batch_engine.search_many(
                panels,
                12,
                true,
                80,
                150
            );


        expect(
            batch_results.size() == 3,
            "search_many returns three results"
        );

        expect(
            batch_engine
                .cache()
                .load_count() == 2,
            "Shard-major loads two shards"
        );

        expect(
            batch_engine
                .cache()
                .eviction_count() == 1,
            "Shard-major performs one eviction"
        );


        /*
         * Query-major oracle:
         *
         * 3 queries x 2 chromosomes
         * => 6 loads, 5 evictions.
         */
        PersistentMultiplexPrimerSearchEngineV2
            scalar_engine(
                manifest,
                1,
                8
            );


        std::vector<
            PersistentMultiplexSearchResultV2
        > scalar_results;


        scalar_results.reserve(
            panels.size()
        );


        for (
            const auto& current :
            panels
        ) {
            scalar_results.push_back(
                scalar_engine.search(
                    current,
                    12,
                    true,
                    80,
                    150
                )
            );
        }


        expect(
            scalar_engine
                .cache()
                .load_count() == 6,
            "Query-major loads six shards"
        );

        expect(
            scalar_engine
                .cache()
                .eviction_count() == 5,
            "Query-major performs five evictions"
        );


        for (
            std::size_t i = 0;
            i < panels.size();
            ++i
        ) {
            expect_equivalent(
                batch_results.at(i),
                scalar_results.at(i)
            );

            expect(
                batch_results
                    .at(i)
                    .total_primer_hits() == 4,
                "Each query has four hits"
            );

            expect(
                batch_results
                    .at(i)
                    .total_window_candidates() == 2,
                "Each query has two window candidates"
            );

            expect(
                batch_results
                    .at(i)
                    .total_cross_amplicons() == 2,
                "Each query has two cross-amplicons"
            );
        }


        std::cout
            << "batch_loads\t"
            << batch_engine
                .cache()
                .load_count()
            << '\n';

        std::cout
            << "scalar_loads\t"
            << scalar_engine
                .cache()
                .load_count()
            << '\n';

        std::cout
            << "batch_evictions\t"
            << batch_engine
                .cache()
                .eviction_count()
            << '\n';

        std::cout
            << "scalar_evictions\t"
            << scalar_engine
                .cache()
                .eviction_count()
            << '\n';

        std::cout
            << "SEARCH_MANY_EXACT_EQUIVALENCE\tYES\n";

        std::cout
            << "SHARD_MAJOR_LOAD_REDUCTION\tYES\n";

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
