#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "primerpair/grch38_genome_loader.hpp"

namespace {


void expect(
    const bool condition,
    const std::string& message
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: "
            +
            message
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}


}  // namespace


int main() {
    try {
        const auto path =
            std::filesystem::
                temp_directory_path()
            /
            "primerpair_grch38_loader_test.fa";


        {
            std::ofstream output(
                path
            );

            output
                << ">NC_000021.9 Homo sapiens chromosome 21, "
                   "GRCh38.p14 Primary Assembly\n"
                << std::string(
                    1200,
                    'A'
                )
                << '\n'

                << ">NT_999999.1 Homo sapiens chromosome 21 "
                   "unlocalized genomic scaffold, "
                   "GRCh38.p14 Primary Assembly\n"
                << std::string(
                    5000,
                    'C'
                )
                << '\n'

                << ">NC_000022.11 Homo sapiens chromosome 22, "
                   "GRCh38.p14 Primary Assembly\n"
                << std::string(
                    1400,
                    'G'
                )
                << '\n'

                << ">NC_000023.11 Homo sapiens chromosome X, "
                   "GRCh38.p14 Primary Assembly\n"
                << std::string(
                    1600,
                    'T'
                )
                << '\n';
        }


        primerpair::GenomeSearchEngine
            genome(
                8
            );


        const auto summary =
            primerpair::
                load_grch38_primary_shards(
                    genome,
                    path.string(),
                    {
                        "chr21",
                        "chr22"
                    }
                );


        expect(
            summary.records_loaded == 2,
            "Two requested chromosomes loaded"
        );

        expect(
            summary.bases_loaded == 2600,
            "Loaded base count correct"
        );

        expect(
            summary.chromosomes.size() == 2,
            "Load summary contains two aliases"
        );

        expect(
            summary.chromosomes.at(0) ==
                "chr21",
            "chr21 loaded first"
        );

        expect(
            summary.chromosomes.at(1) ==
                "chr22",
            "chr22 loaded second"
        );

        expect(
            genome.shard_count() == 2,
            "Genome contains only selected shards"
        );

        expect(
            genome.shard(0).chromosome() ==
                "chr21",
            "Shard 0 is chr21"
        );

        expect(
            genome.shard(1).chromosome() ==
                "chr22",
            "Shard 1 is chr22"
        );

        expect(
            genome.shard(0)
                .sequence_length()
                ==
                1200,
            "chr21 sequence length preserved "
            "after source release"
        );

        expect(
            genome.shard(1)
                .sequence_length()
                ==
                1400,
            "chr22 sequence length preserved "
            "after source release"
        );


        bool duplicate_request_rejected =
            false;

        try {
            primerpair::GenomeSearchEngine
                duplicate_genome(
                    8
                );

            static_cast<void>(
                primerpair::
                    load_grch38_primary_shards(
                        duplicate_genome,
                        path.string(),
                        {
                            "chr21",
                            "chr21"
                        }
                    )
            );

        } catch (
            const std::invalid_argument&
        ) {
            duplicate_request_rejected =
                true;
        }

        expect(
            duplicate_request_rejected,
            "Duplicate requested chromosome rejected"
        );


        bool invalid_alias_rejected =
            false;

        try {
            primerpair::GenomeSearchEngine
                invalid_genome(
                    8
                );

            static_cast<void>(
                primerpair::
                    load_grch38_primary_shards(
                        invalid_genome,
                        path.string(),
                        {
                            "chr23"
                        }
                    )
            );

        } catch (
            const std::invalid_argument&
        ) {
            invalid_alias_rejected =
                true;
        }

        expect(
            invalid_alias_rejected,
            "Invalid canonical chromosome rejected"
        );


        bool missing_rejected =
            false;

        try {
            primerpair::GenomeSearchEngine
                missing_genome(
                    8
                );

            static_cast<void>(
                primerpair::
                    load_grch38_primary_shards(
                        missing_genome,
                        path.string(),
                        {
                            "chr1"
                        }
                    )
            );

        } catch (
            const std::runtime_error&
        ) {
            missing_rejected =
                true;
        }

        expect(
            missing_rejected,
            "Missing requested chromosome detected"
        );


        std::filesystem::remove(
            path
        );


        std::cout
            << "GRCh38 selective shard loader tests passed.\n";

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
