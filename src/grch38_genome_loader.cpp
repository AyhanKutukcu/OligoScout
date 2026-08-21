#include "primerpair/grch38_genome_loader.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "primerpair/fasta_reader.hpp"
#include "primerpair/grch38_reference.hpp"

namespace primerpair {

namespace {


bool canonical_alias(
    const std::string& chromosome
) {
    if (
        chromosome == "chrX"
        ||
        chromosome == "chrY"
    ) {
        return true;
    }

    if (
        !chromosome.starts_with(
            "chr"
        )
    ) {
        return false;
    }

    const std::string number =
        chromosome.substr(
            3
        );

    if (
        number.empty()
        ||
        number.size() > 2
    ) {
        return false;
    }

    unsigned int value = 0;

    for (const char c : number) {
        if (
            c < '0'
            ||
            c > '9'
        ) {
            return false;
        }

        value =
            value * 10u
            +
            static_cast<unsigned int>(
                c - '0'
            );
    }

    return
        value >= 1
        &&
        value <= 22;
}


}  // namespace


Grch38ShardLoadSummary
load_grch38_primary_shards(
    GenomeSearchEngine& genome,
    const std::string& fasta_path,
    const std::vector<std::string>& chromosomes
) {
    if (chromosomes.empty()) {
        throw std::invalid_argument(
            "At least one GRCh38 chromosome "
            "must be requested."
        );
    }


    std::unordered_set<std::string>
        requested;

    requested.reserve(
        chromosomes.size() * 2
    );


    for (const auto& chromosome : chromosomes) {

        if (
            !canonical_alias(
                chromosome
            )
        ) {
            throw std::invalid_argument(
                "Invalid canonical chromosome alias: "
                +
                chromosome
            );
        }

        if (
            !requested.insert(
                chromosome
            ).second
        ) {
            throw std::invalid_argument(
                "Duplicate requested chromosome: "
                +
                chromosome
            );
        }
    }


    Grch38ShardLoadSummary
        summary;


    const std::size_t selected =
        stream_selected_fasta_records(

            fasta_path,

            [&requested](
                const std::string_view name,
                const std::string_view description
            ) {
                const auto alias =
                    grch38_primary_chromosome_alias(
                        name,
                        description
                    );

                return
                    alias
                    &&
                    requested.contains(
                        *alias
                    );
            },

            [&](
                FastaRecord&& record
            ) {
                const auto alias =
                    grch38_primary_chromosome_alias(
                        record.name,
                        record.description
                    );

                if (!alias) {
                    throw std::logic_error(
                        "Selected GRCh38 record "
                        "lost chromosome alias."
                    );
                }

                const std::uint64_t length =
                    static_cast<std::uint64_t>(
                        record.sequence.size()
                    );

                genome.add_shard(
                    *alias,
                    std::move(
                        record.sequence
                    )
                );

                summary.records_loaded +=
                    1;

                summary.bases_loaded +=
                    length;

                summary.chromosomes.push_back(
                    *alias
                );
            }
        );


    if (
        selected !=
        chromosomes.size()
        ||
        summary.records_loaded !=
        chromosomes.size()
    ) {
        throw std::runtime_error(
            "GRCh38 FASTA did not contain every "
            "requested canonical chromosome."
        );
    }


    std::unordered_set<std::string>
        loaded(
            summary.chromosomes.begin(),
            summary.chromosomes.end()
        );


    for (const auto& chromosome : chromosomes) {

        if (
            !loaded.contains(
                chromosome
            )
        ) {
            throw std::runtime_error(
                "Requested GRCh38 chromosome "
                "was not loaded: "
                +
                chromosome
            );
        }
    }


    return summary;
}


}  // namespace primerpair
