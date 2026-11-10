#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "primerpair/genome_search.hpp"

namespace primerpair {


struct Grch38ShardLoadSummary {
    std::size_t
        records_loaded{0};

    std::uint64_t
        bases_loaded{0};

    std::vector<std::string>
        chromosomes;
};


/*
 * Stream GRCh38.p14 FASTA and build only the
 * requested canonical chromosome shards.
 *
 * Example:
 *
 *   {"chr21", "chr22"}
 *
 * Unlocalized/alternate/patch records are ignored.
 *
 * Important:
 * rejected FASTA records are never accumulated
 * as sequence strings in RAM.
 */
[[nodiscard]]
Grch38ShardLoadSummary
load_grch38_primary_shards(
    GenomeSearchEngine& genome,
    const std::string& fasta_path,
    const std::vector<std::string>& chromosomes
);


}  // namespace primerpair
