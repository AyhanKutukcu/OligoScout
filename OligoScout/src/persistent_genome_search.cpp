#include "primerpair/persistent_genome_search.hpp"

#include <stdexcept>
#include <utility>

namespace primerpair {

PersistentGenomeSearchEngine::
PersistentGenomeSearchEngine(
    PpfmManifest manifest,
    const std::size_t cache_capacity,
    const std::size_t suffix_array_sample_rate
)
    : manifest_(
          std::move(
              manifest
          )
      ),
      cache_(
          manifest_,
          cache_capacity,
          suffix_array_sample_rate
      ) {

    if (manifest_.empty()) {
        throw std::invalid_argument(
            "Persistent genome manifest "
            "cannot be empty."
        );
    }
}


GenomePairSearchResult
PersistentGenomeSearchEngine::search_pair(
    const std::string_view primer1,
    const std::string_view primer2,
    const std::size_t max_mismatches,
    const std::uint64_t min_amplicon_length,
    const std::uint64_t max_amplicon_length,
    const SensitivePairAnchorPolicy anchor_policy
) {
    GenomePairSearchResult result;

    result.shards.reserve(
        manifest_.size()
    );

    /*
     * Critical biological invariant:
     *
     * load chromosome
     * -> search chromosome
     * -> assemble pairs inside chromosome
     * -> store result
     * -> continue
     *
     * Primer hits from different chromosomes are
     * never pooled before pair assembly.
     */
    for (
        const auto& entry :
        manifest_.entries()
    ) {
        const GenomeShard& shard =
            cache_.get(
                entry.chromosome
            );

        auto shard_result =
            shard.search_pair(
                primer1,
                primer2,
                max_mismatches,
                min_amplicon_length,
                max_amplicon_length,
                anchor_policy
            );

        result.shards.push_back(
            GenomeShardPairSearchResult{
                shard.id(),
                shard.chromosome(),
                std::move(
                    shard_result
                )
            }
        );
    }

    return result;
}

}  // namespace primerpair
