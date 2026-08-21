#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "primerpair/genome_search.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"

namespace primerpair {

/*
 * Whole-genome search over persistent PPFM shards.
 *
 * Shards are loaded lazily through PpfmShardCache,
 * searched independently, and may then be evicted.
 *
 * Pairing always remains chromosome-local.
 */
class PersistentGenomeSearchEngine {
public:
    PersistentGenomeSearchEngine(
        PpfmManifest manifest,
        std::size_t cache_capacity = 2,
        std::size_t suffix_array_sample_rate = 8
    );

    [[nodiscard]]
    GenomePairSearchResult search_pair(
        std::string_view primer1,
        std::string_view primer2,
        std::size_t max_mismatches = 3,
        std::uint64_t min_amplicon_length = 50,
        std::uint64_t max_amplicon_length = 3000,
        SensitivePairAnchorPolicy anchor_policy =
            SensitivePairAnchorPolicy::AdaptiveLowCost
    );

    [[nodiscard]]
    const PpfmManifest& manifest() const noexcept {
        return manifest_;
    }

    [[nodiscard]]
    const PpfmShardCache& cache() const noexcept {
        return cache_;
    }

private:
    PpfmManifest manifest_;
    PpfmShardCache cache_;
};

}  // namespace primerpair
