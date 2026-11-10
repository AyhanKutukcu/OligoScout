#pragma once

#include <cstddef>
#include <cstdint>
#include <list>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>

#include "primerpair/genome_shard.hpp"
#include "primerpair/ppfm_manifest.hpp"

namespace primerpair {

/*
 * Lazy bounded LRU cache for persistent PPFM shards.
 *
 * Important:
 * A reference returned by get() remains valid only
 * while that chromosome remains resident in the
 * cache. Do not retain the reference across later
 * get() calls that may trigger eviction.
 *
 * Memory policy:
 * On a cache miss when the cache is full, the LRU
 * shard is evicted before the replacement PPFM shard
 * is loaded. This bounds chromosome-scale transient
 * memory usage and avoids old-shard + new-shard peak
 * residency.
 *
 * Consequently, if replacement loading fails after
 * eviction, the previous LRU shard is not restored
 * automatically.
 */
class PpfmShardCache {
public:
    PpfmShardCache(
        PpfmManifest manifest,
        std::size_t capacity,
        std::size_t suffix_array_sample_rate = 8
    );

    [[nodiscard]]
    const GenomeShard& get(
        std::string_view chromosome
    );

    [[nodiscard]]
    bool contains(
        std::string_view chromosome
    ) const noexcept;

    [[nodiscard]]
    std::size_t size() const noexcept {
        return cache_.size();
    }

    [[nodiscard]]
    std::size_t capacity() const noexcept {
        return capacity_;
    }

    [[nodiscard]]
    std::uint64_t hit_count() const noexcept {
        return hit_count_;
    }

    [[nodiscard]]
    std::uint64_t load_count() const noexcept {
        return load_count_;
    }

    [[nodiscard]]
    std::uint64_t eviction_count() const noexcept {
        return eviction_count_;
    }

private:
    using LruList =
        std::list<std::string>;

    struct CachedShard {
        std::unique_ptr<GenomeShard> shard;

        LruList::iterator
            lru_position;
    };

    void touch(
        std::unordered_map<
            std::string,
            CachedShard
        >::iterator it
    );

    PpfmManifest manifest_;

    std::size_t capacity_{0};

    std::size_t
        suffix_array_sample_rate_{8};

    /*
     * Front = most recently used.
     * Back  = least recently used.
     */
    LruList lru_;

    std::unordered_map<
        std::string,
        CachedShard
    > cache_;

    std::uint64_t hit_count_{0};

    std::uint64_t load_count_{0};

    std::uint64_t eviction_count_{0};
};

}  // namespace primerpair
