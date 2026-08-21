#include "primerpair/ppfm_shard_cache.hpp"

#include "primerpair/ppfm_io.hpp"

#include <stdexcept>
#include <utility>

namespace primerpair {

PpfmShardCache::PpfmShardCache(
    PpfmManifest manifest,
    const std::size_t capacity,
    const std::size_t suffix_array_sample_rate
)
    : manifest_(
          std::move(
              manifest
          )
      ),
      capacity_(
          capacity
      ),
      suffix_array_sample_rate_(
          suffix_array_sample_rate
      ) {

    if (capacity_ == 0) {
        throw std::invalid_argument(
            "PPFM shard cache capacity "
            "must be > 0."
        );
    }

    if (
        suffix_array_sample_rate_ == 0
    ) {
        throw std::invalid_argument(
            "PPFM shard cache SA sample "
            "rate must be > 0."
        );
    }

    if (manifest_.empty()) {
        throw std::invalid_argument(
            "PPFM shard cache manifest "
            "cannot be empty."
        );
    }
}


void PpfmShardCache::touch(
    std::unordered_map<
        std::string,
        CachedShard
    >::iterator it
) {
    lru_.erase(
        it->second.lru_position
    );

    lru_.push_front(
        it->first
    );

    it->second.lru_position =
        lru_.begin();
}


bool PpfmShardCache::contains(
    const std::string_view chromosome
) const noexcept {
    return
        cache_.find(
            std::string(
                chromosome
            )
        )
        !=
        cache_.end();
}


const GenomeShard&
PpfmShardCache::get(
    const std::string_view chromosome
) {
    const std::string key(
        chromosome
    );

    auto cached =
        cache_.find(
            key
        );

    if (
        cached !=
        cache_.end()
    ) {
        ++hit_count_;

        touch(
            cached
        );

        return
            *cached->second.shard;
    }


    const PpfmManifestEntry* entry =
        manifest_.find(
            chromosome
        );

    if (entry == nullptr) {
        throw std::out_of_range(
            "Chromosome not present in "
            "PPFM manifest: " +
            key
        );
    }


    /*
     * Memory-bounded miss handling.
     *
     * Evict the least-recently-used resident shard
     * BEFORE loading the replacement.
     *
     * This is important for chromosome-scale PPFM
     * shards: load-before-evict temporarily keeps both
     * the old and new chromosome resident and can nearly
     * double process peak RSS even when capacity_ == 1.
     *
     * Contract note:
     * if the subsequent PPFM load fails, the evicted
     * shard is not restored automatically.
     */
    if (
        cache_.size() >=
        capacity_
    ) {
        if (lru_.empty()) {
            throw std::logic_error(
                "PPFM LRU cache is empty while "
                "resident cache is at capacity."
            );
        }

        const std::string victim =
            lru_.back();

        lru_.pop_back();

        const std::size_t erased =
            cache_.erase(
                victim
            );

        if (erased != 1) {
            throw std::logic_error(
                "PPFM LRU cache internal "
                "state is inconsistent."
            );
        }

        ++eviction_count_;
    }


    auto loaded =
        PpfmIO::load_shard(
            entry->path
        );

    if (
        loaded.chromosome !=
        entry->chromosome
    ) {
        throw std::runtime_error(
            "PPFM internal chromosome does "
            "not match manifest: " +
            entry->chromosome
        );
    }


    const std::size_t forward_rate =
        loaded.index
            .forward_index()
            .suffix_array_sample_rate();

    const std::size_t reverse_rate =
        loaded.index
            .reverse_index()
            .suffix_array_sample_rate();

    if (
        forward_rate !=
            suffix_array_sample_rate_
        ||
        reverse_rate !=
            suffix_array_sample_rate_
    ) {
        throw std::runtime_error(
            "PPFM shard SA sample rate "
            "does not match cache."
        );
    }


    const auto& entries =
        manifest_.entries();

    const std::size_t shard_id =
        static_cast<std::size_t>(
            entry -
            entries.data()
        );


    auto new_shard =
        std::make_unique<GenomeShard>(
            shard_id,
            std::move(
                loaded.chromosome
            ),
            std::move(
                loaded.reference
            ),
            std::move(
                loaded.index
            )
        );


    lru_.push_front(
        key
    );

    auto [
        inserted,
        success
    ] =
        cache_.emplace(
            key,
            CachedShard{
                std::move(
                    new_shard
                ),
                lru_.begin()
            }
        );

    if (!success) {
        lru_.pop_front();

        throw std::logic_error(
            "PPFM cache insertion failed."
        );
    }

    ++load_count_;

    return
        *inserted->second.shard;
}

}  // namespace primerpair
