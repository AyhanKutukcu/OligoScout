#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "primerpair/global_multiplex_cross_join.hpp"
#include "primerpair/multiplex_primer_search.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"
#include "primerpair/primer_pair_search.hpp"

namespace primerpair {

/*
 * Statistics for one persistent chromosome shard.
 *
 * Primer searching and all pair assembly are completed
 * inside this shard before the next chromosome is loaded.
 */
struct PersistentMultiplexShardStatsV2 {
    std::size_t pair_requests{0};
    std::size_t total_primer_slots{0};

    std::size_t unique_primer_queries{0};
    std::size_t reused_primer_slots{0};

    std::size_t total_primer_hits{0};
    std::size_t forward_primer_hits{0};
    std::size_t reverse_primer_hits{0};

    std::size_t intended_join_computations{0};

    /*
     * Logical V1 comparison count:
     *
     *   4 * C(P, 2)
     *
     * No such enumeration is performed by V2.
     */
    std::size_t logical_cross_slot_requests{0};

    std::size_t cross_amplicon_records{0};
};


struct PersistentMultiplexShardResultV2 {
    std::size_t shard_id{0};

    std::string chromosome;

    std::uint64_t sequence_length{0};

    std::vector<
        PrimerPairSearchResult
    > intended_pairs;

    std::vector<
        MultiplexCrossAmplicon
    > cross_amplicons;

    PersistentMultiplexShardStatsV2 stats{};

    GlobalMultiplexCrossJoinStats
        global_cross_stats{};
};


struct PersistentMultiplexSearchResultV2 {
    std::vector<
        PersistentMultiplexShardResultV2
    > shards;


    [[nodiscard]]
    std::size_t
    total_cross_amplicons() const noexcept {
        std::size_t total = 0;

        for (const auto& shard : shards) {
            total +=
                shard.cross_amplicons.size();
        }

        return total;
    }


    [[nodiscard]]
    std::size_t
    total_primer_hits() const noexcept {
        std::size_t total = 0;

        for (const auto& shard : shards) {
            total +=
                shard.stats.total_primer_hits;
        }

        return total;
    }


    [[nodiscard]]
    std::size_t
    total_window_candidates() const noexcept {
        std::size_t total = 0;

        for (const auto& shard : shards) {
            total +=
                shard.global_cross_stats
                    .window_candidates;
        }

        return total;
    }
};


/*
 * Whole-GRCh38 multiplex search over persistent PPFM
 * chromosome shards.
 *
 * Biological invariant:
 *
 *   load chromosome
 *   -> search all unique primers
 *   -> assemble intended pairs
 *   -> global multiplex sweep
 *   -> materialize chromosome-labelled result
 *   -> continue to next chromosome
 *
 * Hits from different chromosomes are NEVER pooled.
 */
class PersistentMultiplexPrimerSearchEngineV2 {
public:
    PersistentMultiplexPrimerSearchEngineV2(
        PpfmManifest manifest,
        std::size_t cache_capacity = 2,
        std::size_t suffix_array_sample_rate = 8
    );


    [[nodiscard]]
    PersistentMultiplexSearchResultV2
    search(
        const std::vector<
            MultiplexPrimerPairRequest
        >& requests,

        std::size_t anchor_length = 12,

        bool include_cross_pairs = true,

        std::uint64_t
            cross_min_amplicon_length = 50,

        std::uint64_t
            cross_max_amplicon_length = 3000
    );


    [[nodiscard]]
    const PpfmManifest&
    manifest() const noexcept {
        return manifest_;
    }


    [[nodiscard]]
    const PpfmShardCache&
    cache() const noexcept {
        return cache_;
    }



    /*
     * Shard-major multi-panel execution.
     *
     * Repeated search() calls are query-major:
     *
     *   query1 -> chr1..chrY
     *   query2 -> chr1..chrY
     *
     * search_many() is shard-major:
     *
     *   chr1 -> all panels
     *   chr2 -> all panels
     *   ...
     *
     * This amortizes persistent PPFM loading across
     * queued whole-genome multiplex queries.
     *
     * Each output element corresponds to the input
     * panel at the same index.
     */
    [[nodiscard]]
    std::vector<
        PersistentMultiplexSearchResultV2
    >
    search_many(
        const std::vector<
            std::vector<
                MultiplexPrimerPairRequest
            >
        >& panels,

        std::size_t anchor_length = 12,

        bool include_cross_pairs = true,

        std::uint64_t
            cross_min_amplicon_length = 50,

        std::uint64_t
            cross_max_amplicon_length = 3000
    );


private:
    PpfmManifest manifest_;
    PpfmShardCache cache_;
};

}  // namespace primerpair
