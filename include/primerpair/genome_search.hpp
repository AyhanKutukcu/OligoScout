#pragma once

#include <filesystem>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/fm_index.hpp"
#include "primerpair/genome_shard.hpp"

namespace primerpair {


struct GenomeShardPairSearchResult {
    std::size_t
        shard_id{0};

    std::string
        chromosome;

    SensitivePairConstrainedSearchResult
        search_result;
};


struct GenomePairSearchResult {
    std::vector<
        GenomeShardPairSearchResult
    > shards;

    [[nodiscard]]
    std::size_t
    shard_count() const noexcept {
        return shards.size();
    }

    [[nodiscard]]
    std::uint64_t
    total_amplicon_count() const noexcept;
};


class GenomeSearchEngine {
public:
    explicit GenomeSearchEngine(
        std::size_t suffix_array_sample_rate =
            FMIndex::kDefaultSuffixArraySampleRate
    );


    /*
     * Adds one independent chromosome/contig shard.
     *
     * Returns its stable zero-based shard id.
     */
    std::size_t add_shard(
        std::string chromosome,
        std::string sequence
    );


    [[nodiscard]]
    std::size_t
    shard_count() const noexcept {
        return shards_.size();
    }


    [[nodiscard]]
    std::size_t
    suffix_array_sample_rate() const noexcept {
        return suffix_array_sample_rate_;
    }
    /*
     * Load an already-built persistent chromosome
     * shard. No reference/FM-index reconstruction.
     */
    std::size_t load_ppfm_shard(
        const std::filesystem::path& path
    );




    [[nodiscard]]
    const GenomeShard&
    shard(
        std::size_t shard_id
    ) const;


    /*
     * Executes a complete pair search independently
     * on every shard.
     *
     * This is the key chromosome-isolation rule:
     * hits from different shards are never passed
     * into one common pairing kernel.
     */
    [[nodiscard]]
    GenomePairSearchResult search_pair(
        std::string_view primer1,
        std::string_view primer2,
        std::size_t max_mismatches = 3,
        std::uint64_t min_amplicon_length = 50,
        std::uint64_t max_amplicon_length = 3000,
        SensitivePairAnchorPolicy anchor_policy =
            SensitivePairAnchorPolicy::
                AdaptiveLowCost
    ) const;


private:
    std::size_t
        suffix_array_sample_rate_;

    std::vector<
        std::unique_ptr<GenomeShard>
    > shards_;
};

}  // namespace primerpair
