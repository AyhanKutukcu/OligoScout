#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/sensitive_pair_constrained_search.hpp"

namespace primerpair {

/*
 * One independently indexed genomic sequence.
 *
 * Typical production use:
 *   shard 0  -> chr1
 *   shard 1  -> chr2
 *   ...
 *
 * Search/pairing never crosses shard boundaries.
 */
class GenomeShard {
public:
    GenomeShard(
        std::size_t id,
        std::string chromosome,
        std::string sequence,
        std::size_t suffix_array_sample_rate =
            FMIndex::kDefaultSuffixArraySampleRate
    );

    GenomeShard(
        const GenomeShard&
    ) = delete;

    GenomeShard& operator=(
        const GenomeShard&
    ) = delete;

    GenomeShard(
        GenomeShard&&
    ) = delete;

    GenomeShard& operator=(
        GenomeShard&&
    ) = delete;


    [[nodiscard]]
    std::size_t id() const noexcept {
        return id_;
    }


    [[nodiscard]]
    const std::string&
    chromosome() const noexcept {
        return chromosome_;
    }


    [[nodiscard]]
    std::uint64_t
    sequence_length() const noexcept {
        return sequence_length_;
    }


    [[nodiscard]]
    std::size_t
    suffix_array_sample_rate() const noexcept {
        return index_
            .forward_index()
            .suffix_array_sample_rate();
    }


    [[nodiscard]]
    const PackedReference&
    reference() const noexcept {
        return reference_;
    }


    [[nodiscard]]
    const BidirectionalFMIndex&
    index() const noexcept {
        return index_;
    }
    /*
     * Persistent-index constructor.
     *
     * The packed reference and both FM indexes are
     * already constructed. They are moved directly
     * into the shard; no FASTA/SA/BWT rebuild occurs.
     */
    GenomeShard(
        std::size_t id,
        std::string chromosome,
        PackedReference reference,
        BidirectionalFMIndex index
    );







    [[nodiscard]]
    SensitivePairConstrainedSearchResult
    search_pair(
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
        id_{0};

    std::string
        chromosome_;

    /*
     * Keep the normalized source sequence alive.
     * This also gives us explicit per-shard length
     * metadata for later persistent-index work.
     */
    std::string
        sequence_;

    std::uint64_t
        sequence_length_{0};

    PackedReference
        reference_;

    BidirectionalFMIndex
        index_;

    SensitivePairConstrainedSearchEngine
        pair_searcher_;
};

}  // namespace primerpair
