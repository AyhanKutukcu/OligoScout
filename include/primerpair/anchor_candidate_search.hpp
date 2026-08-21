#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"

namespace primerpair {

struct AnchorCandidateHit {
    std::uint64_t position{0};
    std::size_t mismatches{0};
};

struct AnchorCandidateSearchResult {
    std::size_t primer_length{0};
    std::size_t anchor_length{0};
    std::size_t max_mismatches{0};

    std::uint64_t anchor_occurrences{0};
    std::uint64_t candidates_verified{0};

    std::vector<AnchorCandidateHit> hits{};

    [[nodiscard]]
    bool empty() const noexcept {
        return hits.empty();
    }

    [[nodiscard]]
    std::size_t hit_count() const noexcept {
        return hits.size();
    }
};

class AnchorCandidateSearcher {
public:
    AnchorCandidateSearcher(
        const BidirectionalFMIndex& index,
        const PackedReference& reference
    ) noexcept;

    /*
     * STRICT profile:
     *
     * - 3' anchor must match exactly.
     * - mismatches are allowed only to the left
     *   of the exact anchor.
     * - Hamming distance only.
     * - max mismatch <= 3.
     */
    [[nodiscard]]
    AnchorCandidateSearchResult search(
        std::string_view primer,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 3
    ) const;


    /*
     * Candidate verification from an anchor hit
     * list produced by an external exact-search
     * backend such as IP-BWT.
     *
     * anchor_occurrences preserves the original
     * exact interval size before defensive
     * coordinate deduplication.
     */
    [[nodiscard]]
    AnchorCandidateSearchResult
    verify_from_anchor_positions(
        std::string_view primer,
        std::vector<std::uint64_t> anchor_positions,
        std::uint64_t anchor_occurrences,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 3
    ) const;

private:
    const BidirectionalFMIndex& index_;
    const PackedReference& reference_;
};

}  // namespace primerpair
