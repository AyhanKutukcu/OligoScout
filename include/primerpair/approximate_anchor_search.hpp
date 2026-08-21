#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"

namespace primerpair {

struct ApproximateAnchorHit {
    BidirectionalInterval state{};

    std::size_t mismatches{0};

    /*
     * Primer üzerindeki 0-based pozisyonlar.
     *
     * 0 = primerin 5' ucundaki ilk baz.
     *
     * Bu ilk sürümde mismatch'ler yalnızca
     * exact 3' anchor'ın solunda olabilir.
     */
    std::vector<std::size_t>
        mismatch_positions{};

    [[nodiscard]]
    bool empty() const noexcept {
        return state.empty();
    }

    [[nodiscard]]
    std::uint64_t match_count() const noexcept {
        return state.size();
    }
};

struct ApproximateAnchorSearchResult {
    std::size_t primer_length{0};
    std::size_t anchor_length{0};
    std::size_t max_mismatches{0};

    std::vector<ApproximateAnchorHit>
        hits{};

    [[nodiscard]]
    bool empty() const noexcept {
        return hits.empty();
    }

    [[nodiscard]]
    std::uint64_t total_match_count() const noexcept {
        std::uint64_t total = 0;

        for (const auto& hit : hits) {
            total += hit.match_count();
        }

        return total;
    }
};

class ApproximateAnchorSearcher {
public:
    explicit ApproximateAnchorSearcher(
        const BidirectionalFMIndex& index
    ) noexcept;

    /*
     * Exact 3' anchor +
     * 5' tarafta Hamming mismatch branching.
     *
     * max_mismatches:
     * 0..3
     */
    [[nodiscard]]
    ApproximateAnchorSearchResult
    search_5prime_mismatches(
        std::string_view primer,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 1
    ) const;

    [[nodiscard]]
    std::vector<std::uint64_t> locate(
        const ApproximateAnchorHit& hit
    ) const;

private:
    const BidirectionalFMIndex& index_;
};

}  // namespace primerpair
