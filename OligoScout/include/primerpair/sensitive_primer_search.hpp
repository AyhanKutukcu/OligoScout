#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {

/*
 * ------------------------------------------------------------
 * SENSITIVE reference search
 *
 * Purpose:
 *
 * Exhaustive Hamming-distance candidate generation across
 * the ENTIRE primer.
 *
 * No position, including the biological 3-prime end,
 * is forced to be exact.
 *
 * This is initially a correctness/reference backend.
 * Performance optimization comes later.
 * ------------------------------------------------------------
 */

struct SensitivePrimerSearchResult {
    std::size_t primer_length{0};

    std::size_t max_mismatches{0};

    std::vector<
        OrientedPrimerSearchHit
    > hits{};

    [[nodiscard]]
    bool empty() const noexcept {
        return hits.empty();
    }

    [[nodiscard]]
    std::size_t hit_count() const noexcept {
        return hits.size();
    }
};

class SensitivePrimerSearchEngine {
public:
    explicit SensitivePrimerSearchEngine(
        const BidirectionalFMIndex& index
    );

    [[nodiscard]]
    SensitivePrimerSearchResult search(
        std::string_view primer,
        std::size_t max_mismatches = 3
    ) const;

private:
    struct BranchState {
        BidirectionalInterval interval{};

        std::size_t mismatches{0};

        /*
         * Coordinates are always expressed in the
         * ORIGINAL primer's 5-prime -> 3-prime system.
         */
        std::uint64_t mismatch_mask{0};
    };

    [[nodiscard]]
    std::vector<OrientedPrimerSearchHit>
    search_reference_oriented_query(
        std::string_view query,
        PrimerOrientation orientation,
        bool reverse_to_original,
        std::size_t max_mismatches
    ) const;

    const BidirectionalFMIndex& index_;
};

}  // namespace primerpair
