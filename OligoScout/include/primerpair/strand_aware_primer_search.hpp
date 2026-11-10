#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/packed_reference.hpp"
#include "primerpair/search_strategy.hpp"
#include "primerpair/single_primer_search.hpp"

namespace primerpair {

/*
 * Orientation is defined relative to the reference sequence.
 *
 * Forward:
 *   input primer sequence matches reference left-to-right.
 *   Polymerase extension proceeds toward increasing coordinates.
 *
 * Reverse:
 *   reverse-complement of the input primer matches reference.
 *   Polymerase extension proceeds toward decreasing coordinates.
 */
enum class PrimerOrientation {
    Forward,
    Reverse
};

struct OrientedPrimerSearchHit {
    std::uint64_t position{0};

    std::size_t mismatches{0};

    PrimerOrientation orientation{
        PrimerOrientation::Forward
    };

    /*
     * Mismatch positions expressed in the ORIGINAL
     * primer's 5-prime -> 3-prime coordinates.
     *
     * bit 0 = primer position 0 (5-prime-most base)
     *
     * This convention is identical for Forward and
     * Reverse genomic orientations.
     *
     * OligoScout MVP primers are 18-35 nt,
     * therefore uint64_t is sufficient.
     */
    std::uint64_t mismatch_mask{0};

    bool operator==(
        const OrientedPrimerSearchHit&
    ) const = default;
};

struct StrandAwarePrimerSearchResult {
    std::size_t primer_length{0};

    SearchStrategyDecision forward_decision{};
    SearchStrategyDecision reverse_decision{};

    std::vector<OrientedPrimerSearchHit> hits{};

    [[nodiscard]]
    bool empty() const noexcept {
        return hits.empty();
    }

    [[nodiscard]]
    std::size_t hit_count() const noexcept {
        return hits.size();
    }
};

[[nodiscard]]
std::string reverse_complement(
    std::string_view sequence
);

[[nodiscard]]
const char* to_string(
    PrimerOrientation orientation
) noexcept;

class StrandAwarePrimerSearchEngine {
public:
    StrandAwarePrimerSearchEngine(
        const BidirectionalFMIndex& index,
        const PackedReference& reference,
        SearchDifficultyThresholds thresholds = {}
    );

    [[nodiscard]]
    StrandAwarePrimerSearchResult search(
        std::string_view primer,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 3
    ) const;

private:
    const BidirectionalFMIndex& index_;
    const PackedReference& reference_;

    SearchStrategyRouter router_;

    /*
     * Existing validated engine handles the
     * Forward orientation.
     */
    SinglePrimerSearchEngine forward_engine_;

    [[nodiscard]]
    SearchStrategyDecision
    make_reverse_decision(
        std::string_view reverse_query,
        std::size_t anchor_length,
        std::size_t max_mismatches
    ) const;

    [[nodiscard]]
    std::vector<OrientedPrimerSearchHit>
    execute_reverse_branching(
        std::string_view reverse_query,
        std::size_t anchor_length,
        std::size_t max_mismatches,
        const BidirectionalInterval& anchor_state
    ) const;

    [[nodiscard]]
    std::vector<OrientedPrimerSearchHit>
    execute_reverse_candidate(
        std::string_view reverse_query,
        std::size_t anchor_length,
        std::size_t max_mismatches,
        const BidirectionalInterval& anchor_state
    ) const;
};

}  // namespace primerpair
