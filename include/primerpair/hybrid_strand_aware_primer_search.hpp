#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_strategy.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {

struct HybridStrandAwarePrimerRequest {
    std::string_view primer{};
    std::size_t max_mismatches{3};
};


struct HybridStrandAwarePrimerResult {
    std::size_t primer_length{0};

    SearchStrategy forward_strategy{
        SearchStrategy::DirectBranching
    };

    SearchStrategy reverse_strategy{
        SearchStrategy::DirectBranching
    };

    bool forward_candidate_backend{false};

    bool reverse_candidate_backend{false};

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


/*
 * Production both-strand hybrid layer.
 *
 * Forward:
 *   existing HybridBatchedPrimerSearchEngine
 *
 * Reverse:
 *   reverse_complement(original)
 *
 *   biological 3-prime anchor =
 *   PREFIX of reverse_query.
 *
 * Candidate route:
 *   Prefix IP-BWT batch anchor
 *       +
 *   PackedReference verification
 *
 * DirectBranching route:
 *   validated StrandAwarePrimerSearchEngine
 *   fallback.
 *
 * Existing engines are deliberately not modified.
 */
class HybridStrandAwarePrimerSearchEngine {
public:
    HybridStrandAwarePrimerSearchEngine(
        const HybridBatchedPrimerSearchEngine&
            forward_hybrid,
        const BatchedAnchorLookup&
            anchor_lookup,
        const IPBWTIndex&
            ipbwt,
        const PackedReference&
            reference,
        const StrandAwarePrimerSearchEngine&
            legacy_strand_engine
    ) noexcept;

    [[nodiscard]]
    std::vector<
        HybridStrandAwarePrimerResult
    >
    search(
        const std::vector<
            HybridStrandAwarePrimerRequest
        >& requests,
        std::size_t anchor_length = 12
    ) const;

private:
    const HybridBatchedPrimerSearchEngine&
        forward_hybrid_;

    const BatchedAnchorLookup&
        anchor_lookup_;

    const IPBWTIndex&
        ipbwt_;

    const PackedReference&
        reference_;

    const StrandAwarePrimerSearchEngine&
        legacy_strand_engine_;

    [[nodiscard]]
    std::vector<
        OrientedPrimerSearchHit
    >
    verify_reverse_candidates(
        std::string_view reverse_query,
        std::vector<std::uint64_t> positions,
        std::size_t anchor_length,
        std::size_t max_mismatches
    ) const;
};

}  // namespace primerpair
