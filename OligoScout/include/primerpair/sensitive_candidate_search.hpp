#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/sensitive_primer_search.hpp"

namespace primerpair {

/*
 * Lossless two-part candidate generator for SENSITIVE search.
 *
 * If the full-primer Hamming distance is <= k, then after
 * splitting the primer into two contiguous parts, at least
 * one part has Hamming distance <= floor(k / 2).
 *
 * Both parts are searched with that reduced budget.
 * Candidate full-primer starts are then verified exactly.
 */
class SensitiveCandidateSearchEngine {
    friend struct SensitiveCandidateSearchAccess;

    /*
     * Pair-constrained search reuses the validated
     * half-seed candidate generator and full-primer
     * verifier directly.
     */
    friend class SensitivePairConstrainedSearchEngine;

public:
    SensitiveCandidateSearchEngine(
        const BidirectionalFMIndex& index,
        const PackedReference& reference
    );

    [[nodiscard]]
    SensitivePrimerSearchResult search(
        std::string_view primer,
        std::size_t max_mismatches = 3
    ) const;

private:
    [[nodiscard]]
    std::vector<std::uint64_t>
    search_seed(
        std::string_view seed,
        std::size_t max_mismatches
    ) const;

    [[nodiscard]]
    std::vector<OrientedPrimerSearchHit>
    search_oriented(
        std::string_view query,
        PrimerOrientation orientation,
        bool reverse_to_original,
        std::size_t max_mismatches
    ) const;

    [[nodiscard]]
    bool verify_candidate(
        std::uint64_t start,
        std::string_view query,
        bool reverse_to_original,
        std::size_t max_mismatches,
        std::size_t& mismatches,
        std::uint64_t& mismatch_mask
    ) const;

    const BidirectionalFMIndex& index_;
    const PackedReference& reference_;
};

}  // namespace primerpair
