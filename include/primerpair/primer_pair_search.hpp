#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/search_profile.hpp"
#include "primerpair/sensitive_adaptive_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {

enum class PrimerIdentity {
    Primer1,
    Primer2
};

struct PrimerPairHit {
    /*
     * Primer identity on the LEFT side of the
     * amplicon. This hit is always FORWARD.
     */
    PrimerIdentity left_primer{
        PrimerIdentity::Primer1
    };

    /*
     * Primer identity on the RIGHT side of the
     * amplicon. This hit is always REVERSE.
     */
    PrimerIdentity right_primer{
        PrimerIdentity::Primer2
    };

    std::uint64_t left_position{0};
    std::uint64_t right_position{0};

    std::size_t left_mismatches{0};
    std::size_t right_mismatches{0};

    /*
     * Half-open genomic amplicon:
     *
     * [amplicon_start, amplicon_end_exclusive)
     */
    std::uint64_t amplicon_start{0};

    std::uint64_t
        amplicon_end_exclusive{0};

    std::uint64_t amplicon_length{0};

    /*
     * Mismatch positions in each ORIGINAL primer's
     * 5-prime -> 3-prime coordinate system.
     */
    std::uint64_t left_mismatch_mask{0};
    std::uint64_t right_mismatch_mask{0};

    [[nodiscard]]
    std::size_t total_mismatches()
        const noexcept {
        return
            left_mismatches +
            right_mismatches;
    }

    bool operator==(
        const PrimerPairHit&
    ) const = default;
};


struct PrimerPairSearchResult {
    std::size_t primer1_length{0};
    std::size_t primer2_length{0};

    std::size_t primer1_single_hit_count{0};
    std::size_t primer2_single_hit_count{0};

    std::uint64_t min_amplicon_length{0};
    std::uint64_t max_amplicon_length{0};

    std::vector<PrimerPairHit> amplicons{};

    [[nodiscard]]
    bool empty() const noexcept {
        return amplicons.empty();
    }

    [[nodiscard]]
    std::size_t amplicon_count()
        const noexcept {
        return amplicons.size();
    }
};


[[nodiscard]]
const char* to_string(
    PrimerIdentity identity
) noexcept;


class SensitivePairConstrainedSearchEngine;

class PrimerPairSearchEngine {
public:
    PrimerPairSearchEngine(
        const BidirectionalFMIndex& index,
        const PackedReference& reference,
        SearchDifficultyThresholds thresholds = {}
    );

    /*
     * Legacy STRICT API.
     *
     * Kept unchanged for backward compatibility.
     */
    [[nodiscard]]
    PrimerPairSearchResult search(
        std::string_view primer1,
        std::string_view primer2,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 3,
        std::uint64_t min_amplicon_length = 50,
        std::uint64_t max_amplicon_length = 3000
    ) const;

    /*
     * Profile-aware API.
     *
     * STRICT:
     *     exactly the same validated path as the
     *     legacy search() overload.
     *
     * SENSITIVE:
     *     validated full-primer adaptive search
     *     for each primer, followed by the SAME
     *     pair sweep-line and merge kernel.
     *
     * anchor_length is ignored only for SENSITIVE.
     */
    [[nodiscard]]
    PrimerPairSearchResult search(
        std::string_view primer1,
        std::string_view primer2,
        SearchProfile profile,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 3,
        std::uint64_t min_amplicon_length = 50,
        std::uint64_t max_amplicon_length = 3000
    ) const;

private:
    /*
     * Allows the pair-constrained SENSITIVE engine
     * to reuse the exact same validated pairing
     * kernel without duplicating PCR geometry.
     */
    friend class SensitivePairConstrainedSearchEngine;

    StrandAwarePrimerSearchEngine
        strand_searcher_;

    SensitiveAdaptiveSearchEngine
        sensitive_searcher_;

    static void append_pairs(
        const std::vector<
            OrientedPrimerSearchHit
        >& forward_hits,

        std::size_t forward_primer_length,
        PrimerIdentity forward_identity,

        const std::vector<
            OrientedPrimerSearchHit
        >& reverse_hits,

        std::size_t reverse_primer_length,
        PrimerIdentity reverse_identity,

        std::uint64_t min_amplicon_length,
        std::uint64_t max_amplicon_length,

        std::vector<PrimerPairHit>& output
    );

    /*
     * Shared PCR pairing kernel.
     *
     * STRICT and SENSITIVE differ only in how
     * single-primer hits are generated.
     */
    static PrimerPairSearchResult assemble_pairs(
        std::string_view primer1,

        const std::vector<
            OrientedPrimerSearchHit
        >& primer1_hits,

        std::string_view primer2,

        const std::vector<
            OrientedPrimerSearchHit
        >& primer2_hits,

        std::uint64_t min_amplicon_length,
        std::uint64_t max_amplicon_length
    );
};

}  // namespace primerpair
