#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/sensitive_adaptive_search.hpp"
#include "primerpair/sensitive_candidate_search.hpp"
#include "primerpair/sensitive_cost_estimator.hpp"

namespace primerpair {

/*
 * Prototype v1:
 *
 * - 20-nt primers
 * - k = 3
 * - one primer searched globally
 * - mate searched only inside PCR-compatible
 *   genomic windows
 *
 * Lossless target:
 *
 * pair-constrained SENSITIVE amplicons
 * ==
 * global-global SENSITIVE amplicons
 */
struct SensitivePairConstrainedTiming {
    double cost_estimator_us{0.0};
    double anchor_search_us{0.0};
    double window_build_merge_us{0.0};

    /*
     * Includes half-seed FM search, locate(),
     * candidate-start construction and dedup.
     */
    double mate_seed_generation_us{0.0};

    /*
     * Includes PCR-window intersection,
     * full-primer verification, hit materialization
     * and final hit normalization.
     */
    double mate_filter_verify_us{0.0};

    double pair_assembly_us{0.0};

    [[nodiscard]]
    double accounted_total_us() const noexcept {
        return
            cost_estimator_us +
            anchor_search_us +
            window_build_merge_us +
            mate_seed_generation_us +
            mate_filter_verify_us +
            pair_assembly_us;
    }
};


enum class SensitivePairAnchorPolicy {
    AdaptiveLowCost,
    ForcePrimer1,
    ForcePrimer2
};


struct SensitivePairConstrainedSearchResult {
    PrimerIdentity anchor_primer{
        PrimerIdentity::Primer1
    };

    std::size_t anchor_global_hit_count{0};

    /*
     * Mate hits found only inside PCR-compatible
     * windows. This is deliberately NOT the mate's
     * genome-wide single-primer hit count.
     */
    std::size_t mate_local_hit_count{0};

    std::size_t
        merged_forward_mate_windows{0};

    std::size_t
        merged_reverse_mate_windows{0};

    /*
     * V1: every genomic start scanned inside windows.
     *
     * V2: only half-seed candidate starts surviving
     * the PCR-window filter and reaching full-primer
     * verification.
     *
     * Field name retained for benchmark/API
     * continuity.
     */
    std::uint64_t
        scanned_mate_start_positions{0};

    SensitiveCandidateCostEstimate
        primer1_cost{};

    SensitiveCandidateCostEstimate
        primer2_cost{};

    SensitivePairConstrainedTiming
        timing{};

    PrimerPairSearchResult
        pair_result{};
};


class SensitivePairConstrainedSearchEngine {
public:
    SensitivePairConstrainedSearchEngine(
        const BidirectionalFMIndex& index,
        const PackedReference& reference
    );

    [[nodiscard]]
    SensitivePairConstrainedSearchResult search(
        std::string_view primer1,
        std::string_view primer2,
        std::size_t max_mismatches = 3,
        std::uint64_t min_amplicon_length = 50,
        std::uint64_t max_amplicon_length = 3000,
        SensitivePairAnchorPolicy anchor_policy =
            SensitivePairAnchorPolicy::AdaptiveLowCost
    ) const;

private:
    SensitiveAdaptiveSearchEngine
        anchor_searcher_;

    /*
     * Reuses the lossless k=3 half-seed candidate
     * generator for the locally constrained mate.
     */
    SensitiveCandidateSearchEngine
        candidate_searcher_;

    SensitiveCandidateCostEstimator
        cost_estimator_;

    const PackedReference&
        reference_;
};

}  // namespace primerpair
