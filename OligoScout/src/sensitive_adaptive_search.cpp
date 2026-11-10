#include "primerpair/sensitive_adaptive_search.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string_view>

namespace primerpair {

SensitiveAdaptiveSearchEngine::
SensitiveAdaptiveSearchEngine(
    const BidirectionalFMIndex& index,
    const PackedReference& reference,
    const std::uint64_t
        k3_max_seed_threshold
)
    : exhaustive_(
          index
      ),
      candidate_(
          index,
          reference
      ),
      estimator_(
          index
      ),
      k3_max_seed_threshold_(
          k3_max_seed_threshold
      ) {
}

SensitiveAdaptiveSearchResult
SensitiveAdaptiveSearchEngine::search(
    const std::string_view primer,
    const std::size_t max_mismatches
) const {
    if (max_mismatches > 3) {
        throw std::invalid_argument(
            "Adaptive SENSITIVE MVP supports "
            "at most 3 mismatches."
        );
    }

    SensitiveAdaptiveSearchResult
        result;

    /*
     * k=0,1,2:
     *
     * chr22 benchmarks showed exhaustive search
     * is faster than the half-seed candidate
     * backend for these mismatch budgets.
     *
     * Do not pay estimator overhead here.
     */
    if (max_mismatches < 3) {

        result.search_result =
            exhaustive_.search(
                primer,
                max_mismatches
            );

        result.backend =
            SensitiveAdaptiveBackend::
                Exhaustive;

        result.estimator_used =
            false;

        return result;
    }

    /*
     * k=3 outside the OligoScout MVP
     * primer-length range:
     *
     * preserve a conservative correctness-first
     * fallback instead of applying a router that
     * has not been benchmarked for these lengths.
     */
    if (
        primer.size() < 18
        ||
        primer.size() > 35
    ) {
        result.search_result =
            exhaustive_.search(
                primer,
                max_mismatches
            );

        result.backend =
            SensitiveAdaptiveBackend::
                Exhaustive;

        result.estimator_used =
            false;

        return result;
    }

    /*
     * k=3, primer length 30..35:
     *
     * GRCh38 chr22 SA8 differential benchmarking
     * showed the lossless half-seed Candidate
     * backend overwhelmingly dominates exhaustive
     * search at these lengths.
     *
     * The estimator itself becomes measurable
     * overhead once half-seeds are this specific,
     * so bypass it entirely.
     */
    if (primer.size() >= 30) {

        result.search_result =
            candidate_.search(
                primer,
                max_mismatches
            );

        result.backend =
            SensitiveAdaptiveBackend::
                Candidate;

        result.estimator_used =
            false;

        return result;
    }

    /*
     * k=3, primer length 18..29:
     *
     * This is the occurrence-sensitive region.
     * Estimate repeat/candidate pressure using
     * the maximum <=1-mismatch half-seed
     * occurrence count.
     *
     * The existing threshold 12275 is retained:
     * a learned universal short-primer threshold
     * improved mean CV runtime by only ~0.6%,
     * which does not justify extra policy
     * complexity at this stage.
     */
    result.cost_estimate =
        estimator_.estimate_k3(
            primer
        );

    result.estimator_used =
        true;

    if (
        result.cost_estimate
            .max_seed_occurrences
        <=
        k3_max_seed_threshold_
    ) {
        result.search_result =
            candidate_.search(
                primer,
                max_mismatches
            );

        result.backend =
            SensitiveAdaptiveBackend::
                Candidate;

    } else {

        result.search_result =
            exhaustive_.search(
                primer,
                max_mismatches
            );

        result.backend =
            SensitiveAdaptiveBackend::
                Exhaustive;
    }

    return result;
}

}  // namespace primerpair
