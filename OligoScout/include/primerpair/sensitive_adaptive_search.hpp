#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/sensitive_candidate_search.hpp"
#include "primerpair/sensitive_cost_estimator.hpp"
#include "primerpair/sensitive_primer_search.hpp"

namespace primerpair {

/*
 * Provisional v1 threshold derived from chr22
 * 2-fold cross-validation.
 *
 * This constant must remain explicit because it
 * still requires independent validation before
 * SearchProfile::Sensitive integration.
 */
inline constexpr std::uint64_t
    kDefaultSensitiveK3MaxSeedThreshold =
        12275ULL;

enum class SensitiveAdaptiveBackend {
    Exhaustive,
    Candidate
};

struct SensitiveAdaptiveSearchResult {
    SensitivePrimerSearchResult search_result{};

    SensitiveAdaptiveBackend backend{
        SensitiveAdaptiveBackend::Exhaustive
    };

    SensitiveCandidateCostEstimate
        cost_estimate{};

    bool estimator_used{false};

    [[nodiscard]]
    std::size_t hit_count() const noexcept {
        return search_result.hit_count();
    }

    [[nodiscard]]
    bool empty() const noexcept {
        return search_result.empty();
    }
};

class SensitiveAdaptiveSearchEngine {
public:
    SensitiveAdaptiveSearchEngine(
        const BidirectionalFMIndex& index,
        const PackedReference& reference,
        std::uint64_t k3_max_seed_threshold =
            kDefaultSensitiveK3MaxSeedThreshold
    );

    [[nodiscard]]
    SensitiveAdaptiveSearchResult search(
        std::string_view primer,
        std::size_t max_mismatches = 3
    ) const;

    [[nodiscard]]
    std::uint64_t
    k3_max_seed_threshold() const noexcept {
        return k3_max_seed_threshold_;
    }

private:
    SensitivePrimerSearchEngine
        exhaustive_;

    SensitiveCandidateSearchEngine
        candidate_;

    SensitiveCandidateCostEstimator
        estimator_;

    std::uint64_t
        k3_max_seed_threshold_;
};

}  // namespace primerpair
