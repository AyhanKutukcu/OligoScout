#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/primer_pair_search.hpp"

namespace primerpair {


struct HybridPrimerPairRequest {
    std::string_view primer1{};
    std::string_view primer2{};

    std::size_t max_mismatches{3};

    std::uint64_t min_amplicon_length{50};
    std::uint64_t max_amplicon_length{3000};
};


struct HybridPrimerPairResult {
    HybridStrandAwarePrimerResult
        primer1_search{};

    HybridStrandAwarePrimerResult
        primer2_search{};

    PrimerPairSearchResult
        pair_result{};
};


class HybridPrimerPairSearchEngine {
public:
    explicit
    HybridPrimerPairSearchEngine(
        const HybridStrandAwarePrimerSearchEngine&
            primer_engine
    ) noexcept;

    /*
     * Tüm pair üyeleri tek primer batch'e flatten edilir:
     *
     * P1a, P1b, P2a, P2b, ...
     *
     * Ardından precomputed hit listeleri reusable
     * O(F + R + K) primer-pair assembler'a verilir.
     */
    [[nodiscard]]
    std::vector<HybridPrimerPairResult>
    search(
        const std::vector<
            HybridPrimerPairRequest
        >& requests,
        std::size_t anchor_length = 12
    ) const;

private:
    const HybridStrandAwarePrimerSearchEngine&
        primer_engine_;
};


}  // namespace primerpair
