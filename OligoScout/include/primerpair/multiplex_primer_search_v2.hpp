#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "primerpair/global_multiplex_cross_join.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/multiplex_primer_search.hpp"

namespace primerpair {


struct MultiplexPrimerSearchResultV2 {
    std::vector<PrimerPairSearchResult>
        intended_pairs{};

    std::vector<MultiplexCrossAmplicon>
        cross_amplicons{};

    MultiplexSearchStats stats{};

    GlobalMultiplexCrossJoinStats
        global_cross_stats{};
};


class MultiplexPrimerSearchEngineV2 {
public:
    explicit
    MultiplexPrimerSearchEngineV2(
        const HybridStrandAwarePrimerSearchEngine&
            primer_engine
    ) noexcept;


    [[nodiscard]]
    MultiplexPrimerSearchResultV2
    search(
        const std::vector<
            MultiplexPrimerPairRequest
        >& requests,

        std::size_t anchor_length = 12,

        bool include_cross_pairs = true,

        std::uint64_t
            cross_min_amplicon_length = 50,

        std::uint64_t
            cross_max_amplicon_length = 3000
    ) const;


private:
    const HybridStrandAwarePrimerSearchEngine&
        primer_engine_;
};


}  // namespace primerpair
