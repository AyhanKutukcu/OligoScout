#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "primerpair/global_multiplex_cross_join.hpp"
#include "primerpair/multiplex_primer_search.hpp"
#include "primerpair/multiplex_primer_search_v2.hpp"
#include "primerpair/multiplex_strategy_router.hpp"


namespace primerpair {


/*
 * Unified production result for automatically routed
 * multiplex execution.
 *
 * Both V1 and V2 already share:
 *
 *   intended_pairs
 *   cross_amplicons
 *   MultiplexSearchStats
 *
 * V2 additionally exposes global sweep statistics.
 *
 * global_cross_stats_available is therefore false when
 * Production V1 executed.
 */
struct AdaptiveMultiplexPrimerSearchResult {

    MultiplexStrategyDecision decision{};


    std::vector<
        PrimerPairSearchResult
    > intended_pairs{};


    std::vector<
        MultiplexCrossAmplicon
    > cross_amplicons{};


    MultiplexSearchStats stats{};


    GlobalMultiplexCrossJoinStats
        global_cross_stats{};


    bool global_cross_stats_available{
        false
    };


    [[nodiscard]]
    bool used_production_v1() const noexcept {
        return
            decision.strategy ==
            MultiplexExecutionStrategy::
                ProductionV1;
    }


    [[nodiscard]]
    bool used_global_sweep_v2() const noexcept {
        return
            decision.strategy ==
            MultiplexExecutionStrategy::
                GlobalSweepV2;
    }
};


/*
 * Production adaptive multiplex execution layer.
 *
 * Routing decision:
 *
 *   small panel
 *       -> MultiplexPrimerSearchEngine
 *          Production V1
 *
 *   large panel
 *       -> MultiplexPrimerSearchEngineV2
 *          Global Sweep V2
 *
 * Default crossover is supplied by
 * MultiplexStrategyRouter.
 *
 * Current validated production default:
 *
 *   pair_count < 12  -> V1
 *   pair_count >= 12 -> V2
 *
 * The adaptive engine does not alter biological search
 * semantics. It only selects the cross-pair execution
 * backend.
 */
class AdaptiveMultiplexPrimerSearchEngine {
public:

    AdaptiveMultiplexPrimerSearchEngine(
        const MultiplexPrimerSearchEngine&
            production_v1,

        const MultiplexPrimerSearchEngineV2&
            global_sweep_v2,

        MultiplexStrategyRouter router =
            MultiplexStrategyRouter{}
    ) noexcept;


    [[nodiscard]]
    AdaptiveMultiplexPrimerSearchResult
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


    [[nodiscard]]
    const MultiplexStrategyRouter&
    router() const noexcept;


private:

    const MultiplexPrimerSearchEngine&
        production_v1_;


    const MultiplexPrimerSearchEngineV2&
        global_sweep_v2_;


    MultiplexStrategyRouter
        router_;
};


}  // namespace primerpair
