#include "primerpair/adaptive_multiplex_primer_search.hpp"

#include <utility>


namespace primerpair {


AdaptiveMultiplexPrimerSearchEngine::
AdaptiveMultiplexPrimerSearchEngine(
    const MultiplexPrimerSearchEngine&
        production_v1,

    const MultiplexPrimerSearchEngineV2&
        global_sweep_v2,

    MultiplexStrategyRouter router
) noexcept
    :
    production_v1_(
        production_v1
    ),
    global_sweep_v2_(
        global_sweep_v2
    ),
    router_(
        std::move(
            router
        )
    ) {
}


AdaptiveMultiplexPrimerSearchResult
AdaptiveMultiplexPrimerSearchEngine::search(
    const std::vector<
        MultiplexPrimerPairRequest
    >& requests,

    const std::size_t anchor_length,

    const bool include_cross_pairs,

    const std::uint64_t
        cross_min_amplicon_length,

    const std::uint64_t
        cross_max_amplicon_length
) const {

    const MultiplexStrategyDecision
        decision =
            router_.decide(
                requests.size()
            );


    AdaptiveMultiplexPrimerSearchResult
        output;


    output.decision =
        decision;


    if (
        decision.strategy ==
        MultiplexExecutionStrategy::
            ProductionV1
    ) {
        auto selected =
            production_v1_.search(
                requests,
                anchor_length,
                include_cross_pairs,
                cross_min_amplicon_length,
                cross_max_amplicon_length
            );


        output.intended_pairs =
            std::move(
                selected.intended_pairs
            );


        output.cross_amplicons =
            std::move(
                selected.cross_amplicons
            );


        output.stats =
            selected.stats;


        output.global_cross_stats_available =
            false;


        return output;
    }


    auto selected =
        global_sweep_v2_.search(
            requests,
            anchor_length,
            include_cross_pairs,
            cross_min_amplicon_length,
            cross_max_amplicon_length
        );


    output.intended_pairs =
        std::move(
            selected.intended_pairs
        );


    output.cross_amplicons =
        std::move(
            selected.cross_amplicons
        );


    output.stats =
        selected.stats;


    output.global_cross_stats =
        selected.global_cross_stats;


    /*
     * When cross-pair analysis is disabled the V2 engine
     * legitimately returns default/empty global sweep
     * statistics.
     */
    output.global_cross_stats_available =
        include_cross_pairs;


    return output;
}


const MultiplexStrategyRouter&
AdaptiveMultiplexPrimerSearchEngine::
router() const noexcept {

    return router_;
}


}  // namespace primerpair
