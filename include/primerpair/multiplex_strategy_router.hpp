#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>


namespace primerpair {


enum class MultiplexExecutionStrategy :
    std::uint8_t {

    ProductionV1 = 0,

    GlobalSweepV2 = 1
};


struct MultiplexStrategyDecision {

    MultiplexExecutionStrategy strategy{
        MultiplexExecutionStrategy::
            ProductionV1
    };


    std::size_t pair_count{0};


    std::size_t crossover_pair_count{
        12
    };


    bool operator==(
        const MultiplexStrategyDecision&
    ) const = default;
};


/*
 * Empirical production default.
 *
 * Derived from:
 *
 *   GRCh38 chr22
 *   benchmark_grch38_shard_multiplex_v2
 *   panel sizes 8/12/16/20/24/32
 *   3 outer runs per panel size
 *   10 timing repetitions per outer run
 *
 * Selection rule:
 *
 *   median>=1.05_and_all_larger>=1.00
 *
 * This is NOT a biological constant.
 *
 * Re-benchmark after material hardware/compiler/backend
 * changes.
 */

/*
 * Empirical crossover interpretation:
 *
 * Observed median crossover: 8 primer pairs.
 *
 * At 8 pairs one independent timing run still favored
 * Production V1, so 8 is considered a marginal/noisy
 * crossover rather than the production boundary.
 *
 * Robust production default: 12 primer pairs.
 *
 * From 12 through 32 pairs every measured outer run
 * favored Global Sweep V2.
 */

inline constexpr std::size_t
kDefaultMultiplexV2CrossoverPairCount =
    12;


class MultiplexStrategyRouter {
public:

    explicit
    MultiplexStrategyRouter(
        const std::size_t
            crossover_pair_count =
                kDefaultMultiplexV2CrossoverPairCount
    )
        :
        crossover_pair_count_(
            crossover_pair_count
        )
    {
        if (
            crossover_pair_count_ <
            2
        ) {
            throw std::invalid_argument(
                "Multiplex crossover pair count "
                "must be >= 2."
            );
        }
    }


    [[nodiscard]]
    MultiplexStrategyDecision
    decide(
        const std::size_t pair_count
    ) const noexcept {

        const auto strategy =
            (
                pair_count >=
                    crossover_pair_count_
                &&
                pair_count >= 2
            )
            ?
            MultiplexExecutionStrategy::
                GlobalSweepV2
            :
            MultiplexExecutionStrategy::
                ProductionV1;


        return MultiplexStrategyDecision{
            strategy,
            pair_count,
            crossover_pair_count_
        };
    }


    [[nodiscard]]
    std::size_t
    crossover_pair_count() const noexcept {

        return crossover_pair_count_;
    }


private:

    std::size_t
        crossover_pair_count_;
};


}  // namespace primerpair
