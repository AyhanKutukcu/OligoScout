#include "primerpair/multiplex_strategy_router.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>


namespace {


void expect(
    const bool condition,
    const char* message
) {
    if (!condition) {
        throw std::runtime_error(
            message
        );
    }


    std::cout
        << "[PASS] "
        << message
        << '\n';
}


}  // namespace


int main() {
    using namespace primerpair;


    try {

        const MultiplexStrategyRouter
            router;


        const std::size_t threshold =
            router
                .crossover_pair_count();


        expect(
            threshold ==
                kDefaultMultiplexV2CrossoverPairCount,
            "Default router uses measured crossover"
        );


        expect(
            threshold >= 2,
            "Measured crossover is valid"
        );


        expect(
            router
                .decide(0)
                .strategy
            ==
            MultiplexExecutionStrategy::
                ProductionV1,
            "Zero pairs route to Production V1"
        );


        expect(
            router
                .decide(1)
                .strategy
            ==
            MultiplexExecutionStrategy::
                ProductionV1,
            "One pair routes to Production V1"
        );


        expect(
            router
                .decide(
                    threshold - 1
                )
                .strategy
            ==
            MultiplexExecutionStrategy::
                ProductionV1,
            "Below crossover routes to Production V1"
        );


        expect(
            router
                .decide(
                    threshold
                )
                .strategy
            ==
            MultiplexExecutionStrategy::
                GlobalSweepV2,
            "Crossover boundary routes to Global Sweep V2"
        );


        expect(
            router
                .decide(
                    threshold + 1
                )
                .strategy
            ==
            MultiplexExecutionStrategy::
                GlobalSweepV2,
            "Above crossover routes to Global Sweep V2"
        );


        const MultiplexStrategyRouter
            custom_router(
                100
            );


        expect(
            custom_router
                .decide(99)
                .strategy
            ==
            MultiplexExecutionStrategy::
                ProductionV1,
            "Custom threshold below boundary routes to V1"
        );


        expect(
            custom_router
                .decide(100)
                .strategy
            ==
            MultiplexExecutionStrategy::
                GlobalSweepV2,
            "Custom threshold boundary routes to V2"
        );


        bool rejected = false;


        try {

            const MultiplexStrategyRouter
                invalid(
                    1
                );

            (void)invalid;

        } catch (
            const std::invalid_argument&
        ) {

            rejected = true;
        }


        expect(
            rejected,
            "Invalid threshold is rejected"
        );


        std::cout
            << "measured_threshold\t"
            << threshold
            << '\n';


        std::cout
            << "MULTIPLEX_ROUTER_BOUNDARY\tYES\n";


        std::cout
            << "MULTIPLEX_ROUTER_OVERRIDE\tYES\n";


        std::cout
            << "ALL_CHECKS\tYES\n";


        return 0;

    } catch (
        const std::exception& error
    ) {

        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';


        return 1;
    }
}
