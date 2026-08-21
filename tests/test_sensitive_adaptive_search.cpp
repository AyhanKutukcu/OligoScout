#include "primerpair/sensitive_adaptive_search.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect(
    const bool condition,
    const std::string& name
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " + name
        );
    }

    std::cout
        << "[PASS] "
        << name
        << '\n';
}

char mutate_base(
    const char base
) {
    return
        base == 'A'
            ? 'C'
            : 'A';
}

std::string mutate(
    std::string sequence,
    const std::vector<std::size_t>& positions
) {
    for (
        const std::size_t position :
        positions
    ) {
        sequence.at(position) =
            mutate_base(
                sequence.at(position)
            );
    }

    return sequence;
}

}  // namespace

int main() {
    try {
        const std::string primer =
            "ACGTTGCAAGTCCTGAACGA";

        const std::vector<
            std::vector<std::size_t>
        > mutation_sets{
            {},
            {0},
            {19},
            {2, 12},
            {1, 3, 15},
            {2, 12, 16}
        };

        std::string reference(
            20,
            'N'
        );

        for (
            const auto& positions :
            mutation_sets
        ) {
            const std::string mutated =
                mutate(
                    primer,
                    positions
                );

            reference += mutated;

            reference.append(
                20,
                'N'
            );

            reference +=
                primerpair::
                    reverse_complement(
                        mutated
                    );

            reference.append(
                20,
                'N'
            );
        }

        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::
            SensitivePrimerSearchEngine
                reference_engine(
                    index
                );

        /*
         * UINT64_MAX forces the candidate route
         * for k=3.
         */
        const primerpair::
            SensitiveAdaptiveSearchEngine
                force_candidate(
                    index,
                    packed,
                    std::numeric_limits<
                        std::uint64_t
                    >::max()
                );

        for (
            std::size_t budget = 0;
            budget <= 3;
            ++budget
        ) {
            const auto reference_result =
                reference_engine.search(
                    primer,
                    budget
                );

            const auto adaptive_result =
                force_candidate.search(
                    primer,
                    budget
                );

            expect(
                adaptive_result
                    .search_result
                    .hits
                ==
                reference_result.hits,
                "Adaptive equals exhaustive reference k=" +
                    std::to_string(
                        budget
                    )
            );

            if (budget < 3) {

                expect(
                    adaptive_result.backend ==
                        primerpair::
                            SensitiveAdaptiveBackend::
                                Exhaustive,
                    "k<3 uses exhaustive backend k=" +
                        std::to_string(
                            budget
                        )
                );

                expect(
                    !adaptive_result
                        .estimator_used,
                    "k<3 skips estimator k=" +
                        std::to_string(
                            budget
                        )
                );

            } else {

                expect(
                    adaptive_result.backend ==
                        primerpair::
                            SensitiveAdaptiveBackend::
                                Candidate,
                    "Forced k3 candidate route"
                );

                expect(
                    adaptive_result
                        .estimator_used,
                    "k3 uses estimator"
                );
            }
        }

        /*
         * Threshold zero forces exhaustive for this
         * primer because its seed occurrence count
         * is non-zero.
         */
        const primerpair::
            SensitiveAdaptiveSearchEngine
                force_exhaustive(
                    index,
                    packed,
                    0
                );

        const auto reference_k3 =
            reference_engine.search(
                primer,
                3
            );

        const auto exhaustive_k3 =
            force_exhaustive.search(
                primer,
                3
            );

        expect(
            exhaustive_k3
                .search_result
                .hits
            ==
            reference_k3.hits,
            "Forced exhaustive k3 equals reference"
        );

        expect(
            exhaustive_k3.backend ==
                primerpair::
                    SensitiveAdaptiveBackend::
                        Exhaustive,
            "Forced k3 exhaustive route"
        );

        expect(
            exhaustive_k3.estimator_used,
            "Forced exhaustive k3 still uses estimator"
        );

        expect(
            exhaustive_k3
                .cost_estimate
                .max_seed_occurrences
            >
            0,
            "Synthetic k3 estimate is non-zero"
        );

        /*
         * Default threshold must remain explicit
         * and inspectable.
         */
        const primerpair::
            SensitiveAdaptiveSearchEngine
                default_engine(
                    index,
                    packed
                );

        expect(
            default_engine
                .k3_max_seed_threshold()
            ==
            primerpair::
                kDefaultSensitiveK3MaxSeedThreshold,
            "Default adaptive threshold"
        );

        /*
         * Lowercase must preserve exact hit set.
         */
        std::string lowercase =
            primer;

        std::transform(
            lowercase.begin(),
            lowercase.end(),
            lowercase.begin(),
            [](
                const unsigned char base
            ) {
                return
                    static_cast<char>(
                        std::tolower(
                            base
                        )
                    );
            }
        );

        const auto lowercase_result =
            force_candidate.search(
                lowercase,
                3
            );

        expect(
            lowercase_result
                .search_result
                .hits
            ==
            reference_k3.hits,
            "Adaptive lowercase search"
        );

        bool excessive_budget_rejected =
            false;

        try {
            static_cast<void>(
                default_engine.search(
                    primer,
                    4
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            excessive_budget_rejected =
                true;
        }

        expect(
            excessive_budget_rejected,
            "Adaptive mismatch budget >3 rejected"
        );

        std::cout
            << "default_k3_threshold\t"
            << primerpair::
                kDefaultSensitiveK3MaxSeedThreshold
            << '\n';

        std::cout
            << "All sensitive adaptive search tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
