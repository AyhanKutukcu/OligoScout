#include "primerpair/sensitive_candidate_search.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>
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

        /*
         * Mutations intentionally cover:
         *
         * k=0
         * k=1
         * k=2
         * k=3 with 2+1 distribution
         * k=3 with 1+2 distribution
         * terminal 3-prime mismatch
         */
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

        const primerpair::
            SensitiveCandidateSearchEngine
                candidate_engine(
                    index,
                    packed
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

            const auto candidate_result =
                candidate_engine.search(
                    primer,
                    budget
                );

            expect(
                candidate_result.hits ==
                    reference_result.hits,
                "Candidate equals exhaustive reference k=" +
                    std::to_string(
                        budget
                    )
            );
        }

        std::cout
            << "All SENSITIVE candidate differential "
            << "tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
