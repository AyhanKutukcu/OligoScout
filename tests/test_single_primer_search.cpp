#include "primerpair/single_primer_search.hpp"
#include "primerpair/approximate_anchor_search.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
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

std::vector<primerpair::PrimerSearchHit>
branching_baseline(
    const primerpair::BidirectionalFMIndex& index,
    const primerpair::ApproximateAnchorSearcher& searcher,
    const std::string& primer,
    const std::size_t anchor_length,
    const std::size_t budget
) {
    const auto result =
        searcher.search_5prime_mismatches(
            primer,
            anchor_length,
            budget
        );

    std::vector<
        primerpair::PrimerSearchHit
    > hits;

    for (const auto& branch : result.hits) {

        const auto positions =
            index.locate(
                branch.state
            );

        for (const auto position : positions) {

            hits.push_back(
                {
                    position,
                    branch.mismatches
                }
            );
        }
    }

    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const auto& lhs,
            const auto& rhs
        ) {
            if (
                lhs.position !=
                rhs.position
            ) {
                return
                    lhs.position <
                    rhs.position;
            }

            return
                lhs.mismatches <
                rhs.mismatches;
        }
    );

    hits.erase(
        std::unique(
            hits.begin(),
            hits.end()
        ),
        hits.end()
    );

    return hits;
}

}  // namespace

int main() {
    try {
        /*
         * --------------------------------------------------
         * Controlled reference
         * --------------------------------------------------
         */

        const std::string primer =
            "ACGTTGCAACGTACGT";

        const std::string reference =
            "TTT"
            "ACGTTGCAACGTACGT"
            "GGG"
            "CCGTTGCAACGTACGT"
            "CCC"
            "TAGTTGCAACGTACGT"
            "GGG"
            "ACGTTGCAACGTACGA"
            "AAA";

        const primerpair::PackedReference
            packed_reference(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::ApproximateAnchorSearcher
            baseline_searcher(
                index
            );

        const primerpair::SinglePrimerSearchEngine
            engine(
                index,
                packed_reference
            );

        /*
         * k0 should route through direct branching
         * for this EASY anchor.
         */
        const auto k0 =
            engine.search(
                primer,
                12,
                0
            );

        expect(
            k0.decision.recommended_strategy ==
                primerpair::
                    SearchStrategy::
                    DirectBranching,
            "Controlled k0 uses direct backend"
        );

        /*
         * k3 should route through candidate
         * verification for EASY.
         */
        const auto k3 =
            engine.search(
                primer,
                12,
                3
            );

        expect(
            k3.decision.recommended_strategy ==
                primerpair::
                    SearchStrategy::
                    AnchorCandidateVerification,
            "Controlled k3 uses candidate backend"
        );

        for (std::size_t budget = 0;
             budget <= 3;
             ++budget) {

            const auto routed =
                engine.search(
                    primer,
                    12,
                    budget
                );

            const auto baseline =
                branching_baseline(
                    index,
                    baseline_searcher,
                    primer,
                    12,
                    budget
                );

            expect(
                routed.hits ==
                    baseline,
                "Controlled routed equivalence k=" +
                    std::to_string(
                        budget
                    )
            );
        }

        expect(
            k3.hits ==
                std::vector<
                    primerpair::PrimerSearchHit
                >{
                    {3, 0},
                    {22, 1},
                    {41, 2}
                },
            "Controlled final hit list"
        );

        /*
         * --------------------------------------------------
         * Deterministic randomized differential test
         * --------------------------------------------------
         */

        std::string random_reference(
            4096,
            'A'
        );

        constexpr char dna[] = {
            'A',
            'C',
            'G',
            'T'
        };

        std::uint64_t rng =
            0xA0761D6478BD642FULL;

        const auto next_random =
            [&rng]() -> std::uint64_t {

                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;

                return rng;
            };

        for (std::size_t i = 0;
             i < random_reference.size();
             ++i) {

            random_reference.at(i) =
                dna[
                    next_random() & 3ULL
                ];
        }

        const primerpair::PackedReference
            random_packed(
                random_reference
            );

        const primerpair::BidirectionalFMIndex
            random_index(
                random_reference
            );

        const primerpair::ApproximateAnchorSearcher
            random_baseline(
                random_index
            );

        const primerpair::SinglePrimerSearchEngine
            random_engine(
                random_index,
                random_packed
            );

        constexpr std::size_t
            primer_length = 20;

        constexpr std::size_t
            anchor_length = 12;

        bool randomized_ok = true;

        for (std::size_t trial = 0;
             trial < 250;
             ++trial) {

            const std::size_t available =
                random_reference.size() -
                primer_length +
                1;

            const std::size_t position =
                static_cast<std::size_t>(
                    next_random() %
                    available
                );

            const std::string current_primer =
                random_reference.substr(
                    position,
                    primer_length
                );

            for (std::size_t budget = 0;
                 budget <= 3;
                 ++budget) {

                const auto routed =
                    random_engine.search(
                        current_primer,
                        anchor_length,
                        budget
                    );

                const auto baseline =
                    branching_baseline(
                        random_index,
                        random_baseline,
                        current_primer,
                        anchor_length,
                        budget
                    );

                if (
                    routed.hits !=
                    baseline
                ) {
                    randomized_ok =
                        false;

                    break;
                }
            }

            if (!randomized_ok) {
                break;
            }
        }

        expect(
            randomized_ok,
            "250 randomized routed differential trials"
        );

        /*
         * Lowercase normalization.
         */
        const auto lowercase =
            engine.search(
                "acgttgcaacgtacgt",
                12,
                3
            );

        expect(
            lowercase.hits ==
                k3.hits,
            "Lowercase routed search"
        );

        bool excessive_budget_rejected =
            false;

        try {
            static_cast<void>(
                engine.search(
                    primer,
                    12,
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
            "Dispatcher mismatch >3 rejection"
        );

        std::cout
            << "All single-primer dispatcher tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
