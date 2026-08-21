#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/approximate_anchor_search.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
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

using PositionMismatch =
    std::pair<
        std::uint64_t,
        std::size_t
    >;

std::vector<PositionMismatch>
collect_branching(
    const primerpair::BidirectionalFMIndex& index,
    const primerpair::ApproximateAnchorSearchResult& result
) {
    std::vector<PositionMismatch> output;

    for (const auto& hit : result.hits) {

        const auto positions =
            index.locate(
                hit.state
            );

        for (const auto position : positions) {
            output.emplace_back(
                position,
                hit.mismatches
            );
        }
    }

    std::sort(
        output.begin(),
        output.end()
    );

    output.erase(
        std::unique(
            output.begin(),
            output.end()
        ),
        output.end()
    );

    return output;
}

std::vector<PositionMismatch>
collect_candidates(
    const primerpair::AnchorCandidateSearchResult& result
) {
    std::vector<PositionMismatch> output;

    for (const auto& hit : result.hits) {
        output.emplace_back(
            hit.position,
            hit.mismatches
        );
    }

    std::sort(
        output.begin(),
        output.end()
    );

    return output;
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

        const std::string exact =
            "ACGTTGCAACGTACGT";

        const std::string one_mismatch =
            "CCGTTGCAACGTACGT";

        const std::string two_mismatch =
            "TAGTTGCAACGTACGT";

        /*
         * Last anchor base differs, therefore STRICT
         * backend must never return this sequence.
         */
        const std::string anchor_mismatch =
            "ACGTTGCAACGTACGA";

        const std::string reference_text =
            "TTT" +
            exact +
            "GGG" +
            one_mismatch +
            "CCC" +
            two_mismatch +
            "GGG" +
            anchor_mismatch +
            "AAA";

        const primerpair::BidirectionalFMIndex
            index(
                reference_text
            );

        const primerpair::PackedReference
            packed_reference(
                reference_text
            );

        const primerpair::ApproximateAnchorSearcher
            branching(
                index
            );

        const primerpair::AnchorCandidateSearcher
            candidates(
                index,
                packed_reference
            );

        for (std::size_t budget = 0;
             budget <= 3;
             ++budget) {

            const auto branching_result =
                branching
                    .search_5prime_mismatches(
                        primer,
                        12,
                        budget
                    );

            const auto candidate_result =
                candidates.search(
                    primer,
                    12,
                    budget
                );

            expect(
                collect_candidates(
                    candidate_result
                ) ==
                collect_branching(
                    index,
                    branching_result
                ),
                "Controlled differential k=" +
                    std::to_string(budget)
            );
        }

        const auto three =
            candidates.search(
                primer,
                12,
                3
            );

        const auto controlled =
            collect_candidates(
                three
            );

        expect(
            controlled ==
                std::vector<PositionMismatch>{
                    {3, 0},
                    {22, 1},
                    {41, 2}
                },
            "Controlled candidate positions/mismatches"
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
            0xD1B54A32D192ED03ULL;

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

        const primerpair::BidirectionalFMIndex
            random_index(
                random_reference
            );

        const primerpair::PackedReference
            random_packed(
                random_reference
            );

        const primerpair::ApproximateAnchorSearcher
            random_branching(
                random_index
            );

        const primerpair::AnchorCandidateSearcher
            random_candidates(
                random_index,
                random_packed
            );

        bool randomized_ok = true;

        constexpr std::size_t
            random_primer_length = 20;

        constexpr std::size_t
            random_anchor_length = 12;

        for (std::size_t trial = 0;
             trial < 250;
             ++trial) {

            const std::size_t available =
                random_reference.size() -
                random_primer_length +
                1;

            const std::size_t position =
                static_cast<std::size_t>(
                    next_random() %
                    available
                );

            const std::string current_primer =
                random_reference.substr(
                    position,
                    random_primer_length
                );

            for (std::size_t budget = 0;
                 budget <= 3;
                 ++budget) {

                const auto branch_result =
                    random_branching
                        .search_5prime_mismatches(
                            current_primer,
                            random_anchor_length,
                            budget
                        );

                const auto candidate_result =
                    random_candidates.search(
                        current_primer,
                        random_anchor_length,
                        budget
                    );

                if (
                    collect_candidates(
                        candidate_result
                    ) !=
                    collect_branching(
                        random_index,
                        branch_result
                    )
                ) {
                    randomized_ok = false;
                    break;
                }
            }

            if (!randomized_ok) {
                break;
            }
        }

        expect(
            randomized_ok,
            "250 randomized branching/candidate trials"
        );

        /*
         * Lowercase.
         */
        const auto lowercase =
            candidates.search(
                "acgttgcaacgtacgt",
                12,
                2
            );

        expect(
            collect_candidates(
                lowercase
            ) ==
                std::vector<PositionMismatch>{
                    {3, 0},
                    {22, 1},
                    {41, 2}
                },
            "Lowercase candidate search"
        );

        bool excessive_budget_rejected = false;

        try {
            static_cast<void>(
                candidates.search(
                    primer,
                    12,
                    4
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            excessive_budget_rejected = true;
        }

        expect(
            excessive_budget_rejected,
            "Candidate mismatch >3 rejection"
        );

        std::cout
            << "All anchor candidate verification tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
