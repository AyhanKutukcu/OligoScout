#include "primerpair/single_primer_search.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace primerpair {

namespace {

std::string normalize_primer(
    const std::string_view primer
) {
    if (primer.empty()) {
        throw std::invalid_argument(
            "Primer cannot be empty."
        );
    }

    std::string normalized;
    normalized.reserve(
        primer.size()
    );

    for (const char raw : primer) {

        const char base =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        raw
                    )
                )
            );

        switch (base) {

            case 'A':
            case 'C':
            case 'G':
            case 'T':
                normalized.push_back(base);
                break;

            default:
                throw std::invalid_argument(
                    "Primer must contain only A/C/G/T."
                );
        }
    }

    return normalized;
}

void normalize_hits(
    std::vector<PrimerSearchHit>& hits
) {
    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const PrimerSearchHit& lhs,
            const PrimerSearchHit& rhs
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
}

struct BranchState {
    BidirectionalInterval state{};
    std::size_t mismatches{0};
};

}  // namespace

SinglePrimerSearchEngine::
SinglePrimerSearchEngine(
    const BidirectionalFMIndex& index,
    const PackedReference& reference,
    SearchDifficultyThresholds thresholds
)
    : index_(index),
      reference_(reference),
      router_(
          index,
          thresholds
      ) {
}

std::vector<PrimerSearchHit>
SinglePrimerSearchEngine::
execute_branching(
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches,
    const BidirectionalInterval& anchor_state
) const {
    const std::string normalized =
        normalize_primer(
            primer
        );

    if (
        anchor_length == 0 ||
        anchor_length >
            normalized.size()
    ) {
        throw std::invalid_argument(
            "Invalid anchor length."
        );
    }

    if (anchor_state.empty()) {
        return {};
    }

    const std::size_t prefix_length =
        normalized.size() -
        anchor_length;

    std::vector<BranchState> current;

    current.push_back(
        BranchState{
            anchor_state,
            0
        }
    );

    constexpr std::array<char, 4>
        alphabet{
            'A',
            'C',
            'G',
            'T'
        };

    /*
     * Extend from the nucleotide immediately
     * left of the exact 3-prime anchor toward
     * the 5-prime end.
     */
    for (std::size_t remaining =
             prefix_length;
         remaining > 0;
         --remaining) {

        const std::size_t primer_position =
            remaining - 1;

        const char expected =
            normalized.at(
                primer_position
            );

        std::vector<BranchState> next;

        next.reserve(
            current.size() * 4
        );

        for (const auto& branch : current) {

            /*
             * Mismatch budget tamamen dolmuşsa yalnızca
             * primerde beklenen baz legal child olabilir.
             *
             * Bu durumda dört child üretmek yerine eski
             * scalar extension yolu daha ucuzdur.
             */
            if (
                branch.mismatches ==
                max_mismatches
            ) {
                const BidirectionalInterval
                    next_state =
                        index_.extend_left(
                            branch.state,
                            expected
                        );

                if (!next_state.empty()) {
                    next.push_back(
                        BranchState{
                            next_state,
                            branch.mismatches
                        }
                    );
                }

                continue;
            }

            /*
             * Burada:
             *
             *   branch.mismatches < max_mismatches
             *
             * olduğundan A/C/G/T çocuklarının tamamı
             * mismatch budget açısından potansiyel olarak
             * geçerlidir.
             *
             * Dört ayrı extend_left() yerine parent için
             * tek toplu rank/extension çağrısı kullanılır.
             *
             * children sırası:
             *   0=A, 1=C, 2=G, 3=T
             */
            const auto children =
                index_.extend_left_all(
                    branch.state
                );

            for (
                std::size_t base_index = 0;
                base_index < alphabet.size();
                ++base_index
            ) {
                const char base =
                    alphabet.at(
                        base_index
                    );

                const std::size_t mismatches =
                    branch.mismatches +
                    (
                        base == expected
                            ? 0
                            : 1
                    );

                /*
                 * branch.mismatches < max_mismatches
                 * olduğu için normalde bu kontrol her
                 * zaman geçer. Defensive invariant olarak
                 * tutuluyor.
                 */
                if (
                    mismatches >
                    max_mismatches
                ) {
                    continue;
                }

                const BidirectionalInterval&
                    next_state =
                        children.at(
                            base_index
                        );

                if (next_state.empty()) {
                    continue;
                }

                next.push_back(
                    BranchState{
                        next_state,
                        mismatches
                    }
                );
            }
        }


        current =
            std::move(
                next
            );

        if (current.empty()) {
            break;
        }
    }

    std::vector<PrimerSearchHit> hits;

    for (const auto& branch : current) {

        const auto positions =
            index_.locate_unsorted(
                branch.state
            );

        for (const auto position : positions) {

            hits.push_back(
                PrimerSearchHit{
                    position,
                    branch.mismatches
                }
            );
        }
    }

    normalize_hits(
        hits
    );

    return hits;
}

std::vector<PrimerSearchHit>
SinglePrimerSearchEngine::
execute_candidate(
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches,
    const BidirectionalInterval& anchor_state
) const {
    const std::string normalized =
        normalize_primer(
            primer
        );

    if (
        anchor_length == 0 ||
        anchor_length >
            normalized.size()
    ) {
        throw std::invalid_argument(
            "Invalid anchor length."
        );
    }

    if (anchor_state.empty()) {
        return {};
    }

    const std::size_t prefix_length =
        normalized.size() -
        anchor_length;

    std::vector<std::uint64_t>
        anchor_positions =
            index_.locate(
                anchor_state
            );

    /*
     * Defensive coordinate normalization.
     */
    std::sort(
        anchor_positions.begin(),
        anchor_positions.end()
    );

    anchor_positions.erase(
        std::unique(
            anchor_positions.begin(),
            anchor_positions.end()
        ),
        anchor_positions.end()
    );

    const std::string_view prefix(
        normalized.data(),
        prefix_length
    );

    std::vector<PrimerSearchHit> hits;

    hits.reserve(
        anchor_positions.size()
    );

    for (
        const std::uint64_t anchor_position :
        anchor_positions
    ) {

        if (
            anchor_position <
            prefix_length
        ) {
            continue;
        }

        const std::uint64_t candidate_start =
            anchor_position -
            static_cast<std::uint64_t>(
                prefix_length
            );

        if (
            candidate_start >
            reference_.size()
        ) {
            continue;
        }

        if (
            normalized.size() >
            reference_.size() -
                candidate_start
        ) {
            continue;
        }

        std::size_t mismatches = 0;

        if (prefix_length != 0) {

            mismatches =
                reference_
                    .bounded_hamming_distance(
                        candidate_start,
                        prefix,
                        max_mismatches
                    );

            if (
                mismatches >
                max_mismatches
            ) {
                continue;
            }
        }

        hits.push_back(
            PrimerSearchHit{
                candidate_start,
                mismatches
            }
        );
    }

    normalize_hits(
        hits
    );

    return hits;
}

SinglePrimerSearchResult
SinglePrimerSearchEngine::search(
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    const SearchStrategyDecision decision =
        router_.decide(
            primer,
            anchor_length,
            max_mismatches
        );

    /*
     * This interval was already computed once by
     * SearchDifficultyEstimator.
     *
     * Both routed backends reuse it.
     */
    const BidirectionalInterval&
        anchor_state =
            decision
                .difficulty_profile
                .anchor_state;

    std::vector<PrimerSearchHit> hits;

    switch (
        decision.recommended_strategy
    ) {

        case SearchStrategy::DirectBranching:

            hits =
                execute_branching(
                    primer,
                    anchor_length,
                    max_mismatches,
                    anchor_state
                );

            break;

        case SearchStrategy::
            AnchorCandidateVerification:

            hits =
                execute_candidate(
                    primer,
                    anchor_length,
                    max_mismatches,
                    anchor_state
                );

            break;

        case SearchStrategy::
            SplitSeedCandidate:

            throw std::logic_error(
                "SplitSeedCandidate backend "
                "is not implemented."
            );
    }

    return SinglePrimerSearchResult{
        decision,
        std::move(
            hits
        )
    };
}

}  // namespace primerpair
