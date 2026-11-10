#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <cctype>
#include <stdexcept>
#include <utility>
#include <vector>

namespace primerpair {

namespace {

struct BranchState {
    BidirectionalInterval state{};
    std::size_t mismatches{0};
};

/*
 * Build a mismatch bit mask directly against the
 * packed reference.
 *
 * query is always represented left-to-right on the
 * reference.
 *
 * If reverse_to_original is true, query is the
 * reverse-complement representation of the original
 * primer. Query coordinate i must therefore map to:
 *
 * original_position = length - 1 - i
 */
std::uint64_t build_mismatch_mask(
    const PackedReference& reference,
    const std::uint64_t start,
    const std::string_view query,
    const bool reverse_to_original
) {
    if (query.size() > 64) {
        throw std::invalid_argument(
            "Mismatch mask supports primers up to 64 nt."
        );
    }

    if (
        start >
        reference.size()
    ) {
        throw std::out_of_range(
            "Mismatch-mask start outside reference."
        );
    }

    if (
        query.size() >
        reference.size() - start
    ) {
        throw std::out_of_range(
            "Mismatch-mask query exceeds reference."
        );
    }

    std::uint64_t mask = 0;

    for (
        std::size_t query_position = 0;
        query_position < query.size();
        ++query_position
    ) {
        const char query_base =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        query.at(
                            query_position
                        )
                    )
                )
            );

        const char reference_base =
            reference.base_at(
                start +
                static_cast<std::uint64_t>(
                    query_position
                )
            );

        if (
            reference_base ==
            query_base
        ) {
            continue;
        }

        const std::size_t
            original_position =
                reverse_to_original
                    ? query.size() -
                        1 -
                        query_position
                    : query_position;

        mask |=
            (
                std::uint64_t{1}
                <<
                original_position
            );
    }

    return mask;
}

void validate_mismatch_count(
    const std::uint64_t mask,
    const std::size_t expected
) {
    const std::size_t observed =
        static_cast<std::size_t>(
            std::popcount(
                mask
            )
        );

    if (observed != expected) {
        throw std::logic_error(
            "Mismatch count and mismatch mask disagree."
        );
    }
}

char complement(
    const char raw
) {
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
            return 'T';

        case 'C':
            return 'G';

        case 'G':
            return 'C';

        case 'T':
            return 'A';

        default:
            throw std::invalid_argument(
                "Primer must contain only A/C/G/T."
            );
    }
}

void normalize_hits(
    std::vector<OrientedPrimerSearchHit>& hits
) {
    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const OrientedPrimerSearchHit& lhs,
            const OrientedPrimerSearchHit& rhs
        ) {
            if (
                lhs.position !=
                rhs.position
            ) {
                return
                    lhs.position <
                    rhs.position;
            }

            if (
                lhs.orientation !=
                rhs.orientation
            ) {
                return
                    static_cast<int>(
                        lhs.orientation
                    ) <
                    static_cast<int>(
                        rhs.orientation
                    );
            }

            if (
                lhs.mismatches !=
                rhs.mismatches
            ) {
                return
                    lhs.mismatches <
                    rhs.mismatches;
            }

            return
                lhs.mismatch_mask <
                rhs.mismatch_mask;
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

}  // namespace

std::string reverse_complement(
    const std::string_view sequence
) {
    if (sequence.empty()) {
        throw std::invalid_argument(
            "Primer cannot be empty."
        );
    }

    std::string output;

    output.reserve(
        sequence.size()
    );

    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {
        output.push_back(
            complement(
                *it
            )
        );
    }

    return output;
}

const char* to_string(
    const PrimerOrientation orientation
) noexcept {
    switch (orientation) {

        case PrimerOrientation::Forward:
            return "FORWARD";

        case PrimerOrientation::Reverse:
            return "REVERSE";
    }

    return "UNKNOWN";
}

StrandAwarePrimerSearchEngine::
StrandAwarePrimerSearchEngine(
    const BidirectionalFMIndex& index,
    const PackedReference& reference,
    SearchDifficultyThresholds thresholds
)
    : index_(index),
      reference_(reference),
      router_(
          index,
          thresholds
      ),
      forward_engine_(
          index,
          reference,
          thresholds
      ) {
}

SearchStrategyDecision
StrandAwarePrimerSearchEngine::
make_reverse_decision(
    const std::string_view reverse_query,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    if (
        anchor_length == 0 ||
        anchor_length >
            reverse_query.size()
    ) {
        throw std::invalid_argument(
            "Invalid anchor length."
        );
    }

    /*
     * CRITICAL:
     *
     * reverse_query = reverse-complement(original primer)
     *
     * Therefore the biological 3-prime end of the
     * original primer is at the LEFT side of
     * reverse_query.
     *
     * Exact anchor = PREFIX, not suffix.
     */
    const std::string_view anchor(
        reverse_query.data(),
        anchor_length
    );

    const BidirectionalInterval anchor_state =
        index_.search(
            anchor
        );

    const std::uint64_t occurrences =
        anchor_state.size();

    const SearchDifficulty difficulty =
        router_
            .estimator()
            .classify(
                occurrences
            );

    SearchDifficultyProfile profile{
        difficulty,
        anchor_state,
        occurrences,
        reverse_query.size(),
        anchor_length
    };

    return SearchStrategyDecision{
        std::move(
            profile
        ),
        SearchStrategyRouter::choose(
            difficulty,
            max_mismatches
        ),
        max_mismatches
    };
}

std::vector<OrientedPrimerSearchHit>
StrandAwarePrimerSearchEngine::
execute_reverse_branching(
    const std::string_view reverse_query,
    const std::size_t anchor_length,
    const std::size_t max_mismatches,
    const BidirectionalInterval& anchor_state
) const {
    if (anchor_state.empty()) {
        return {};
    }

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
     * The prefix anchor represents the original
     * primer's exact 3-prime region.
     *
     * Move toward the original primer's 5-prime end
     * by extending RIGHT through reverse_query.
     */
    for (
        std::size_t position = anchor_length;
        position < reverse_query.size();
        ++position
    ) {
        const char expected =
            reverse_query.at(
                position
            );

        std::vector<BranchState> next;

        next.reserve(
            current.size() * 4
        );

        for (const auto& branch : current) {

            for (const char base : alphabet) {

                const std::size_t mismatches =
                    branch.mismatches +
                    (
                        base == expected
                            ? 0
                            : 1
                    );

                if (
                    mismatches >
                    max_mismatches
                ) {
                    continue;
                }

                const BidirectionalInterval
                    next_state =
                        index_.extend_right(
                            branch.state,
                            base
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

    std::vector<OrientedPrimerSearchHit> hits;

    for (const auto& branch : current) {

        const auto positions =
            index_.locate_unsorted(
                branch.state
            );

        for (const auto position : positions) {

            const std::uint64_t mismatch_mask =
                build_mismatch_mask(
                    reference_,
                    position,
                    reverse_query,
                    true
                );

            validate_mismatch_count(
                mismatch_mask,
                branch.mismatches
            );

            hits.push_back(
                OrientedPrimerSearchHit{
                    position,
                    branch.mismatches,
                    PrimerOrientation::Reverse,
                    mismatch_mask
                }
            );
        }
    }

    normalize_hits(
        hits
    );

    return hits;
}

std::vector<OrientedPrimerSearchHit>
StrandAwarePrimerSearchEngine::
execute_reverse_candidate(
    const std::string_view reverse_query,
    const std::size_t anchor_length,
    const std::size_t max_mismatches,
    const BidirectionalInterval& anchor_state
) const {
    if (anchor_state.empty()) {
        return {};
    }

    std::vector<std::uint64_t>
        anchor_positions =
            index_.locate(
                anchor_state
            );

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

    const std::size_t suffix_length =
        reverse_query.size() -
        anchor_length;

    const std::string_view suffix(
        reverse_query.data() +
            anchor_length,
        suffix_length
    );

    std::vector<OrientedPrimerSearchHit> hits;

    hits.reserve(
        anchor_positions.size()
    );

    for (
        const std::uint64_t candidate_start :
        anchor_positions
    ) {
        if (
            candidate_start >
            reference_.size()
        ) {
            continue;
        }

        if (
            reverse_query.size() >
            reference_.size() -
                candidate_start
        ) {
            continue;
        }

        std::size_t mismatches = 0;

        if (suffix_length != 0) {

            const std::uint64_t suffix_start =
                candidate_start +
                static_cast<std::uint64_t>(
                    anchor_length
                );

            mismatches =
                reference_
                    .bounded_hamming_distance(
                        suffix_start,
                        suffix,
                        max_mismatches
                    );

            if (
                mismatches >
                max_mismatches
            ) {
                continue;
            }
        }

        const std::uint64_t mismatch_mask =
            build_mismatch_mask(
                reference_,
                candidate_start,
                reverse_query,
                true
            );

        validate_mismatch_count(
            mismatch_mask,
            mismatches
        );

        hits.push_back(
            OrientedPrimerSearchHit{
                candidate_start,
                mismatches,
                PrimerOrientation::Reverse,
                mismatch_mask
            }
        );
    }

    normalize_hits(
        hits
    );

    return hits;
}

StrandAwarePrimerSearchResult
StrandAwarePrimerSearchEngine::search(
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    /*
     * Forward orientation:
     * existing validated suffix-anchor engine.
     */
    const SinglePrimerSearchResult forward =
        forward_engine_.search(
            primer,
            anchor_length,
            max_mismatches
        );

    const std::string reverse_query =
        reverse_complement(
            primer
        );

    /*
     * Reverse orientation requires a PREFIX anchor.
     */
    const SearchStrategyDecision reverse_decision =
        make_reverse_decision(
            reverse_query,
            anchor_length,
            max_mismatches
        );

    std::vector<OrientedPrimerSearchHit> hits;

    hits.reserve(
        forward.hits.size()
    );

    for (const auto& hit : forward.hits) {

        const std::uint64_t mismatch_mask =
            build_mismatch_mask(
                reference_,
                hit.position,
                primer,
                false
            );

        validate_mismatch_count(
            mismatch_mask,
            hit.mismatches
        );

        hits.push_back(
            OrientedPrimerSearchHit{
                hit.position,
                hit.mismatches,
                PrimerOrientation::Forward,
                mismatch_mask
            }
        );
    }

    std::vector<OrientedPrimerSearchHit>
        reverse_hits;

    switch (
        reverse_decision.recommended_strategy
    ) {

        case SearchStrategy::DirectBranching:

            reverse_hits =
                execute_reverse_branching(
                    reverse_query,
                    anchor_length,
                    max_mismatches,
                    reverse_decision
                        .difficulty_profile
                        .anchor_state
                );

            break;

        case SearchStrategy::
            AnchorCandidateVerification:

            reverse_hits =
                execute_reverse_candidate(
                    reverse_query,
                    anchor_length,
                    max_mismatches,
                    reverse_decision
                        .difficulty_profile
                        .anchor_state
                );

            break;

        case SearchStrategy::
            SplitSeedCandidate:

            throw std::logic_error(
                "SplitSeedCandidate backend "
                "is not implemented."
            );
    }

    hits.insert(
        hits.end(),
        reverse_hits.begin(),
        reverse_hits.end()
    );

    normalize_hits(
        hits
    );

    return StrandAwarePrimerSearchResult{
        primer.size(),
        forward.decision,
        reverse_decision,
        std::move(
            hits
        )
    };
}

}  // namespace primerpair
