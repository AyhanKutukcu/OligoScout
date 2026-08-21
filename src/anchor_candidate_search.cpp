#include "primerpair/anchor_candidate_search.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <utility>

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

}  // namespace

AnchorCandidateSearcher::
AnchorCandidateSearcher(
    const BidirectionalFMIndex& index,
    const PackedReference& reference
) noexcept
    : index_(index),
      reference_(reference) {
}

AnchorCandidateSearchResult
AnchorCandidateSearcher::search(
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    if (
        max_mismatches >
        3
    ) {
        throw std::invalid_argument(
            "Current MVP supports at most "
            "3 mismatches per primer."
        );
    }


    const std::string normalized =
        normalize_primer(
            primer
        );


    if (
        anchor_length ==
        0
    ) {
        throw std::invalid_argument(
            "Anchor length must be > 0."
        );
    }


    if (
        anchor_length >
        normalized.size()
    ) {
        throw std::invalid_argument(
            "Anchor length cannot exceed "
            "primer length."
        );
    }


    const std::size_t prefix_length =
        normalized.size()
        -
        anchor_length;


    const std::string_view anchor(
        normalized.data()
            +
            prefix_length,
        anchor_length
    );


    const BidirectionalInterval
        anchor_state =
            index_.search(
                anchor
            );


    std::vector<std::uint64_t>
        anchor_positions;


    if (
        !anchor_state.empty()
    ) {
        anchor_positions =
            index_.locate(
                anchor_state
            );
    }


    /*
     * Both the legacy BiFM path and the new
     * IP-BWT batch path now converge on the
     * exact same candidate-verification code.
     */
    return
        verify_from_anchor_positions(
            normalized,
            std::move(
                anchor_positions
            ),
            anchor_state.size(),
            anchor_length,
            max_mismatches
        );
}


AnchorCandidateSearchResult
AnchorCandidateSearcher::
verify_from_anchor_positions(
    const std::string_view primer,
    std::vector<std::uint64_t> anchor_positions,
    const std::uint64_t anchor_occurrences,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    if (
        max_mismatches >
        3
    ) {
        throw std::invalid_argument(
            "Current MVP supports at most "
            "3 mismatches per primer."
        );
    }


    const std::string normalized =
        normalize_primer(
            primer
        );


    if (
        anchor_length ==
        0
    ) {
        throw std::invalid_argument(
            "Anchor length must be > 0."
        );
    }


    if (
        anchor_length >
        normalized.size()
    ) {
        throw std::invalid_argument(
            "Anchor length cannot exceed "
            "primer length."
        );
    }


    const std::size_t prefix_length =
        normalized.size()
        -
        anchor_length;


    AnchorCandidateSearchResult result;


    result.primer_length =
        normalized.size();

    result.anchor_length =
        anchor_length;

    result.max_mismatches =
        max_mismatches;

    result.anchor_occurrences =
        anchor_occurrences;


    if (
        anchor_positions.empty()
    ) {
        return result;
    }


    /*
     * Preserve the same defensive deduplication
     * performed by the original BiFM-backed path.
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


    for (
        const std::uint64_t anchor_position :
        anchor_positions
    ) {
        /*
         * Primer start would become negative.
         */
        if (
            anchor_position <
            prefix_length
        ) {
            continue;
        }


        const std::uint64_t candidate_start =
            anchor_position
            -
            static_cast<std::uint64_t>(
                prefix_length
            );


        /*
         * Entire primer must fit in the shard.
         */
        if (
            candidate_start >
            reference_.size()
        ) {
            continue;
        }


        if (
            normalized.size() >
            reference_.size()
            -
            candidate_start
        ) {
            continue;
        }


        ++result.candidates_verified;


        std::size_t mismatches = 0;


        /*
         * STRICT profile:
         *
         * exact anchor is not re-verified;
         * mismatches are counted only in the
         * sequence to the left of that anchor.
         */
        if (
            prefix_length !=
            0
        ) {
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


        result.hits.push_back(
            AnchorCandidateHit{
                candidate_start,
                mismatches
            }
        );
    }


    return result;
}



}  // namespace primerpair
