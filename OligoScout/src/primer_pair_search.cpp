#include "primerpair/primer_pair_search.hpp"

#include <algorithm>
#include <limits>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

namespace primerpair {

namespace {

std::vector<OrientedPrimerSearchHit>
extract_orientation(
    const std::vector<
        OrientedPrimerSearchHit
    >& hits,
    const PrimerOrientation orientation
) {
    std::vector<
        OrientedPrimerSearchHit
    > output;

    for (const auto& hit : hits) {

        if (
            hit.orientation ==
            orientation
        ) {
            output.push_back(
                hit
            );
        }
    }

    /*
     * StrandAware engine already returns ordered
     * hits, but keep this helper self-contained.
     */
    std::sort(
        output.begin(),
        output.end(),
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

    return output;
}

std::uint64_t saturating_add(
    const std::uint64_t lhs,
    const std::uint64_t rhs
) noexcept {
    if (
        lhs >
        std::numeric_limits<
            std::uint64_t
        >::max() - rhs
    ) {
        return
            std::numeric_limits<
                std::uint64_t
            >::max();
    }

    return lhs + rhs;
}

bool product_less(
    const PrimerPairHit& lhs,
    const PrimerPairHit& rhs
) noexcept {
    if (
        lhs.amplicon_start !=
        rhs.amplicon_start
    ) {
        return
            lhs.amplicon_start <
            rhs.amplicon_start;
    }

    if (
        lhs.amplicon_end_exclusive !=
        rhs.amplicon_end_exclusive
    ) {
        return
            lhs.amplicon_end_exclusive <
            rhs.amplicon_end_exclusive;
    }

    if (
        lhs.left_primer !=
        rhs.left_primer
    ) {
        return
            static_cast<int>(
                lhs.left_primer
            ) <
            static_cast<int>(
                rhs.left_primer
            );
    }

    if (
        lhs.right_primer !=
        rhs.right_primer
    ) {
        return
            static_cast<int>(
                lhs.right_primer
            ) <
            static_cast<int>(
                rhs.right_primer
            );
    }

    if (
        lhs.left_mismatches !=
        rhs.left_mismatches
    ) {
        return
            lhs.left_mismatches <
            rhs.left_mismatches;
    }

    if (
        lhs.right_mismatches !=
        rhs.right_mismatches
    ) {
        return
            lhs.right_mismatches <
            rhs.right_mismatches;
    }

    if (
        lhs.left_mismatch_mask !=
        rhs.left_mismatch_mask
    ) {
        return
            lhs.left_mismatch_mask <
            rhs.left_mismatch_mask;
    }

    return
        lhs.right_mismatch_mask <
        rhs.right_mismatch_mask;
}

void unique_sorted_products(
    std::vector<PrimerPairHit>& products
) {
    /*
     * append_pairs() is expected to emit each
     * orientation stream in genomic order.
     *
     * Keep an explicit O(K) correctness guard.
     */
    if (
        !std::is_sorted(
            products.begin(),
            products.end(),
            product_less
        )
    ) {
        throw std::logic_error(
            "Primer-pair orientation stream "
            "is not coordinate-sorted."
        );
    }

    products.erase(
        std::unique(
            products.begin(),
            products.end()
        ),
        products.end()
    );
}

std::vector<PrimerPairHit>
merge_product_streams(
    std::vector<PrimerPairHit> first,
    std::vector<PrimerPairHit> second
) {
    unique_sorted_products(
        first
    );

    unique_sorted_products(
        second
    );

    std::vector<PrimerPairHit> merged;

    merged.reserve(
        first.size() +
        second.size()
    );

    /*
     * Both orientation streams are sorted.
     *
     * std::merge:
     *
     * O(K1 + K2)
     */
    std::merge(
        first.begin(),
        first.end(),

        second.begin(),
        second.end(),

        std::back_inserter(
            merged
        ),

        product_less
    );

    /*
     * Any cross-stream duplicates become adjacent
     * after merge and can be removed linearly.
     */
    merged.erase(
        std::unique(
            merged.begin(),
            merged.end()
        ),
        merged.end()
    );

    return merged;
}

}  // namespace

const char* to_string(
    const PrimerIdentity identity
) noexcept {
    switch (identity) {

        case PrimerIdentity::Primer1:
            return "PRIMER1";

        case PrimerIdentity::Primer2:
            return "PRIMER2";
    }

    return "UNKNOWN";
}

PrimerPairSearchEngine::
PrimerPairSearchEngine(
    const BidirectionalFMIndex& index,
    const PackedReference& reference,
    SearchDifficultyThresholds thresholds
)
    : strand_searcher_(
          index,
          reference,
          thresholds
      ),
      sensitive_searcher_(
          index,
          reference
      ) {
}

void PrimerPairSearchEngine::
append_pairs(
    const std::vector<
        OrientedPrimerSearchHit
    >& forward_hits,

    const std::size_t forward_primer_length,
    const PrimerIdentity forward_identity,

    const std::vector<
        OrientedPrimerSearchHit
    >& reverse_hits,

    const std::size_t reverse_primer_length,
    const PrimerIdentity reverse_identity,

    const std::uint64_t min_amplicon_length,
    const std::uint64_t max_amplicon_length,

    std::vector<PrimerPairHit>& output
) {
    if (
        forward_hits.empty() ||
        reverse_hits.empty()
    ) {
        return;
    }

    /*
     * If the maximum allowed product is shorter
     * than the reverse primer itself, no pair is
     * possible.
     */
    if (
        max_amplicon_length <
        reverse_primer_length
    ) {
        return;
    }

    const std::uint64_t reverse_length =
        static_cast<std::uint64_t>(
            reverse_primer_length
        );

    const std::uint64_t forward_length =
        static_cast<std::uint64_t>(
            forward_primer_length
        );

    /*
     * --------------------------------------------------
     * Monotonic sweep-line pairing
     *
     * forward_hits and reverse_hits are both sorted
     * by genomic coordinate.
     *
     * As forward positions increase, both the lower
     * and upper valid reverse-coordinate boundaries
     * can only move to the right.
     *
     * Therefore lower_index and upper_index never
     * move backwards.
     *
     * Complexity:
     *
     * O(F + R + K)
     *
     * where K is the number of emitted amplicons.
     * --------------------------------------------------
     */

    std::size_t lower_index = 0;
    std::size_t upper_index = 0;

    for (
        const auto& forward_hit :
        forward_hits
    ) {
        /*
         * Non-overlapping inward-facing geometry:
         *
         * FORWARD --->        <--- REVERSE
         */
        std::uint64_t minimum_reverse_position =
            saturating_add(
                forward_hit.position,
                forward_length
            );

        /*
         * Enforce minimum amplicon length.
         *
         * amplicon =
         *
         * reverse_position
         * + reverse_length
         * - forward_position
         */
        if (
            min_amplicon_length >
            reverse_length
        ) {
            const std::uint64_t
                by_min_amplicon =
                    saturating_add(
                        forward_hit.position,
                        min_amplicon_length -
                            reverse_length
                    );

            minimum_reverse_position =
                std::max(
                    minimum_reverse_position,
                    by_min_amplicon
                );
        }

        const std::uint64_t
            maximum_reverse_position =
                saturating_add(
                    forward_hit.position,
                    max_amplicon_length -
                        reverse_length
                );

        if (
            minimum_reverse_position >
            maximum_reverse_position
        ) {
            continue;
        }

        /*
         * Advance lower bound.
         *
         * This pointer never moves backwards.
         */
        while (
            lower_index <
                reverse_hits.size() &&
            reverse_hits
                .at(lower_index)
                .position <
                minimum_reverse_position
        ) {
            ++lower_index;
        }

        /*
         * upper_index cannot be left of lower_index.
         */
        if (
            upper_index <
            lower_index
        ) {
            upper_index =
                lower_index;
        }

        /*
         * Advance upper bound to the first reverse
         * hit outside the maximum allowed coordinate.
         *
         * This pointer also never moves backwards.
         */
        while (
            upper_index <
                reverse_hits.size() &&
            reverse_hits
                .at(upper_index)
                .position <=
                maximum_reverse_position
        ) {
            ++upper_index;
        }

        /*
         * Every reverse hit in:
         *
         * [lower_index, upper_index)
         *
         * is now inside the PCR-compatible genomic
         * coordinate window for this forward hit.
         */
        for (
            std::size_t index =
                lower_index;
            index <
                upper_index;
            ++index
        ) {
            const auto& reverse_hit =
                reverse_hits.at(
                    index
                );

            if (
                reverse_hit.position >
                std::numeric_limits<
                    std::uint64_t
                >::max() -
                    reverse_length
            ) {
                continue;
            }

            const std::uint64_t
                amplicon_end =
                    reverse_hit.position +
                    reverse_length;

            if (
                amplicon_end <=
                forward_hit.position
            ) {
                continue;
            }

            const std::uint64_t
                amplicon_length =
                    amplicon_end -
                    forward_hit.position;

            /*
             * Defensive validation.
             *
             * The sweep window should already
             * guarantee these bounds.
             */
            if (
                amplicon_length <
                    min_amplicon_length ||
                amplicon_length >
                    max_amplicon_length
            ) {
                continue;
            }

            output.push_back(
                PrimerPairHit{
                    forward_identity,
                    reverse_identity,

                    forward_hit.position,
                    reverse_hit.position,

                    forward_hit.mismatches,
                    reverse_hit.mismatches,

                    forward_hit.position,
                    amplicon_end,
                    amplicon_length,

                    forward_hit.mismatch_mask,
                    reverse_hit.mismatch_mask
                }
            );
        }
    }
}

PrimerPairSearchResult
PrimerPairSearchEngine::assemble_pairs(
    const std::string_view primer1,

    const std::vector<
        OrientedPrimerSearchHit
    >& primer1_hits,

    const std::string_view primer2,

    const std::vector<
        OrientedPrimerSearchHit
    >& primer2_hits,

    const std::uint64_t min_amplicon_length,
    const std::uint64_t max_amplicon_length
) {
    if (
        min_amplicon_length == 0
    ) {
        throw std::invalid_argument(
            "Minimum amplicon length must be > 0."
        );
    }

    if (
        min_amplicon_length >
        max_amplicon_length
    ) {
        throw std::invalid_argument(
            "Minimum amplicon length cannot "
            "exceed maximum amplicon length."
        );
    }

    const auto primer1_forward =
        extract_orientation(
            primer1_hits,
            PrimerOrientation::Forward
        );

    const auto primer1_reverse =
        extract_orientation(
            primer1_hits,
            PrimerOrientation::Reverse
        );

    const auto primer2_forward =
        extract_orientation(
            primer2_hits,
            PrimerOrientation::Forward
        );

    const auto primer2_reverse =
        extract_orientation(
            primer2_hits,
            PrimerOrientation::Reverse
        );

    /*
     * Keep the two biological configurations
     * as independent coordinate-sorted streams.
     */
    std::vector<PrimerPairHit>
        primer1_to_primer2;

    std::vector<PrimerPairHit>
        primer2_to_primer1;

    /*
     * Primer1 --->      <--- Primer2
     */
    append_pairs(
        primer1_forward,
        primer1.size(),
        PrimerIdentity::Primer1,

        primer2_reverse,
        primer2.size(),
        PrimerIdentity::Primer2,

        min_amplicon_length,
        max_amplicon_length,

        primer1_to_primer2
    );

    /*
     * Primer2 --->      <--- Primer1
     */
    append_pairs(
        primer2_forward,
        primer2.size(),
        PrimerIdentity::Primer2,

        primer1_reverse,
        primer1.size(),
        PrimerIdentity::Primer1,

        min_amplicon_length,
        max_amplicon_length,

        primer2_to_primer1
    );

    /*
     * Same validated O(K1 + K2) merge.
     */
    std::vector<PrimerPairHit> products =
        merge_product_streams(
            std::move(
                primer1_to_primer2
            ),
            std::move(
                primer2_to_primer1
            )
        );

    return PrimerPairSearchResult{
        primer1.size(),
        primer2.size(),

        primer1_hits.size(),
        primer2_hits.size(),

        min_amplicon_length,
        max_amplicon_length,

        std::move(
            products
        )
    };
}


PrimerPairSearchResult
PrimerPairSearchEngine::search(
    const std::string_view primer1,
    const std::string_view primer2,
    const std::size_t anchor_length,
    const std::size_t max_mismatches,
    const std::uint64_t min_amplicon_length,
    const std::uint64_t max_amplicon_length
) const {
    /*
     * Legacy STRICT path.
     *
     * Candidate generation is unchanged.
     */
    const auto primer1_result =
        strand_searcher_.search(
            primer1,
            anchor_length,
            max_mismatches
        );

    const auto primer2_result =
        strand_searcher_.search(
            primer2,
            anchor_length,
            max_mismatches
        );

    return assemble_pairs(
        primer1,
        primer1_result.hits,

        primer2,
        primer2_result.hits,

        min_amplicon_length,
        max_amplicon_length
    );
}


PrimerPairSearchResult
PrimerPairSearchEngine::search(
    const std::string_view primer1,
    const std::string_view primer2,
    const SearchProfile profile,
    const std::size_t anchor_length,
    const std::size_t max_mismatches,
    const std::uint64_t min_amplicon_length,
    const std::uint64_t max_amplicon_length
) const {
    switch (profile) {

        case SearchProfile::Strict:

            /*
             * Delegate to legacy STRICT API so the
             * existing validated behavior remains
             * exactly one implementation.
             */
            return search(
                primer1,
                primer2,
                anchor_length,
                max_mismatches,
                min_amplicon_length,
                max_amplicon_length
            );


        case SearchProfile::Sensitive: {
            /*
             * Full-primer adaptive SENSITIVE search.
             *
             * anchor_length intentionally has no
             * meaning for this profile.
             */
            static_cast<void>(
                anchor_length
            );

            const auto primer1_result =
                sensitive_searcher_.search(
                    primer1,
                    max_mismatches
                );

            const auto primer2_result =
                sensitive_searcher_.search(
                    primer2,
                    max_mismatches
                );

            return assemble_pairs(
                primer1,
                primer1_result
                    .search_result
                    .hits,

                primer2,
                primer2_result
                    .search_result
                    .hits,

                min_amplicon_length,
                max_amplicon_length
            );
        }
    }

    throw std::logic_error(
        "Unknown search profile."
    );
}

}  // namespace primerpair
