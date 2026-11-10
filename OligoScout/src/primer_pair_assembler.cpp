#include "primerpair/primer_pair_assembler.hpp"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
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

    output.reserve(
        hits.size()
    );

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
     * Hybrid strand results are already normalized,
     * but the assembler deliberately guarantees its
     * own monotonic-coordinate precondition.
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

    output.erase(
        std::unique(
            output.begin(),
            output.end()
        ),
        output.end()
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


void append_pairs(
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
     * Monotonic sweep-line.
     *
     * Both lower_index and upper_index only move
     * toward increasing genomic coordinates.
     */
    std::size_t lower_index = 0;
    std::size_t upper_index = 0;


    for (
        const auto& forward_hit :
        forward_hits
    ) {
        /*
         * Prevent overlapping/outward geometry:
         *
         * FORWARD --->       <--- REVERSE
         */
        std::uint64_t
            minimum_reverse_position =
                saturating_add(
                    forward_hit.position,
                    forward_length
                );


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


        if (
            upper_index <
            lower_index
        ) {
            upper_index =
                lower_index;
        }


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


void normalize_products(
    std::vector<
        PrimerPairHit
    >& products
) {
    std::sort(
        products.begin(),
        products.end(),
        product_less
    );

    products.erase(
        std::unique(
            products.begin(),
            products.end()
        ),
        products.end()
    );
}


}  // namespace


PrimerPairSearchResult
assemble_primer_pair_hits(
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


    std::vector<PrimerPairHit>
        products;


    /*
     * Configuration 1:
     *
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

        products
    );


    /*
     * Configuration 2:
     *
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

        products
    );


    normalize_products(
        products
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

}  // namespace primerpair
