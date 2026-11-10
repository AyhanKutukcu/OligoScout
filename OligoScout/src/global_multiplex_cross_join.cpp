#include "primerpair/global_multiplex_cross_join.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace primerpair {

namespace {


struct ForwardHit {
    std::size_t primer_id{0};
    std::size_t primer_length{0};

    std::uint64_t position{0};

    std::size_t mismatches{0};
    std::uint64_t mismatch_mask{0};
};


struct ReverseHit {
    std::size_t primer_id{0};
    std::size_t primer_length{0};

    std::uint64_t position{0};
    std::uint64_t end_exclusive{0};

    std::size_t mismatches{0};
    std::uint64_t mismatch_mask{0};
};


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
    const GlobalMultiplexCrossProduct& lhs,
    const GlobalMultiplexCrossProduct& rhs
) noexcept {
    return
        std::tie(
            lhs.amplicon_start,
            lhs.amplicon_end_exclusive,

            lhs.forward_primer_id,
            lhs.reverse_primer_id,

            lhs.forward_position,
            lhs.reverse_position,

            lhs.forward_mismatches,
            lhs.reverse_mismatches,

            lhs.forward_mismatch_mask,
            lhs.reverse_mismatch_mask
        )
        <
        std::tie(
            rhs.amplicon_start,
            rhs.amplicon_end_exclusive,

            rhs.forward_primer_id,
            rhs.reverse_primer_id,

            rhs.forward_position,
            rhs.reverse_position,

            rhs.forward_mismatches,
            rhs.reverse_mismatches,

            rhs.forward_mismatch_mask,
            rhs.reverse_mismatch_mask
        );
}


}  // namespace


GlobalMultiplexCrossJoinResult
global_multiplex_cross_join(
    const std::vector<
        GlobalMultiplexPrimerHits
    >& primers,

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
            "Minimum amplicon length cannot exceed maximum."
        );
    }


    GlobalMultiplexCrossJoinResult result;

    result.stats.unique_primers =
        primers.size();


    std::vector<ForwardHit>
        forward_hits;

    std::vector<ReverseHit>
        reverse_hits;


    /*
     * ------------------------------------------------
     * Flatten all unique-primer hit lists.
     * ------------------------------------------------
     */
    for (
        const auto& primer :
        primers
    ) {
        if (
            primer.primer_length == 0
        ) {
            throw std::invalid_argument(
                "Primer length must be > 0."
            );
        }


        for (
            const auto& hit :
            primer.hits
        ) {
            if (
                hit.orientation ==
                PrimerOrientation::Forward
            ) {
                forward_hits.push_back(
                    ForwardHit{
                        primer.primer_id,
                        primer.primer_length,

                        hit.position,

                        hit.mismatches,
                        hit.mismatch_mask
                    }
                );

                continue;
            }


            if (
                hit.position >
                std::numeric_limits<
                    std::uint64_t
                >::max() -
                    primer.primer_length
            ) {
                continue;
            }


            reverse_hits.push_back(
                ReverseHit{
                    primer.primer_id,
                    primer.primer_length,

                    hit.position,

                    hit.position +
                        static_cast<std::uint64_t>(
                            primer.primer_length
                        ),

                    hit.mismatches,
                    hit.mismatch_mask
                }
            );
        }
    }


    result.stats.forward_hits =
        forward_hits.size();

    result.stats.reverse_hits =
        reverse_hits.size();


    /*
     * Forward hits genomic start'a göre.
     */
    std::sort(
        forward_hits.begin(),
        forward_hits.end(),
        [](
            const ForwardHit& lhs,
            const ForwardHit& rhs
        ) {
            return
                std::tie(
                    lhs.position,
                    lhs.primer_id,
                    lhs.mismatches,
                    lhs.mismatch_mask
                )
                <
                std::tie(
                    rhs.position,
                    rhs.primer_id,
                    rhs.mismatches,
                    rhs.mismatch_mask
                );
        }
    );


    /*
     * Reverse hits genomic END coordinate'e göre.
     *
     * Böylece amplicon-length window monoton
     * iki-pointer sweep ile yürütülebilir.
     */
    std::sort(
        reverse_hits.begin(),
        reverse_hits.end(),
        [](
            const ReverseHit& lhs,
            const ReverseHit& rhs
        ) {
            return
                std::tie(
                    lhs.end_exclusive,
                    lhs.position,
                    lhs.primer_id,
                    lhs.mismatches,
                    lhs.mismatch_mask
                )
                <
                std::tie(
                    rhs.end_exclusive,
                    rhs.position,
                    rhs.primer_id,
                    rhs.mismatches,
                    rhs.mismatch_mask
                );
        }
    );


    std::size_t lower_index = 0;
    std::size_t upper_index = 0;


    /*
     * ------------------------------------------------
     * GLOBAL MONOTONIC SWEEP
     * ------------------------------------------------
     */
    for (
        const auto& forward :
        forward_hits
    ) {
        const std::uint64_t
            minimum_reverse_end =
                saturating_add(
                    forward.position,
                    min_amplicon_length
                );


        const std::uint64_t
            maximum_reverse_end =
                saturating_add(
                    forward.position,
                    max_amplicon_length
                );


        while (
            lower_index <
                reverse_hits.size()
            &&
            reverse_hits
                .at(lower_index)
                .end_exclusive <
                minimum_reverse_end
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
                reverse_hits.size()
            &&
            reverse_hits
                .at(upper_index)
                .end_exclusive <=
                maximum_reverse_end
        ) {
            ++upper_index;
        }


        const std::uint64_t
            minimum_reverse_position =
                saturating_add(
                    forward.position,
                    static_cast<std::uint64_t>(
                        forward.primer_length
                    )
                );


        for (
            std::size_t index =
                lower_index;

            index <
                upper_index;

            ++index
        ) {
            const auto& reverse =
                reverse_hits.at(
                    index
                );


            ++result
                .stats
                .window_candidates;


            /*
             * Non-overlapping inward-facing PCR:
             *
             * FORWARD --->       <--- REVERSE
             */
            if (
                reverse.position <
                minimum_reverse_position
            ) {
                continue;
            }


            if (
                reverse.end_exclusive <=
                forward.position
            ) {
                continue;
            }


            const std::uint64_t
                amplicon_length =
                    reverse.end_exclusive -
                    forward.position;


            /*
             * Defensive validation.
             */
            if (
                amplicon_length <
                    min_amplicon_length
                ||
                amplicon_length >
                    max_amplicon_length
            ) {
                continue;
            }


            result.products.push_back(
                GlobalMultiplexCrossProduct{
                    forward.primer_id,
                    reverse.primer_id,

                    forward.position,
                    reverse.position,

                    forward.position,
                    reverse.end_exclusive,
                    amplicon_length,

                    forward.mismatches,
                    reverse.mismatches,

                    forward.mismatch_mask,
                    reverse.mismatch_mask
                }
            );


            ++result
                .stats
                .emitted_products;
        }
    }


    std::sort(
        result.products.begin(),
        result.products.end(),
        product_less
    );


    result.products.erase(
        std::unique(
            result.products.begin(),
            result.products.end()
        ),
        result.products.end()
    );


    result.stats.unique_products =
        result.products.size();


    return result;
}


}  // namespace primerpair
