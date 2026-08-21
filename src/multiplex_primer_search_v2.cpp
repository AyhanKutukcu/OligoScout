#include "primerpair/multiplex_primer_search_v2.hpp"

#include "primerpair/primer_pair_assembler.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace primerpair {

namespace {


struct SlotMembership {
    std::size_t pair_index{0};

    MultiplexPrimerSlot slot{
        MultiplexPrimerSlot::Primer1
    };
};


std::string canonical_primer(
    const std::string_view primer
) {
    std::string output;

    output.reserve(
        primer.size()
    );

    for (const char raw : primer) {
        output.push_back(
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        raw
                    )
                )
            )
        );
    }

    return output;
}


bool cross_less(
    const MultiplexCrossAmplicon& lhs,
    const MultiplexCrossAmplicon& rhs
) noexcept {
    return
        std::tie(
            lhs.forward_pair_index,
            lhs.forward_slot,
            lhs.reverse_pair_index,
            lhs.reverse_slot,

            lhs.amplicon_start,
            lhs.amplicon_end_exclusive,
            lhs.amplicon_length,

            lhs.forward_mismatches,
            lhs.reverse_mismatches,

            lhs.forward_mismatch_mask,
            lhs.reverse_mismatch_mask
        )
        <
        std::tie(
            rhs.forward_pair_index,
            rhs.forward_slot,
            rhs.reverse_pair_index,
            rhs.reverse_slot,

            rhs.amplicon_start,
            rhs.amplicon_end_exclusive,
            rhs.amplicon_length,

            rhs.forward_mismatches,
            rhs.reverse_mismatches,

            rhs.forward_mismatch_mask,
            rhs.reverse_mismatch_mask
        );
}


void normalize_cross_products(
    std::vector<
        MultiplexCrossAmplicon
    >& products
) {
    std::sort(
        products.begin(),
        products.end(),
        cross_less
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


MultiplexPrimerSearchEngineV2::
MultiplexPrimerSearchEngineV2(
    const HybridStrandAwarePrimerSearchEngine&
        primer_engine
) noexcept
    : primer_engine_(
          primer_engine
      ) {
}


MultiplexPrimerSearchResultV2
MultiplexPrimerSearchEngineV2::search(
    const std::vector<
        MultiplexPrimerPairRequest
    >& requests,

    const std::size_t anchor_length,

    const bool include_cross_pairs,

    const std::uint64_t
        cross_min_amplicon_length,

    const std::uint64_t
        cross_max_amplicon_length
) const {
    MultiplexPrimerSearchResultV2 result;


    result.stats.pair_requests =
        requests.size();

    result.stats.total_primer_slots =
        requests.size() * 2;


    if (requests.empty()) {
        return result;
    }


    if (
        include_cross_pairs &&
        cross_min_amplicon_length == 0
    ) {
        throw std::invalid_argument(
            "Cross-pair minimum amplicon "
            "length must be > 0."
        );
    }


    if (
        include_cross_pairs &&
        cross_min_amplicon_length >
            cross_max_amplicon_length
    ) {
        throw std::invalid_argument(
            "Cross-pair minimum amplicon "
            "length exceeds maximum."
        );
    }


    using PrimerKey =
        std::pair<
            std::string,
            std::size_t
        >;


    std::map<
        PrimerKey,
        std::size_t
    > unique_lookup;


    std::vector<std::string>
        unique_primers;


    std::vector<std::size_t>
        unique_mismatch_budgets;


    std::vector<
        std::array<std::size_t, 2>
    > slot_unique_ids(
        requests.size()
    );


    unique_primers.reserve(
        requests.size() * 2
    );

    unique_mismatch_budgets.reserve(
        requests.size() * 2
    );


    /*
     * ================================================
     * Unique primer table
     * ================================================
     */
    for (
        std::size_t pair_index = 0;
        pair_index < requests.size();
        ++pair_index
    ) {
        const auto& request =
            requests.at(
                pair_index
            );


        const std::array<
            std::string_view,
            2
        > sequences{
            request.primer1,
            request.primer2
        };


        for (
            std::size_t slot = 0;
            slot < 2;
            ++slot
        ) {
            const std::string sequence =
                canonical_primer(
                    sequences.at(slot)
                );


            const PrimerKey key{
                sequence,
                request.max_mismatches
            };


            const auto found =
                unique_lookup.find(
                    key
                );


            if (
                found !=
                unique_lookup.end()
            ) {
                slot_unique_ids
                    .at(pair_index)
                    .at(slot) =
                        found->second;

                continue;
            }


            const std::size_t uid =
                unique_primers.size();


            unique_primers.push_back(
                sequence
            );


            unique_mismatch_budgets
                .push_back(
                    request.max_mismatches
                );


            unique_lookup.emplace(
                key,
                uid
            );


            slot_unique_ids
                .at(pair_index)
                .at(slot) =
                    uid;
        }
    }


    result.stats.unique_primer_queries =
        unique_primers.size();


    result.stats.reused_primer_slots =
        result.stats.total_primer_slots -
        result.stats.unique_primer_queries;


    /*
     * ================================================
     * ONE hybrid both-strand batch
     * ================================================
     */
    std::vector<
        HybridStrandAwarePrimerRequest
    > unique_requests;


    unique_requests.reserve(
        unique_primers.size()
    );


    for (
        std::size_t uid = 0;
        uid < unique_primers.size();
        ++uid
    ) {
        unique_requests.push_back(
            HybridStrandAwarePrimerRequest{
                unique_primers.at(uid),
                unique_mismatch_budgets.at(uid)
            }
        );
    }


    const auto unique_results =
        primer_engine_.search(
            unique_requests,
            anchor_length
        );


    if (
        unique_results.size() !=
        unique_primers.size()
    ) {
        throw std::logic_error(
            "Multiplex V2 unique-primer "
            "result-count mismatch."
        );
    }


    /*
     * ================================================
     * Intended pairs
     * ================================================
     */
    result.intended_pairs.reserve(
        requests.size()
    );


    for (
        std::size_t pair_index = 0;
        pair_index < requests.size();
        ++pair_index
    ) {
        const auto& request =
            requests.at(
                pair_index
            );


        const std::size_t uid1 =
            slot_unique_ids
                .at(pair_index)
                .at(0);


        const std::size_t uid2 =
            slot_unique_ids
                .at(pair_index)
                .at(1);


        result.intended_pairs.push_back(
            assemble_primer_pair_hits(
                request.primer1,
                unique_results.at(uid1).hits,

                request.primer2,
                unique_results.at(uid2).hits,

                request.min_amplicon_length,
                request.max_amplicon_length
            )
        );


        ++result
            .stats
            .intended_join_computations;
    }


    if (!include_cross_pairs) {
        return result;
    }


    /*
     * Logical v1 cross-slot request count.
     *
     * Kept for direct statistics comparison.
     */
    result.stats.cross_slot_pair_requests =
        (
            requests.size() *
            (
                requests.size() - 1
            ) /
            2
        ) *
        4;


    /*
     * ================================================
     * Unique primer → panel slot membership
     * ================================================
     */
    std::vector<
        std::vector<SlotMembership>
    > memberships(
        unique_primers.size()
    );


    for (
        std::size_t pair_index = 0;
        pair_index < requests.size();
        ++pair_index
    ) {
        memberships
            .at(
                slot_unique_ids
                    .at(pair_index)
                    .at(0)
            )
            .push_back(
                SlotMembership{
                    pair_index,
                    MultiplexPrimerSlot::
                        Primer1
                }
            );


        memberships
            .at(
                slot_unique_ids
                    .at(pair_index)
                    .at(1)
            )
            .push_back(
                SlotMembership{
                    pair_index,
                    MultiplexPrimerSlot::
                        Primer2
                }
            );
    }


    /*
     * ================================================
     * Global hit-level sweep
     * ================================================
     */
    std::vector<
        GlobalMultiplexPrimerHits
    > global_inputs;


    global_inputs.reserve(
        unique_primers.size()
    );


    for (
        std::size_t uid = 0;
        uid < unique_primers.size();
        ++uid
    ) {
        global_inputs.push_back(
            GlobalMultiplexPrimerHits{
                uid,
                unique_primers.at(uid).size(),
                unique_results.at(uid).hits
            }
        );
    }


    const auto global =
        global_multiplex_cross_join(
            global_inputs,
            cross_min_amplicon_length,
            cross_max_amplicon_length
        );


    result.global_cross_stats =
        global.stats;


    /*
     * ================================================
     * Physical product → panel slot memberships
     *
     * No P^2 primer-pair enumeration.
     * ================================================
     */
    for (
        const auto& product :
        global.products
    ) {
        const auto&
            forward_memberships =
                memberships.at(
                    product.forward_primer_id
                );


        const auto&
            reverse_memberships =
                memberships.at(
                    product.reverse_primer_id
                );


        for (
            const auto& forward :
            forward_memberships
        ) {
            for (
                const auto& reverse :
                reverse_memberships
            ) {
                /*
                 * Intended pair products do not belong
                 * in cross_amplicons.
                 *
                 * V1 only compares pair_a != pair_b.
                 */
                if (
                    forward.pair_index ==
                    reverse.pair_index
                ) {
                    continue;
                }


                result.cross_amplicons.push_back(
                    MultiplexCrossAmplicon{
                        forward.pair_index,
                        forward.slot,

                        reverse.pair_index,
                        reverse.slot,

                        product.amplicon_start,
                        product.amplicon_end_exclusive,
                        product.amplicon_length,

                        product.forward_mismatches,
                        product.reverse_mismatches,

                        product.forward_mismatch_mask,
                        product.reverse_mismatch_mask
                    }
                );
            }
        }
    }


    normalize_cross_products(
        result.cross_amplicons
    );


    result.stats.cross_amplicon_records =
        result.cross_amplicons.size();


    return result;
}


}  // namespace primerpair
