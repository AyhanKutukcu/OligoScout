#include "primerpair/multiplex_primer_search.hpp"

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


MultiplexPrimerSearchEngine::
MultiplexPrimerSearchEngine(
    const HybridStrandAwarePrimerSearchEngine&
        primer_engine
) noexcept
    : primer_engine_(
          primer_engine
      ) {
}


MultiplexPrimerSearchResult
MultiplexPrimerSearchEngine::search(
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
    MultiplexPrimerSearchResult result;

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


    /*
     * ------------------------------------------------
     * Unique primer table.
     *
     * Key:
     *   (canonical uppercase sequence,
     *    max_mismatches)
     * ------------------------------------------------
     */

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


    unique_primers.reserve(
        requests.size() * 2
    );

    unique_mismatch_budgets.reserve(
        requests.size() * 2
    );


    std::vector<
        std::array<std::size_t, 2>
    > slot_unique_ids(
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
            const std::string canonical =
                canonical_primer(
                    sequences.at(
                        slot
                    )
                );


            const PrimerKey key{
                canonical,
                request.max_mismatches
            };


            const auto existing =
                unique_lookup.find(
                    key
                );


            if (
                existing !=
                unique_lookup.end()
            ) {
                slot_unique_ids
                    .at(pair_index)
                    .at(slot) =
                        existing->second;

                continue;
            }


            const std::size_t id =
                unique_primers.size();


            unique_primers.push_back(
                canonical
            );


            unique_mismatch_budgets
                .push_back(
                    request.max_mismatches
                );


            unique_lookup.emplace(
                key,
                id
            );


            slot_unique_ids
                .at(pair_index)
                .at(slot) =
                    id;
        }
    }


    result.stats.unique_primer_queries =
        unique_primers.size();


    result.stats.reused_primer_slots =
        result.stats.total_primer_slots -
        result.stats.unique_primer_queries;


    /*
     * ------------------------------------------------
     * One both-strand hybrid search for ALL unique
     * primers.
     * ------------------------------------------------
     */

    std::vector<
        HybridStrandAwarePrimerRequest
    > unique_requests;


    unique_requests.reserve(
        unique_primers.size()
    );


    for (
        std::size_t id = 0;
        id < unique_primers.size();
        ++id
    ) {
        unique_requests.push_back(
            HybridStrandAwarePrimerRequest{
                unique_primers.at(
                    id
                ),
                unique_mismatch_budgets.at(
                    id
                )
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
            "Multiplex unique-primer "
            "result-count mismatch."
        );
    }


    /*
     * ------------------------------------------------
     * Intended pair assembly.
     *
     * Search tekrar edilmez. Yalnız hit listeleri
     * reusable sweep assembler'a verilir.
     * ------------------------------------------------
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


        const std::size_t id1 =
            slot_unique_ids
                .at(pair_index)
                .at(0);


        const std::size_t id2 =
            slot_unique_ids
                .at(pair_index)
                .at(1);


        result.intended_pairs.push_back(
            assemble_primer_pair_hits(
                request.primer1,
                unique_results
                    .at(id1)
                    .hits,

                request.primer2,
                unique_results
                    .at(id2)
                    .hits,

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
     * ------------------------------------------------
     * Cross-pair unique join cache.
     *
     * Key primer IDs are canonicalized:
     *
     *   low_id <= high_id
     *
     * Therefore the same physical primer combination
     * is sweep-joined only once.
     * ------------------------------------------------
     */

    using CrossJoinKey =
        std::tuple<
            std::size_t,
            std::size_t,
            std::uint64_t,
            std::uint64_t
        >;


    std::map<
        CrossJoinKey,
        std::vector<PrimerPairHit>
    > cross_join_cache;


    /*
     * Pair_i x Pair_j, i < j.
     *
     * Four primer-slot combinations:
     *
     *   P1_i / P1_j
     *   P1_i / P2_j
     *   P2_i / P1_j
     *   P2_i / P2_j
     */
    for (
        std::size_t pair_a = 0;
        pair_a < requests.size();
        ++pair_a
    ) {
        for (
            std::size_t pair_b =
                pair_a + 1;

            pair_b < requests.size();
            ++pair_b
        ) {
            for (
                std::size_t slot_a = 0;
                slot_a < 2;
                ++slot_a
            ) {
                for (
                    std::size_t slot_b = 0;
                    slot_b < 2;
                    ++slot_b
                ) {
                    ++result
                        .stats
                        .cross_slot_pair_requests;


                    const std::size_t uid_a =
                        slot_unique_ids
                            .at(pair_a)
                            .at(slot_a);


                    const std::size_t uid_b =
                        slot_unique_ids
                            .at(pair_b)
                            .at(slot_b);


                    const std::size_t low_uid =
                        std::min(
                            uid_a,
                            uid_b
                        );


                    const std::size_t high_uid =
                        std::max(
                            uid_a,
                            uid_b
                        );


                    const CrossJoinKey key{
                        low_uid,
                        high_uid,
                        cross_min_amplicon_length,
                        cross_max_amplicon_length
                    };


                    auto cached =
                        cross_join_cache.find(
                            key
                        );


                    if (
                        cached ==
                        cross_join_cache.end()
                    ) {
                        auto assembled =
                            assemble_primer_pair_hits(
                                unique_primers.at(
                                    low_uid
                                ),

                                unique_results
                                    .at(low_uid)
                                    .hits,

                                unique_primers.at(
                                    high_uid
                                ),

                                unique_results
                                    .at(high_uid)
                                    .hits,

                                cross_min_amplicon_length,
                                cross_max_amplicon_length
                            );


                        cached =
                            cross_join_cache
                                .emplace(
                                    key,
                                    std::move(
                                        assembled
                                            .amplicons
                                    )
                                )
                                .first;


                        ++result
                            .stats
                            .unique_cross_join_computations;

                    } else {
                        ++result
                            .stats
                            .reused_cross_join_requests;
                    }


                    const auto slot_enum_a =
                        slot_a == 0
                            ? MultiplexPrimerSlot::
                                Primer1
                            : MultiplexPrimerSlot::
                                Primer2;


                    const auto slot_enum_b =
                        slot_b == 0
                            ? MultiplexPrimerSlot::
                                Primer1
                            : MultiplexPrimerSlot::
                                Primer2;


                    for (
                        const auto& product :
                        cached->second
                    ) {
                        std::size_t
                            forward_pair = 0;

                        MultiplexPrimerSlot
                            forward_slot =
                                MultiplexPrimerSlot::
                                    Primer1;

                        std::size_t
                            reverse_pair = 0;

                        MultiplexPrimerSlot
                            reverse_slot =
                                MultiplexPrimerSlot::
                                    Primer1;


                        if (
                            uid_a !=
                            uid_b
                        ) {
                            const std::size_t
                                forward_uid =
                                    product.left_primer ==
                                        PrimerIdentity::
                                            Primer1
                                        ? low_uid
                                        : high_uid;


                            if (
                                forward_uid ==
                                uid_a
                            ) {
                                forward_pair =
                                    pair_a;

                                forward_slot =
                                    slot_enum_a;

                                reverse_pair =
                                    pair_b;

                                reverse_slot =
                                    slot_enum_b;

                            } else {
                                forward_pair =
                                    pair_b;

                                forward_slot =
                                    slot_enum_b;

                                reverse_pair =
                                    pair_a;

                                reverse_slot =
                                    slot_enum_a;
                            }

                        } else {
                            /*
                             * Aynı sequence iki farklı panel slotunda
                             * kullanılıyorsa local Primer1/Primer2
                             * kimliği slot association için kullanılır.
                             */
                            if (
                                product.left_primer ==
                                PrimerIdentity::
                                    Primer1
                            ) {
                                forward_pair =
                                    pair_a;

                                forward_slot =
                                    slot_enum_a;

                                reverse_pair =
                                    pair_b;

                                reverse_slot =
                                    slot_enum_b;

                            } else {
                                forward_pair =
                                    pair_b;

                                forward_slot =
                                    slot_enum_b;

                                reverse_pair =
                                    pair_a;

                                reverse_slot =
                                    slot_enum_a;
                            }
                        }


                        result.cross_amplicons.push_back(
                            MultiplexCrossAmplicon{
                                forward_pair,
                                forward_slot,

                                reverse_pair,
                                reverse_slot,

                                product
                                    .amplicon_start,

                                product
                                    .amplicon_end_exclusive,

                                product
                                    .amplicon_length,

                                product
                                    .left_mismatches,

                                product
                                    .right_mismatches,

                                product
                                    .left_mismatch_mask,

                                product
                                    .right_mismatch_mask
                            }
                        );
                    }
                }
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
