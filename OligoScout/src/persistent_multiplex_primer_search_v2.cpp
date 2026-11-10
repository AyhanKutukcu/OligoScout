#include "primerpair/persistent_multiplex_primer_search_v2.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "primerpair/primer_pair_assembler.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {

namespace {


struct SlotMembership {
    std::size_t pair_index{0};

    MultiplexPrimerSlot slot{
        MultiplexPrimerSlot::Primer1
    };
};


std::string canonical_primer(
    const std::string_view sequence
) {
    std::string output;

    output.reserve(
        sequence.size()
    );

    for (const char raw : sequence) {
        const char base =
            static_cast<char>(
                std::toupper(
                    static_cast<
                        unsigned char
                    >(raw)
                )
            );

        output.push_back(
            base
        );
    }

    return output;
}


bool cross_less(
    const MultiplexCrossAmplicon& lhs,
    const MultiplexCrossAmplicon& rhs
) noexcept {
    if (
        lhs.forward_pair_index !=
        rhs.forward_pair_index
    ) {
        return
            lhs.forward_pair_index <
            rhs.forward_pair_index;
    }

    if (
        lhs.forward_slot !=
        rhs.forward_slot
    ) {
        return
            static_cast<int>(
                lhs.forward_slot
            )
            <
            static_cast<int>(
                rhs.forward_slot
            );
    }

    if (
        lhs.reverse_pair_index !=
        rhs.reverse_pair_index
    ) {
        return
            lhs.reverse_pair_index <
            rhs.reverse_pair_index;
    }

    if (
        lhs.reverse_slot !=
        rhs.reverse_slot
    ) {
        return
            static_cast<int>(
                lhs.reverse_slot
            )
            <
            static_cast<int>(
                rhs.reverse_slot
            );
    }

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
        lhs.forward_mismatches !=
        rhs.forward_mismatches
    ) {
        return
            lhs.forward_mismatches <
            rhs.forward_mismatches;
    }

    if (
        lhs.reverse_mismatches !=
        rhs.reverse_mismatches
    ) {
        return
            lhs.reverse_mismatches <
            rhs.reverse_mismatches;
    }

    if (
        lhs.forward_mismatch_mask !=
        rhs.forward_mismatch_mask
    ) {
        return
            lhs.forward_mismatch_mask <
            rhs.forward_mismatch_mask;
    }

    return
        lhs.reverse_mismatch_mask <
        rhs.reverse_mismatch_mask;
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


PersistentMultiplexPrimerSearchEngineV2::
PersistentMultiplexPrimerSearchEngineV2(
    PpfmManifest manifest,
    const std::size_t cache_capacity,
    const std::size_t suffix_array_sample_rate
)
    : manifest_(
          std::move(
              manifest
          )
      ),
      cache_(
          manifest_,
          cache_capacity,
          suffix_array_sample_rate
      ) {

    if (manifest_.empty()) {
        throw std::invalid_argument(
            "Persistent multiplex manifest "
            "cannot be empty."
        );
    }
}


PersistentMultiplexSearchResultV2
PersistentMultiplexPrimerSearchEngineV2::
search(
    const std::vector<
        MultiplexPrimerPairRequest
    >& requests,

    const std::size_t anchor_length,

    const bool include_cross_pairs,

    const std::uint64_t
        cross_min_amplicon_length,

    const std::uint64_t
        cross_max_amplicon_length
) {
    PersistentMultiplexSearchResultV2 result;


    /*
     * Empty panels do not require shard loading.
     */
    if (requests.empty()) {
        return result;
    }


    if (anchor_length == 0) {
        throw std::invalid_argument(
            "Anchor length must be > 0."
        );
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
     * =================================================
     * Panel-level unique-primer table.
     *
     * This is created ONCE for the whole query, not once
     * per chromosome.
     * =================================================
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
                    sequences.at(
                        slot
                    )
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


    /*
     * Primer UID -> panel pair/slot memberships.
     *
     * Built once and reused on every chromosome.
     */
    std::vector<
        std::vector<
            SlotMembership
        >
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


    result.shards.reserve(
        manifest_.size()
    );


    /*
     * =================================================
     * CRITICAL WHOLE-GENOME LOOP
     * =================================================
     *
     * Each chromosome is fully completed before the
     * next chromosome is loaded.
     *
     * No hit list survives into another chromosome's
     * global pair sweep.
     */
    for (
        const auto& entry :
        manifest_.entries()
    ) {
        const GenomeShard& shard =
            cache_.get(
                entry.chromosome
            );


        StrandAwarePrimerSearchEngine
            primer_engine(
                shard.index(),
                shard.reference()
            );


        PersistentMultiplexShardResultV2
            shard_result;


        shard_result.shard_id =
            shard.id();


        shard_result.chromosome =
            shard.chromosome();


        shard_result.sequence_length =
            shard.sequence_length();


        shard_result.stats.pair_requests =
            requests.size();


        shard_result.stats.total_primer_slots =
            requests.size() * 2;


        shard_result.stats.unique_primer_queries =
            unique_primers.size();


        shard_result.stats.reused_primer_slots =
            shard_result
                .stats
                .total_primer_slots
            -
            unique_primers.size();


        /*
         * ---------------------------------------------
         * Search every unique primer exactly once in
         * this chromosome.
         * ---------------------------------------------
         */
        std::vector<
            StrandAwarePrimerSearchResult
        > unique_results;


        unique_results.reserve(
            unique_primers.size()
        );


        for (
            std::size_t uid = 0;
            uid < unique_primers.size();
            ++uid
        ) {
            auto primer_result =
                primer_engine.search(
                    unique_primers.at(uid),
                    anchor_length,
                    unique_mismatch_budgets
                        .at(uid)
                );


            for (
                const auto& hit :
                primer_result.hits
            ) {
                ++shard_result
                    .stats
                    .total_primer_hits;


                if (
                    hit.orientation ==
                    PrimerOrientation::
                        Forward
                ) {
                    ++shard_result
                        .stats
                        .forward_primer_hits;

                } else {
                    ++shard_result
                        .stats
                        .reverse_primer_hits;
                }
            }


            unique_results.push_back(
                std::move(
                    primer_result
                )
            );
        }


        /*
         * ---------------------------------------------
         * Intended pairs.
         *
         * Assembly is chromosome-local.
         * ---------------------------------------------
         */
        shard_result.intended_pairs.reserve(
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


            shard_result
                .intended_pairs
                .push_back(
                    assemble_primer_pair_hits(
                        request.primer1,
                        unique_results
                            .at(uid1)
                            .hits,

                        request.primer2,
                        unique_results
                            .at(uid2)
                            .hits,

                        request
                            .min_amplicon_length,

                        request
                            .max_amplicon_length
                    )
                );


            ++shard_result
                .stats
                .intended_join_computations;
        }


        /*
         * Cross-pair analysis optional.
         */
        if (include_cross_pairs) {

            shard_result
                .stats
                .logical_cross_slot_requests =
                    (
                        requests.size() *
                        (
                            requests.size() - 1
                        )
                        /
                        2
                    )
                    *
                    4;


            /*
             * -----------------------------------------
             * ONE global hit-level sweep for this
             * chromosome.
             * -----------------------------------------
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
                        unique_primers
                            .at(uid)
                            .size(),

                        unique_results
                            .at(uid)
                            .hits
                    }
                );
            }


            const auto global =
                global_multiplex_cross_join(
                    global_inputs,
                    cross_min_amplicon_length,
                    cross_max_amplicon_length
                );


            shard_result.global_cross_stats =
                global.stats;


            /*
             * -----------------------------------------
             * Physical product -> panel memberships.
             *
             * No pair_a/pair_b O(P^2) enumeration.
             * -----------------------------------------
             */
            for (
                const auto& product :
                global.products
            ) {
                const auto&
                    forward_memberships =
                        memberships.at(
                            product
                                .forward_primer_id
                        );


                const auto&
                    reverse_memberships =
                        memberships.at(
                            product
                                .reverse_primer_id
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
                         * Product belonging to the same
                         * intended pair is not a multiplex
                         * cross-amplicon.
                         */
                        if (
                            forward.pair_index ==
                            reverse.pair_index
                        ) {
                            continue;
                        }


                        shard_result
                            .cross_amplicons
                            .push_back(
                                MultiplexCrossAmplicon{
                                    forward
                                        .pair_index,

                                    forward
                                        .slot,

                                    reverse
                                        .pair_index,

                                    reverse
                                        .slot,

                                    product
                                        .amplicon_start,

                                    product
                                        .amplicon_end_exclusive,

                                    product
                                        .amplicon_length,

                                    product
                                        .forward_mismatches,

                                    product
                                        .reverse_mismatches,

                                    product
                                        .forward_mismatch_mask,

                                    product
                                        .reverse_mismatch_mask
                                }
                            );
                    }
                }
            }


            normalize_cross_products(
                shard_result
                    .cross_amplicons
            );


            shard_result
                .stats
                .cross_amplicon_records =
                    shard_result
                        .cross_amplicons
                        .size();
        }


        result.shards.push_back(
            std::move(
                shard_result
            )
        );
    }


    return result;
}



std::vector<
    PersistentMultiplexSearchResultV2
>
PersistentMultiplexPrimerSearchEngineV2::
search_many(
    const std::vector<
        std::vector<
            MultiplexPrimerPairRequest
        >
    >& panels,

    const std::size_t anchor_length,

    const bool include_cross_pairs,

    const std::uint64_t
        cross_min_amplicon_length,

    const std::uint64_t
        cross_max_amplicon_length
) {
    if (panels.empty()) {
        return {};
    }


    if (anchor_length == 0) {
        throw std::invalid_argument(
            "Persistent multiplex batch "
            "anchor length must be > 0."
        );
    }


    if (
        cross_min_amplicon_length == 0
        ||
        cross_min_amplicon_length >
            cross_max_amplicon_length
    ) {
        throw std::invalid_argument(
            "Invalid persistent multiplex "
            "batch cross-amplicon range."
        );
    }


    struct SlotMembership {
        std::size_t pair_index{0};

        MultiplexPrimerSlot slot{
            MultiplexPrimerSlot::Primer1
        };
    };


    struct PreparedPanel {
        bool active{false};

        std::vector<
            std::string
        > unique_primers{};

        std::vector<
            std::size_t
        > unique_mismatch_budgets{};

        std::vector<
            std::array<
                std::size_t,
                2
            >
        > slot_unique_ids{};

        std::vector<
            std::vector<
                SlotMembership
            >
        > memberships{};
    };


    const auto canonicalize =
        [](
            const std::string_view primer
        ) {
            if (primer.empty()) {
                throw std::invalid_argument(
                    "Persistent multiplex primer "
                    "cannot be empty."
                );
            }


            std::string result;

            result.reserve(
                primer.size()
            );


            for (const char raw : primer) {

                const char base =
                    static_cast<char>(
                        std::toupper(
                            static_cast<
                                unsigned char
                            >(
                                raw
                            )
                        )
                    );


                if (
                    base != 'A'
                    &&
                    base != 'C'
                    &&
                    base != 'G'
                    &&
                    base != 'T'
                ) {
                    throw std::invalid_argument(
                        "Persistent multiplex primer "
                        "must contain only A/C/G/T."
                    );
                }


                result.push_back(
                    base
                );
            }


            return result;
        };


    /*
     * -------------------------------------------------
     * Prepare panel metadata ONCE.
     *
     * This work is shared by all chromosome shards.
     * -------------------------------------------------
     */
    std::vector<
        PreparedPanel
    > prepared_panels;


    prepared_panels.resize(
        panels.size()
    );


    for (
        std::size_t panel_index = 0;
        panel_index < panels.size();
        ++panel_index
    ) {
        const auto& requests =
            panels.at(
                panel_index
            );


        if (requests.empty()) {
            continue;
        }


        auto& prepared =
            prepared_panels.at(
                panel_index
            );


        prepared.active = true;


        prepared.slot_unique_ids.resize(
            requests.size()
        );


        std::map<
            std::pair<
                std::string,
                std::size_t
            >,
            std::size_t
        > unique_ids;


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
                std::pair<
                    std::string_view,
                    MultiplexPrimerSlot
                >,
                2
            > slots{
                std::pair{
                    request.primer1,
                    MultiplexPrimerSlot::
                        Primer1
                },

                std::pair{
                    request.primer2,
                    MultiplexPrimerSlot::
                        Primer2
                }
            };


            for (
                std::size_t slot_index = 0;
                slot_index < slots.size();
                ++slot_index
            ) {
                std::string canonical =
                    canonicalize(
                        slots
                            .at(slot_index)
                            .first
                    );


                const auto key =
                    std::make_pair(
                        canonical,
                        request
                            .max_mismatches
                    );


                auto found =
                    unique_ids.find(
                        key
                    );


                std::size_t uid = 0;


                if (
                    found ==
                    unique_ids.end()
                ) {
                    uid =
                        prepared
                            .unique_primers
                            .size();


                    unique_ids.emplace(
                        key,
                        uid
                    );


                    prepared
                        .unique_primers
                        .push_back(
                            std::move(
                                canonical
                            )
                        );


                    prepared
                        .unique_mismatch_budgets
                        .push_back(
                            request
                                .max_mismatches
                        );


                    prepared
                        .memberships
                        .emplace_back();

                } else {
                    uid =
                        found->second;
                }


                prepared
                    .slot_unique_ids
                    .at(pair_index)
                    .at(slot_index) =
                        uid;


                prepared
                    .memberships
                    .at(uid)
                    .push_back(
                        SlotMembership{
                            pair_index,
                            slots
                                .at(slot_index)
                                .second
                        }
                    );
            }
        }
    }


    /*
     * One final result object per input panel.
     */
    std::vector<
        PersistentMultiplexSearchResultV2
    > outputs(
        panels.size()
    );


    /*
     * =================================================
     * SHARD-MAJOR EXECUTION.
     *
     * THIS IS THE IMPORTANT DIFFERENCE.
     * =================================================
     */
    for (
        const auto& entry :
        manifest_.entries()
    ) {
        const GenomeShard& shard =
            cache_.get(
                entry.chromosome
            );


        StrandAwarePrimerSearchEngine
            primer_engine(
                shard.index(),
                shard.reference()
            );


        /*
         * Every active panel is completed while this
         * chromosome remains resident.
         */
        for (
            std::size_t panel_index = 0;
            panel_index < panels.size();
            ++panel_index
        ) {
            const auto& requests =
                panels.at(
                    panel_index
                );


            const auto& prepared =
                prepared_panels.at(
                    panel_index
                );


            /*
             * Preserve search(empty-panel) behavior:
             * empty input returns an empty result rather
             * than 24 empty shard records.
             */
            if (!prepared.active) {
                continue;
            }


            PersistentMultiplexShardResultV2
                shard_result;


            shard_result.shard_id =
                shard.id();


            shard_result.chromosome =
                shard.chromosome();


            shard_result.sequence_length =
                shard.sequence_length();


            shard_result.stats.pair_requests =
                requests.size();


            shard_result.stats.total_primer_slots =
                requests.size() * 2;


            shard_result
                .stats
                .unique_primer_queries =
                    prepared
                        .unique_primers
                        .size();


            shard_result
                .stats
                .reused_primer_slots =
                    shard_result
                        .stats
                        .total_primer_slots
                    -
                    prepared
                        .unique_primers
                        .size();


            /*
             * -----------------------------------------
             * Search unique primers once per panel for
             * this chromosome.
             * -----------------------------------------
             */
            std::vector<
                StrandAwarePrimerSearchResult
            > unique_results;


            unique_results.reserve(
                prepared
                    .unique_primers
                    .size()
            );


            for (
                std::size_t uid = 0;
                uid <
                    prepared
                        .unique_primers
                        .size();
                ++uid
            ) {
                auto primer_result =
                    primer_engine.search(
                        prepared
                            .unique_primers
                            .at(uid),

                        anchor_length,

                        prepared
                            .unique_mismatch_budgets
                            .at(uid)
                    );


                for (
                    const auto& hit :
                    primer_result.hits
                ) {
                    ++shard_result
                        .stats
                        .total_primer_hits;


                    if (
                        hit.orientation ==
                        PrimerOrientation::
                            Forward
                    ) {
                        ++shard_result
                            .stats
                            .forward_primer_hits;

                    } else {
                        ++shard_result
                            .stats
                            .reverse_primer_hits;
                    }
                }


                unique_results.push_back(
                    std::move(
                        primer_result
                    )
                );
            }


            /*
             * -----------------------------------------
             * Intended pair assembly.
             * -----------------------------------------
             */
            shard_result
                .intended_pairs
                .reserve(
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
                    prepared
                        .slot_unique_ids
                        .at(pair_index)
                        .at(0);


                const std::size_t uid2 =
                    prepared
                        .slot_unique_ids
                        .at(pair_index)
                        .at(1);


                shard_result
                    .intended_pairs
                    .push_back(
                        assemble_primer_pair_hits(
                            request.primer1,

                            unique_results
                                .at(uid1)
                                .hits,

                            request.primer2,

                            unique_results
                                .at(uid2)
                                .hits,

                            request
                                .min_amplicon_length,

                            request
                                .max_amplicon_length
                        )
                    );


                ++shard_result
                    .stats
                    .intended_join_computations;
            }


            /*
             * -----------------------------------------
             * Global multiplex cross sweep.
             * -----------------------------------------
             */
            if (include_cross_pairs) {

                shard_result
                    .stats
                    .logical_cross_slot_requests =
                        (
                            requests.size()
                            *
                            (
                                requests.size()
                                -
                                1
                            )
                            /
                            2
                        )
                        *
                        4;


                std::vector<
                    GlobalMultiplexPrimerHits
                > global_inputs;


                global_inputs.reserve(
                    prepared
                        .unique_primers
                        .size()
                );


                for (
                    std::size_t uid = 0;
                    uid <
                        prepared
                            .unique_primers
                            .size();
                    ++uid
                ) {
                    global_inputs.push_back(
                        GlobalMultiplexPrimerHits{
                            uid,

                            prepared
                                .unique_primers
                                .at(uid)
                                .size(),

                            unique_results
                                .at(uid)
                                .hits
                        }
                    );
                }


                const auto global =
                    global_multiplex_cross_join(
                        global_inputs,

                        cross_min_amplicon_length,

                        cross_max_amplicon_length
                    );


                shard_result.global_cross_stats =
                    global.stats;


                for (
                    const auto& product :
                    global.products
                ) {
                    const auto&
                        forward_memberships =
                            prepared
                                .memberships
                                .at(
                                    product
                                        .forward_primer_id
                                );


                    const auto&
                        reverse_memberships =
                            prepared
                                .memberships
                                .at(
                                    product
                                        .reverse_primer_id
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
                             * Same intended pair is
                             * not a multiplex cross
                             * product.
                             */
                            if (
                                forward.pair_index ==
                                reverse.pair_index
                            ) {
                                continue;
                            }


                            shard_result
                                .cross_amplicons
                                .push_back(
                                    MultiplexCrossAmplicon{
                                        forward
                                            .pair_index,

                                        forward
                                            .slot,

                                        reverse
                                            .pair_index,

                                        reverse
                                            .slot,

                                        product
                                            .amplicon_start,

                                        product
                                            .amplicon_end_exclusive,

                                        product
                                            .amplicon_length,

                                        product
                                            .forward_mismatches,

                                        product
                                            .reverse_mismatches,

                                        product
                                            .forward_mismatch_mask,

                                        product
                                            .reverse_mismatch_mask
                                    }
                                );
                        }
                    }
                }


                normalize_cross_products(
                    shard_result
                        .cross_amplicons
                );


                shard_result
                    .stats
                    .cross_amplicon_records =
                        shard_result
                            .cross_amplicons
                            .size();
            }


            outputs
                .at(panel_index)
                .shards
                .push_back(
                    std::move(
                        shard_result
                    )
                );
        }
    }


    return outputs;
}


}  // namespace primerpair
