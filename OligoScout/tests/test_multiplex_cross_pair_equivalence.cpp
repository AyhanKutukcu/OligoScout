#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/multiplex_primer_search.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {


std::string make_reference(
    const std::size_t length
) {
    constexpr char bases[] = {
        'A', 'C', 'G', 'T'
    };

    std::uint64_t state =
        0xD1B54A32D192ED03ULL;

    std::string reference;
    reference.reserve(length);

    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        reference.push_back(
            bases[
                static_cast<std::size_t>(
                    (state >> 32U) &
                    3ULL
                )
            ]
        );
    }

    return reference;
}


bool cross_less(
    const primerpair::MultiplexCrossAmplicon& lhs,
    const primerpair::MultiplexCrossAmplicon& rhs
) {
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


void normalize(
    std::vector<
        primerpair::MultiplexCrossAmplicon
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


int main() {
    try {
        using namespace primerpair;

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            primer_length = 24;

        constexpr std::size_t
            pair_count = 8;

        constexpr std::uint64_t
            cross_min = 50;

        constexpr std::uint64_t
            cross_max = 1500;


        const std::string reference =
            make_reference(
                120000
            );


        PackedReference packed(
            reference
        );

        BidirectionalFMIndex bifm(
            reference
        );

        IPBWTIndex ipbwt(
            reference,
            anchor_length
        );

        SearchDifficultyEstimator estimator(
            bifm
        );

        BatchedAnchorLookup anchor_lookup(
            ipbwt,
            estimator
        );

        AnchorCandidateSearcher verifier(
            bifm,
            packed
        );

        BatchedCandidateSearchEngine
            candidate_engine(
                ipbwt,
                anchor_lookup,
                verifier
            );

        SinglePrimerSearchEngine
            single_engine(
                bifm,
                packed
            );

        HybridBatchedPrimerSearchEngine
            forward_hybrid(
                candidate_engine,
                single_engine
            );

        StrandAwarePrimerSearchEngine
            legacy_strand(
                bifm,
                packed
            );

        HybridStrandAwarePrimerSearchEngine
            hybrid_strand(
                forward_hybrid,
                anchor_lookup,
                ipbwt,
                packed,
                legacy_strand
            );

        MultiplexPrimerSearchEngine
            multiplex(
                hybrid_strand
            );

        PrimerPairSearchEngine
            legacy_pair(
                bifm,
                packed
            );


        std::vector<std::string>
            primer1_storage;

        std::vector<std::string>
            primer2_storage;

        primer1_storage.reserve(
            pair_count
        );

        primer2_storage.reserve(
            pair_count
        );


        /*
         * Aynı ~2kb bölge içinde farklı pair'ler:
         * çok sayıda cross-product oluşturması
         * beklenir.
         */
        for (
            std::size_t i = 0;
            i < pair_count;
            ++i
        ) {
            const std::size_t left =
                10000 +
                i * 110;

            const std::size_t right =
                left +
                300 +
                i * 17;


            primer1_storage.push_back(
                reference.substr(
                    left,
                    primer_length
                )
            );


            primer2_storage.push_back(
                reverse_complement(
                    reference.substr(
                        right,
                        primer_length
                    )
                )
            );
        }


        std::vector<
            MultiplexPrimerPairRequest
        > requests;

        requests.reserve(
            pair_count
        );


        for (
            std::size_t i = 0;
            i < pair_count;
            ++i
        ) {
            requests.push_back(
                MultiplexPrimerPairRequest{
                    primer1_storage.at(i),
                    primer2_storage.at(i),
                    3,
                    50,
                    1500
                }
            );
        }


        const auto observed =
            multiplex.search(
                requests,
                anchor_length,
                true,
                cross_min,
                cross_max
            );


        std::vector<
            MultiplexCrossAmplicon
        > expected;


        for (
            std::size_t pair_a = 0;
            pair_a < pair_count;
            ++pair_a
        ) {
            for (
                std::size_t pair_b =
                    pair_a + 1;

                pair_b < pair_count;
                ++pair_b
            ) {
                const std::array<
                    std::string_view,
                    2
                > seq_a{
                    primer1_storage.at(pair_a),
                    primer2_storage.at(pair_a)
                };


                const std::array<
                    std::string_view,
                    2
                > seq_b{
                    primer1_storage.at(pair_b),
                    primer2_storage.at(pair_b)
                };


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
                        const auto products =
                            legacy_pair.search(
                                seq_a.at(slot_a),
                                seq_b.at(slot_b),

                                anchor_length,
                                3,

                                cross_min,
                                cross_max
                            );


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
                            products.amplicons
                        ) {
                            const bool
                                a_is_forward =
                                    product.left_primer ==
                                    PrimerIdentity::
                                        Primer1;


                            expected.push_back(
                                MultiplexCrossAmplicon{
                                    a_is_forward
                                        ? pair_a
                                        : pair_b,

                                    a_is_forward
                                        ? slot_enum_a
                                        : slot_enum_b,

                                    a_is_forward
                                        ? pair_b
                                        : pair_a,

                                    a_is_forward
                                        ? slot_enum_b
                                        : slot_enum_a,

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


        auto actual =
            observed.cross_amplicons;


        normalize(
            expected
        );

        normalize(
            actual
        );


        if (
            expected !=
            actual
        ) {
            std::cerr
                << "CROSS_PAIR_MISMATCH"
                << '\t'
                << "expected="
                << expected.size()
                << '\t'
                << "actual="
                << actual.size()
                << '\n';

            return 1;
        }


        const std::size_t
            expected_slot_requests =
                (
                    pair_count *
                    (
                        pair_count - 1
                    ) /
                    2
                ) *
                4;


        if (
            observed
                .stats
                .cross_slot_pair_requests
            !=
            expected_slot_requests
        ) {
            throw std::runtime_error(
                "Cross slot-pair request count mismatch."
            );
        }


        std::cout
            << "pairs\t"
            << pair_count
            << '\n';

        std::cout
            << "cross_slot_pair_requests\t"
            << observed
                .stats
                .cross_slot_pair_requests
            << '\n';

        std::cout
            << "unique_cross_join_computations\t"
            << observed
                .stats
                .unique_cross_join_computations
            << '\n';

        std::cout
            << "cross_amplicons\t"
            << actual.size()
            << '\n';

        std::cout
            << "ALL_CHECKS\tYES\n";

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';

        return 1;
    }
}
