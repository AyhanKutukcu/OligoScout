#include "primerpair/primer_pair_search.hpp"

#include <algorithm>

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect(
    const bool condition,
    const std::string& name
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " + name
        );
    }

    std::cout
        << "[PASS] "
        << name
        << '\n';
}

}  // namespace

int main() {
    try {
        const std::string primer1 =
            "ACGTTGCAAGTCCTGAACGA";

        const std::string primer2 =
            "TGCAGATCCGTACGATGTCA";

        expect(
            primer1.size() == 20,
            "Primer1 length"
        );

        expect(
            primer2.size() == 20,
            "Primer2 length"
        );

        const std::string rc1 =
            primerpair::reverse_complement(
                primer1
            );

        const std::string rc2 =
            primerpair::reverse_complement(
                primer2
            );

        std::string reference;

        /*
         * --------------------------------------------------
         * Product 1
         *
         * Primer1 --->       <--- Primer2
         *
         * Expected length = 100 bp.
         * --------------------------------------------------
         */

        reference.append(
            10,
            'N'
        );

        const std::uint64_t p1_forward_pos =
            reference.size();

        reference +=
            primer1;

        reference.append(
            60,
            'N'
        );

        const std::uint64_t p2_reverse_pos =
            reference.size();

        reference +=
            rc2;

        const std::uint64_t
            product1_length =
                p2_reverse_pos +
                primer2.size() -
                p1_forward_pos;

        expect(
            product1_length == 100,
            "Synthetic product1 length"
        );

        /*
         * Large separator prevents accidental
         * cross-pairing inside the 50-120 bp window.
         */
        reference.append(
            200,
            'N'
        );

        /*
         * --------------------------------------------------
         * Product 2
         *
         * Primer2 --->       <--- Primer1
         *
         * Expected length = 90 bp.
         * --------------------------------------------------
         */

        const std::uint64_t p2_forward_pos =
            reference.size();

        reference +=
            primer2;

        reference.append(
            50,
            'N'
        );

        const std::uint64_t p1_reverse_pos =
            reference.size();

        reference +=
            rc1;

        const std::uint64_t
            product2_length =
                p1_reverse_pos +
                primer1.size() -
                p2_forward_pos;

        expect(
            product2_length == 90,
            "Synthetic product2 length"
        );

        /*
         * --------------------------------------------------
         * Outward / distant decoy geometry.
         *
         * rc1 comes before primer2.
         *
         * It must not create a valid product inside
         * the 50-120 bp search window.
         * --------------------------------------------------
         */

        reference.append(
            200,
            'N'
        );

        reference +=
            rc1;

        reference.append(
            60,
            'N'
        );

        reference +=
            primer2;

        reference.append(
            20,
            'N'
        );

        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::PrimerPairSearchEngine
            engine(
                index,
                packed
            );

        /*
         * Exact matching only.
         *
         * Amplicon window 50-120 bp should retain
         * exactly the two designed products.
         */
        const auto result =
            engine.search(
                primer1,
                primer2,
                12,
                0,
                50,
                120
            );

        expect(
            result.amplicon_count() == 2,
            "Two valid exact amplicons"
        );

        expect(
            result.amplicons.at(0) ==
                primerpair::PrimerPairHit{
                    primerpair::
                        PrimerIdentity::
                        Primer1,

                    primerpair::
                        PrimerIdentity::
                        Primer2,

                    p1_forward_pos,
                    p2_reverse_pos,

                    0,
                    0,

                    p1_forward_pos,
                    p2_reverse_pos +
                        primer2.size(),

                    100
                },
            "Primer1-forward Primer2-reverse product"
        );

        expect(
            result.amplicons.at(1) ==
                primerpair::PrimerPairHit{
                    primerpair::
                        PrimerIdentity::
                        Primer2,

                    primerpair::
                        PrimerIdentity::
                        Primer1,

                    p2_forward_pos,
                    p1_reverse_pos,

                    0,
                    0,

                    p2_forward_pos,
                    p1_reverse_pos +
                        primer1.size(),

                    90
                },
            "Primer2-forward Primer1-reverse product"
        );

        /*
         * Tight size window:
         * only the 100 bp product remains.
         */
        const auto tight =
            engine.search(
                primer1,
                primer2,
                12,
                0,
                95,
                105
            );

        expect(
            tight.amplicon_count() == 1,
            "Amplicon size filtering"
        );

        expect(
            tight.amplicons.front()
                .amplicon_length == 100,
            "100 bp product retained"
        );

        /*
         * 101-120 excludes both 90 and 100 bp.
         */
        const auto none =
            engine.search(
                primer1,
                primer2,
                12,
                0,
                101,
                120
            );

        expect(
            none.empty(),
            "Out-of-range products rejected"
        );

        expect(
            result.amplicons.at(0)
                .total_mismatches() == 0,
            "Pair mismatch accounting"
        );

        expect(
            std::string(
                primerpair::to_string(
                    primerpair::
                        PrimerIdentity::
                        Primer1
                )
            ) ==
                "PRIMER1",
            "Primer1 identity name"
        );

        expect(
            std::string(
                primerpair::to_string(
                    primerpair::
                        PrimerIdentity::
                        Primer2
                )
            ) ==
                "PRIMER2",
            "Primer2 identity name"
        );

        /*
         * --------------------------------------------------
         * Sorted-stream merge test
         *
         * Put the Primer2->Primer1 product BEFORE the
         * Primer1->Primer2 product genomically.
         *
         * Concatenating the two orientation streams would
         * therefore produce the wrong output order.
         * A correct std::merge must restore genomic order.
         * --------------------------------------------------
         */

        std::string merge_reference(
            10,
            'N'
        );

        const std::uint64_t
            early_p2_forward =
                merge_reference.size();

        merge_reference +=
            primer2;

        merge_reference.append(
            50,
            'N'
        );

        const std::uint64_t
            early_p1_reverse =
                merge_reference.size();

        merge_reference +=
            rc1;

        merge_reference.append(
            200,
            'N'
        );

        const std::uint64_t
            late_p1_forward =
                merge_reference.size();

        merge_reference +=
            primer1;

        merge_reference.append(
            60,
            'N'
        );

        const std::uint64_t
            late_p2_reverse =
                merge_reference.size();

        merge_reference +=
            rc2;

        merge_reference.append(
            10,
            'N'
        );

        const primerpair::PackedReference
            merge_packed(
                merge_reference
            );

        const primerpair::BidirectionalFMIndex
            merge_index(
                merge_reference
            );

        const primerpair::PrimerPairSearchEngine
            merge_engine(
                merge_index,
                merge_packed
            );

        const auto merge_result =
            merge_engine.search(
                primer1,
                primer2,
                12,
                0,
                50,
                120
            );

        expect(
            merge_result.amplicon_count() == 2,
            "Two interleaved orientation products"
        );

        expect(
            merge_result.amplicons.at(0)
                .left_primer ==
                primerpair::
                    PrimerIdentity::
                    Primer2,
            "Earlier Primer2-forward product merged first"
        );

        expect(
            merge_result.amplicons.at(0)
                .amplicon_start ==
                early_p2_forward,
            "Earlier merged product coordinate"
        );

        expect(
            merge_result.amplicons.at(1)
                .left_primer ==
                primerpair::
                    PrimerIdentity::
                    Primer1,
            "Later Primer1-forward product merged second"
        );

        expect(
            merge_result.amplicons.at(1)
                .amplicon_start ==
                late_p1_forward,
            "Later merged product coordinate"
        );

        expect(
            early_p1_reverse >
                early_p2_forward &&
            late_p2_reverse >
                late_p1_forward,
            "Interleaved synthetic geometry valid"
        );


        /*
         * --------------------------------------------------
         * Pair-level mismatch mask propagation
         *
         * Primer1 has one mismatch at ORIGINAL
         * primer position 0 (5-prime-most base).
         *
         * Primer2 reverse binding site is exact.
         * --------------------------------------------------
         */

        std::string primer1_5prime_mismatch =
            primer1;

        primer1_5prime_mismatch.at(0) =
            primer1_5prime_mismatch.at(0) == 'A'
                ? 'C'
                : 'A';

        std::string mask_reference(
            10,
            'N'
        );

        const std::uint64_t
            mask_left_position =
                mask_reference.size();

        mask_reference +=
            primer1_5prime_mismatch;

        mask_reference.append(
            60,
            'N'
        );

        const std::uint64_t
            mask_right_position =
                mask_reference.size();

        mask_reference +=
            rc2;

        mask_reference.append(
            10,
            'N'
        );

        const primerpair::PackedReference
            mask_packed(
                mask_reference
            );

        const primerpair::BidirectionalFMIndex
            mask_index(
                mask_reference
            );

        const primerpair::PrimerPairSearchEngine
            mask_engine(
                mask_index,
                mask_packed
            );

        const auto mask_result =
            mask_engine.search(
                primer1,
                primer2,
                12,
                1,
                50,
                120
            );

        const auto mask_product_it =
            std::find_if(
                mask_result.amplicons.begin(),
                mask_result.amplicons.end(),
                [
                    mask_left_position,
                    mask_right_position
                ](
                    const auto& product
                ) {
                    return
                        product.left_primer ==
                            primerpair::
                                PrimerIdentity::
                                Primer1 &&
                        product.right_primer ==
                            primerpair::
                                PrimerIdentity::
                                Primer2 &&
                        product.left_position ==
                            mask_left_position &&
                        product.right_position ==
                            mask_right_position;
                }
            );

        expect(
            mask_product_it !=
                mask_result.amplicons.end(),
            "Pair-level mismatch product found"
        );

        expect(
            mask_product_it
                ->left_mismatches == 1,
            "Pair-level left mismatch count"
        );

        expect(
            mask_product_it
                ->right_mismatches == 0,
            "Pair-level right mismatch count"
        );

        expect(
            mask_product_it
                ->left_mismatch_mask ==
                std::uint64_t{1},
            "Pair-level left mask preserves original 5-prime position 0"
        );

        expect(
            mask_product_it
                ->right_mismatch_mask == 0,
            "Pair-level exact right mask is zero"
        );


        bool invalid_range_rejected =
            false;

        try {
            static_cast<void>(
                engine.search(
                    primer1,
                    primer2,
                    12,
                    0,
                    3000,
                    50
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            invalid_range_rejected =
                true;
        }

        expect(
            invalid_range_rejected,
            "Invalid amplicon range rejected"
        );

        std::cout
            << "primer1_single_hits\t"
            << result.primer1_single_hit_count
            << '\n';

        std::cout
            << "primer2_single_hits\t"
            << result.primer2_single_hit_count
            << '\n';

        std::cout
            << "valid_amplicons\t"
            << result.amplicon_count()
            << '\n';

        std::cout
            << "All primer-pair search tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
