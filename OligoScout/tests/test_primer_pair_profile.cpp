#include "primerpair/primer_pair_search.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

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

char mutate_base(
    const char base
) {
    return
        base == 'A'
            ? 'C'
            : 'A';
}

}  // namespace


int main() {
    try {
        const std::string primer1 =
            "ACGTTGCAAGTCCTGAACGA";

        const std::string primer2 =
            "TGCAGATCCGTACGATGTCA";

        const std::string rc2 =
            primerpair::reverse_complement(
                primer2
            );

        /*
         * Primer1 ORIGINAL position 19:
         * biological terminal 3-prime mismatch.
         */
        std::string terminal_mismatch =
            primer1;

        terminal_mismatch.back() =
            mutate_base(
                terminal_mismatch.back()
            );

        std::string reference(
            10,
            'N'
        );

        /*
         * Exact 100-bp product.
         */
        const std::uint64_t
            exact_left =
                reference.size();

        reference += primer1;

        reference.append(
            60,
            'N'
        );

        const std::uint64_t
            exact_right =
                reference.size();

        reference += rc2;

        /*
         * Prevent cross-pairing.
         */
        reference.append(
            200,
            'N'
        );

        /*
         * Second 100-bp product.
         *
         * Primer1 binding site contains a
         * terminal 3-prime mismatch.
         */
        const std::uint64_t
            mismatch_left =
                reference.size();

        reference +=
            terminal_mismatch;

        reference.append(
            60,
            'N'
        );

        const std::uint64_t
            mismatch_right =
                reference.size();

        reference += rc2;

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
         * --------------------------------------------------
         * Legacy STRICT == explicit STRICT
         * --------------------------------------------------
         */

        const auto legacy =
            engine.search(
                primer1,
                primer2,
                12,
                1,
                50,
                120
            );

        const auto strict =
            engine.search(
                primer1,
                primer2,
                primerpair::
                    SearchProfile::
                        Strict,
                12,
                1,
                50,
                120
            );

        expect(
            strict.amplicons ==
                legacy.amplicons,
            "Explicit STRICT pair hits equal legacy API"
        );

        expect(
            strict.primer1_single_hit_count ==
                legacy.primer1_single_hit_count,
            "STRICT Primer1 single-hit count unchanged"
        );

        expect(
            strict.primer2_single_hit_count ==
                legacy.primer2_single_hit_count,
            "STRICT Primer2 single-hit count unchanged"
        );

        expect(
            strict.amplicon_count() == 1,
            "STRICT retains only exact-anchor product"
        );

        expect(
            strict.amplicons.front()
                .left_position ==
                exact_left,
            "STRICT exact product left coordinate"
        );

        expect(
            strict.amplicons.front()
                .right_position ==
                exact_right,
            "STRICT exact product right coordinate"
        );


        /*
         * --------------------------------------------------
         * SENSITIVE
         * --------------------------------------------------
         */

        const auto sensitive =
            engine.search(
                primer1,
                primer2,
                primerpair::
                    SearchProfile::
                        Sensitive,
                12,
                1,
                50,
                120
            );

        expect(
            sensitive.amplicon_count() == 2,
            "SENSITIVE recovers both PCR products"
        );

        /*
         * Exact product must remain.
         */
        const auto exact_it =
            std::find_if(
                sensitive.amplicons.begin(),
                sensitive.amplicons.end(),
                [
                    exact_left,
                    exact_right
                ](
                    const auto& product
                ) {
                    return
                        product.left_primer ==
                            primerpair::
                                PrimerIdentity::
                                    Primer1
                        &&
                        product.right_primer ==
                            primerpair::
                                PrimerIdentity::
                                    Primer2
                        &&
                        product.left_position ==
                            exact_left
                        &&
                        product.right_position ==
                            exact_right
                        &&
                        product.left_mismatches == 0
                        &&
                        product.right_mismatches == 0;
                }
            );

        expect(
            exact_it !=
                sensitive.amplicons.end(),
            "SENSITIVE preserves exact pair"
        );

        /*
         * ORIGINAL primer position 19.
         */
        const std::uint64_t
            terminal_mask =
                std::uint64_t{1} << 19;

        const auto mismatch_it =
            std::find_if(
                sensitive.amplicons.begin(),
                sensitive.amplicons.end(),
                [
                    mismatch_left,
                    mismatch_right,
                    terminal_mask
                ](
                    const auto& product
                ) {
                    return
                        product.left_primer ==
                            primerpair::
                                PrimerIdentity::
                                    Primer1
                        &&
                        product.right_primer ==
                            primerpair::
                                PrimerIdentity::
                                    Primer2
                        &&
                        product.left_position ==
                            mismatch_left
                        &&
                        product.right_position ==
                            mismatch_right
                        &&
                        product.left_mismatches == 1
                        &&
                        product.right_mismatches == 0
                        &&
                        product.left_mismatch_mask ==
                            terminal_mask
                        &&
                        product.right_mismatch_mask == 0
                        &&
                        product.amplicon_length == 100;
                }
            );

        expect(
            mismatch_it !=
                sensitive.amplicons.end(),
            "SENSITIVE pair recovers terminal 3-prime mismatch"
        );

        /*
         * STRICT must not contain the same product.
         */
        const auto strict_mismatch_it =
            std::find_if(
                strict.amplicons.begin(),
                strict.amplicons.end(),
                [
                    mismatch_left,
                    mismatch_right
                ](
                    const auto& product
                ) {
                    return
                        product.left_position ==
                            mismatch_left
                        &&
                        product.right_position ==
                            mismatch_right;
                }
            );

        expect(
            strict_mismatch_it ==
                strict.amplicons.end(),
            "STRICT pair rejects terminal 3-prime mismatch"
        );

        /*
         * Genomic ordering must remain intact.
         */
        expect(
            sensitive.amplicons.at(0)
                .amplicon_start <
            sensitive.amplicons.at(1)
                .amplicon_start,
            "SENSITIVE amplicons remain genomic-order sorted"
        );

        /*
         * Invalid range protection must also work
         * through the SENSITIVE path.
         */
        bool invalid_range_rejected =
            false;

        try {
            static_cast<void>(
                engine.search(
                    primer1,
                    primer2,
                    primerpair::
                        SearchProfile::
                            Sensitive,
                    12,
                    1,
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
            "SENSITIVE invalid amplicon range rejected"
        );

        std::cout
            << "strict_amplicons\t"
            << strict.amplicon_count()
            << '\n';

        std::cout
            << "sensitive_amplicons\t"
            << sensitive.amplicon_count()
            << '\n';

        std::cout
            << "All profile-aware primer-pair tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
