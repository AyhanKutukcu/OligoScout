#include "primerpair/sensitive_pair_constrained_search.hpp"

#include <cstddef>
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

        const std::string rc1 =
            primerpair::
                reverse_complement(
                    primer1
                );

        const std::string rc2 =
            primerpair::
                reverse_complement(
                    primer2
                );

        /*
         * Terminal 3-prime mismatch for Primer1.
         */
        std::string primer1_terminal =
            primer1;

        primer1_terminal.back() =
            mutate_base(
                primer1_terminal.back()
            );

        /*
         * Internal mismatches for Primer2 reverse
         * binding site.
         */
        std::string rc2_two_mm =
            rc2;

        rc2_two_mm.at(2) =
            mutate_base(
                rc2_two_mm.at(2)
            );

        rc2_two_mm.at(11) =
            mutate_base(
                rc2_two_mm.at(11)
            );

        std::string reference(
            30,
            'N'
        );

        /*
         * Product 1: exact 100 bp.
         */
        reference += primer1;

        reference.append(
            60,
            'N'
        );

        reference += rc2;

        reference.append(
            200,
            'N'
        );

        /*
         * Product 2:
         * Primer1 terminal 3' mismatch.
         */
        reference += primer1_terminal;

        reference.append(
            60,
            'N'
        );

        reference += rc2;

        reference.append(
            200,
            'N'
        );

        /*
         * Product 3:
         * two mismatches on Primer2 binding site.
         */
        reference += primer1;

        reference.append(
            60,
            'N'
        );

        reference += rc2_two_mm;

        reference.append(
            200,
            'N'
        );

        /*
         * Opposite biological orientation:
         *
         * Primer2 --->       <--- Primer1
         */
        reference += primer2;

        reference.append(
            50,
            'N'
        );

        reference += rc1;

        reference.append(
            30,
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
            global_engine(
                index,
                packed
            );

        const primerpair::
            SensitivePairConstrainedSearchEngine
                constrained_engine(
                    index,
                    packed
                );

        const auto global =
            global_engine.search(
                primer1,
                primer2,
                primerpair::
                    SearchProfile::
                        Sensitive,
                12,
                3,
                50,
                120
            );

        const auto constrained =
            constrained_engine.search(
                primer1,
                primer2,
                3,
                50,
                120
            );

        expect(
            constrained
                .pair_result
                .amplicons ==
            global.amplicons,
            "Pair-constrained k3 equals global-global SENSITIVE"
        );

        expect(
            constrained
                .pair_result
                .amplicon_count() ==
            global.amplicon_count(),
            "Pair-constrained amplicon count"
        );

        expect(
            !constrained
                .pair_result
                .empty(),
            "Pair-constrained finds synthetic products"
        );

        expect(
            constrained
                .anchor_global_hit_count >
            0,
            "Global anchor hits found"
        );

        expect(
            constrained
                .mate_local_hit_count >
            0,
            "Local mate hits found"
        );

        expect(
            constrained
                .scanned_mate_start_positions >
            0,
            "Local mate positions scanned"
        );

        bool budget_rejected =
            false;

        try {
            static_cast<void>(
                constrained_engine.search(
                    primer1,
                    primer2,
                    2,
                    50,
                    120
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            budget_rejected =
                true;
        }

        expect(
            budget_rejected,
            "v1 rejects non-k3 budget"
        );

        bool length17_rejected =
            false;

        try {
            static_cast<void>(
                constrained_engine.search(
                    std::string(
                        17,
                        'A'
                    ),
                    primer2,
                    3,
                    50,
                    120
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            length17_rejected =
                true;
        }

        expect(
            length17_rejected,
            "MVP rejects 17-nt primer"
        );


        bool length36_rejected =
            false;

        try {
            static_cast<void>(
                constrained_engine.search(
                    std::string(
                        36,
                        'A'
                    ),
                    primer2,
                    3,
                    50,
                    120
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            length36_rejected =
                true;
        }

        expect(
            length36_rejected,
            "MVP rejects 36-nt primer"
        );

        std::cout
            << "anchor_primer\t"
            << primerpair::to_string(
                   constrained
                       .anchor_primer
               )
            << '\n';

        std::cout
            << "global_anchor_hits\t"
            << constrained
                .anchor_global_hit_count
            << '\n';

        std::cout
            << "local_mate_hits\t"
            << constrained
                .mate_local_hit_count
            << '\n';

        std::cout
            << "scanned_mate_starts\t"
            << constrained
                .scanned_mate_start_positions
            << '\n';

        std::cout
            << "amplicons\t"
            << constrained
                .pair_result
                .amplicon_count()
            << '\n';

        const auto forced_p1 =
            constrained_engine.search(
                primer1,
                primer2,
                3,
                50,
                120,
                primerpair::
                    SensitivePairAnchorPolicy::
                        ForcePrimer1
            );

        const auto forced_p2 =
            constrained_engine.search(
                primer1,
                primer2,
                3,
                50,
                120,
                primerpair::
                    SensitivePairAnchorPolicy::
                        ForcePrimer2
            );

        expect(
            forced_p1.anchor_primer ==
                primerpair::PrimerIdentity::Primer1,
            "Forced Primer1 anchor selected"
        );

        expect(
            forced_p2.anchor_primer ==
                primerpair::PrimerIdentity::Primer2,
            "Forced Primer2 anchor selected"
        );

        expect(
            forced_p1.pair_result.amplicons ==
                global.amplicons &&
            forced_p2.pair_result.amplicons ==
                global.amplicons,
            "Both forced anchors remain lossless"
        );

        std::cout
            << "All pair-constrained SENSITIVE tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
