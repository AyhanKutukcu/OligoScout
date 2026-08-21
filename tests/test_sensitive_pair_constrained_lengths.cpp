#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/search_profile.hpp"
#include "primerpair/sensitive_adaptive_search.hpp"
#include "primerpair/sensitive_pair_constrained_search.hpp"

namespace {

void expect(
    const bool condition,
    const std::string& message
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " + message
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}


char complement(
    const char base
) {
    switch (base) {
        case 'A':
            return 'T';

        case 'C':
            return 'G';

        case 'G':
            return 'C';

        case 'T':
            return 'A';

        default:
            throw std::invalid_argument(
                "Invalid nucleotide."
            );
    }
}


std::string reverse_complement(
    const std::string_view sequence
) {
    std::string output;

    output.reserve(
        sequence.size()
    );

    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {
        output.push_back(
            complement(
                *it
            )
        );
    }

    return output;
}


/*
 * Deterministic pseudo-random A/C/G/T reference.
 *
 * Fixed seed gives reproducible tests while avoiding
 * a trivial low-complexity ACGTACGT... reference.
 */
std::string make_reference(
    const std::size_t length
) {
    static constexpr
        std::array<char, 4>
            bases{
                'A',
                'C',
                'G',
                'T'
            };

    std::uint32_t state =
        0x9E3779B9u;

    std::string reference;

    reference.reserve(
        length
    );

    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        /*
         * xorshift32
         */
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;

        const std::size_t base_index =
            static_cast<std::size_t>(
                (
                    state >>
                    30
                )
                &
                0x3u
            );

        reference.push_back(
            bases.at(
                base_index
            )
        );
    }

    return reference;
}

}  // namespace


int main() {
    try {
        static constexpr
            std::array<std::size_t, 18>
                primer_lengths{
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                    31,
                    32,
                    33,
                    34,
                    35
                };

        constexpr std::uint64_t
            forward_position = 80;

        constexpr std::uint64_t
            reverse_position = 220;

        constexpr std::uint64_t
            min_amplicon = 50;

        constexpr std::uint64_t
            max_amplicon = 300;

        constexpr std::size_t
            max_mismatches = 3;

        constexpr std::size_t
            legacy_anchor_length = 12;


        const std::string reference =
            make_reference(
                600
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

        const primerpair::
            SensitiveAdaptiveSearchEngine
                adaptive_engine(
                    index,
                    packed
                );


        for (
            const std::size_t primer_length :
            primer_lengths
        ) {
            const std::string primer1 =
                reference.substr(
                    static_cast<std::size_t>(
                        forward_position
                    ),
                    primer_length
                );

            const std::string primer2 =
                reverse_complement(
                    reference.substr(
                        static_cast<std::size_t>(
                            reverse_position
                        ),
                        primer_length
                    )
                );


            /*
             * Independent global-global SENSITIVE
             * reference result.
             */
            const auto global =
                global_engine.search(
                    primer1,
                    primer2,
                    primerpair::
                        SearchProfile::
                            Sensitive,
                    legacy_anchor_length,
                    max_mismatches,
                    min_amplicon,
                    max_amplicon
                );


            /*
             * Default pair-constrained execution.
             */
            const auto constrained =
                constrained_engine.search(
                    primer1,
                    primer2,
                    max_mismatches,
                    min_amplicon,
                    max_amplicon
                );


            /*
             * Both possible anchor identities must
             * remain lossless.
             */
            const auto forced_p1 =
                constrained_engine.search(
                    primer1,
                    primer2,
                    max_mismatches,
                    min_amplicon,
                    max_amplicon,
                    primerpair::
                        SensitivePairAnchorPolicy::
                            ForcePrimer1
                );

            const auto forced_p2 =
                constrained_engine.search(
                    primer1,
                    primer2,
                    max_mismatches,
                    min_amplicon,
                    max_amplicon,
                    primerpair::
                        SensitivePairAnchorPolicy::
                            ForcePrimer2
                );


            expect(
                !global.amplicons.empty(),
                std::to_string(
                    primer_length
                ) +
                "-nt designed amplicon recovered"
            );


            expect(
                constrained
                    .pair_result
                    .amplicons
                ==
                global
                    .amplicons,
                std::to_string(
                    primer_length
                ) +
                "-nt constrained equals "
                "global SENSITIVE"
            );


            expect(
                forced_p1
                    .pair_result
                    .amplicons
                ==
                global
                    .amplicons,
                std::to_string(
                    primer_length
                ) +
                "-nt forced P1 remains lossless"
            );


            expect(
                forced_p2
                    .pair_result
                    .amplicons
                ==
                global
                    .amplicons,
                std::to_string(
                    primer_length
                ) +
                "-nt forced P2 remains lossless"
            );


            expect(
                forced_p1.anchor_primer ==
                    primerpair::
                        PrimerIdentity::
                            Primer1,
                std::to_string(
                    primer_length
                ) +
                "-nt forced P1 selected"
            );


            expect(
                forced_p2.anchor_primer ==
                    primerpair::
                        PrimerIdentity::
                            Primer2,
                std::to_string(
                    primer_length
                ) +
                "-nt forced P2 selected"
            );


            /*
             * Validate the production hybrid k3
             * routing policy.
             *
             * 18..29:
             *   estimator + occurrence router.
             *
             * 30..35:
             *   direct Candidate, estimator skipped.
             */
            const auto adaptive =
                adaptive_engine.search(
                    primer1,
                    max_mismatches
                );

            if (primer_length <= 29) {

                expect(
                    adaptive.estimator_used,
                    std::to_string(
                        primer_length
                    ) +
                    "-nt k3 uses short-primer estimator"
                );

            } else {

                expect(
                    !adaptive.estimator_used,
                    std::to_string(
                        primer_length
                    ) +
                    "-nt k3 skips estimator"
                );

                expect(
                    adaptive.backend ==
                        primerpair::
                            SensitiveAdaptiveBackend::
                                Candidate,
                    std::to_string(
                        primer_length
                    ) +
                    "-nt k3 uses direct Candidate"
                );
            }


            std::cout
                << "length\t"
                << primer_length
                << "\tamplicons\t"
                << global.amplicons.size()
                << "\tverified_candidates\t"
                << constrained
                    .scanned_mate_start_positions
                << '\n';
        }


        /*
         * MVP boundaries themselves.
         */
        bool rejected17 = false;

        try {
            static_cast<void>(
                constrained_engine.search(
                    std::string(
                        17,
                        'A'
                    ),
                    std::string(
                        18,
                        'C'
                    ),
                    3,
                    50,
                    300
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            rejected17 = true;
        }

        expect(
            rejected17,
            "17-nt primer rejected"
        );


        bool rejected36 = false;

        try {
            static_cast<void>(
                constrained_engine.search(
                    std::string(
                        36,
                        'A'
                    ),
                    std::string(
                        18,
                        'C'
                    ),
                    3,
                    50,
                    300
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            rejected36 = true;
        }

        expect(
            rejected36,
            "36-nt primer rejected"
        );


        std::cout
            << "All 18–35 nt pair-constrained "
            << "generalization tests passed.\n";

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
