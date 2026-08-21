#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_assembler.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

std::string make_reference(
    const std::size_t length
) {
    constexpr char bases[] = {
        'A', 'C', 'G', 'T'
    };

    std::uint64_t state =
        0xC6BC279692B5CC83ULL;

    std::string reference;

    reference.reserve(
        length
    );

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

}  // namespace


int main() {
    try {
        using namespace primerpair;

        constexpr std::size_t
            primer_length = 24;

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::uint64_t
            min_amplicon = 50;

        constexpr std::uint64_t
            max_amplicon = 1000;


        const std::string reference =
            make_reference(
                100000
            );


        PackedReference packed(
            reference
        );

        BidirectionalFMIndex index(
            reference
        );

        StrandAwarePrimerSearchEngine
            strand_engine(
                index,
                packed
            );

        PrimerPairSearchEngine
            legacy_pair_engine(
                index,
                packed
            );


        std::size_t checks = 0;


        for (
            std::size_t trial = 0;
            trial < 500;
            ++trial
        ) {
            const std::size_t
                position1 =
                    (
                        trial *
                        3571 +
                        101
                    ) %
                    (
                        reference.size() -
                        1500
                    );

            const std::size_t
                position2 =
                    position1 +
                    100 +
                    (
                        trial %
                        700
                    );


            const std::string primer1 =
                reference.substr(
                    position1,
                    primer_length
                );


            /*
             * Primer2 is chosen from genomic reference
             * then reverse-complemented so the intended
             * pair has inward-facing PCR geometry.
             */
            const std::string genomic2 =
                reference.substr(
                    position2,
                    primer_length
                );

            const std::string primer2 =
                reverse_complement(
                    genomic2
                );


            const auto primer1_result =
                strand_engine.search(
                    primer1,
                    anchor_length,
                    3
                );

            const auto primer2_result =
                strand_engine.search(
                    primer2,
                    anchor_length,
                    3
                );


            const auto observed =
                assemble_primer_pair_hits(
                    primer1,
                    primer1_result.hits,
                    primer2,
                    primer2_result.hits,
                    min_amplicon,
                    max_amplicon
                );


            const auto expected =
                legacy_pair_engine.search(
                    primer1,
                    primer2,
                    anchor_length,
                    3,
                    min_amplicon,
                    max_amplicon
                );


            if (
                observed.amplicons !=
                expected.amplicons
            ) {
                std::cerr
                    << "PAIR_ASSEMBLER_MISMATCH"
                    << '\t'
                    << "trial="
                    << trial
                    << '\t'
                    << "expected="
                    << expected.amplicons.size()
                    << '\t'
                    << "observed="
                    << observed.amplicons.size()
                    << '\n';

                return 1;
            }


            if (
                observed.primer1_single_hit_count !=
                    expected.primer1_single_hit_count ||
                observed.primer2_single_hit_count !=
                    expected.primer2_single_hit_count
            ) {
                throw std::runtime_error(
                    "Primer hit-count mismatch."
                );
            }


            ++checks;
        }


        std::cout
            << "checks\t"
            << checks
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
