#include "primerpair/primer_binding_site.hpp"

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void expect(
    const bool condition,
    const char* message
) {
    if (!condition) {
        throw std::runtime_error(
            message
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}


std::string reverse_complement(
    const std::string& sequence
) {
    std::string result;

    result.reserve(
        sequence.size()
    );

    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {
        switch (*it) {
            case 'A':
                result.push_back('T');
                break;

            case 'C':
                result.push_back('G');
                break;

            case 'G':
                result.push_back('C');
                break;

            case 'T':
                result.push_back('A');
                break;

            default:
                throw std::runtime_error(
                    "Unexpected nucleotide."
                );
        }
    }

    return result;
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string primer1 =
            "ACGTACGTTGCAACGTACGT";

        const std::string primer2 =
            "TTGGAACCTTGGAACCTTGG";


        /*
         * Introduce one mismatch in Primer2-oriented
         * target at zero-based position 1.
         */
        std::string primer2_target =
            primer2;

        primer2_target.at(1) =
            (
                primer2_target.at(1) == 'A'
                ? 'C'
                : 'A'
            );


        const std::string right_genomic =
            reverse_complement(
                primer2_target
            );


        const std::string sequence =
            std::string(
                50,
                'N'
            )
            +
            primer1
            +
            std::string(
                60,
                'N'
            )
            +
            right_genomic
            +
            std::string(
                50,
                'N'
            );


        PackedReference reference(
            sequence
        );


        PrimerPairHit hit;

        hit.left_primer =
            PrimerIdentity::Primer1;

        hit.right_primer =
            PrimerIdentity::Primer2;

        hit.left_position =
            50;

        hit.right_position =
            130;

        hit.amplicon_start =
            50;

        hit.amplicon_end_exclusive =
            150;

        hit.amplicon_length =
            100;

        hit.left_mismatches =
            0;

        hit.right_mismatches =
            1;

        hit.left_mismatch_mask =
            0;

        hit.right_mismatch_mask =
            (
                std::uint64_t{1}
                <<
                1
            );


        const auto sites =
            extract_pair_binding_sites(
                reference,
                hit,
                primer1,
                primer2
            );


        expect(
            sites.left.primer_sequence ==
                primer1,
            "Left primer sequence restored"
        );

        expect(
            sites.left.binding_sequence ==
                primer1,
            "Left genomic site oriented correctly"
        );

        expect(
            !sites.left.reverse_strand,
            "Left site marked forward"
        );

        expect(
            sites.left.mismatch_count() == 0,
            "Left site exact"
        );


        expect(
            sites.right.primer_sequence ==
                primer2,
            "Right primer sequence restored"
        );

        expect(
            sites.right.genomic_sequence ==
                right_genomic,
            "Right raw genomic sequence restored"
        );

        expect(
            sites.right.binding_sequence ==
                primer2_target,
            "Right binding site reverse-complemented"
        );

        expect(
            sites.right.reverse_strand,
            "Right site marked reverse"
        );

        expect(
            sites.right.mismatch_count() == 1,
            "Right site contains one mismatch"
        );

        expect(
            sites.right.mismatches.at(0)
                .primer_position
            ==
            1,
            "Mismatch position is zero-based 1"
        );

        expect(
            sites.right.mismatches.at(0)
                .distance_from_3prime
            ==
            primer2.size() - 2,
            "3-prime mismatch distance correct"
        );

        expect(
            sites.right.mismatch_mask ==
            (
                std::uint64_t{1}
                <<
                1
            ),
            "Observed mismatch mask matches search hit"
        );


        std::cout
            << "left_binding\t"
            << sites.left.binding_sequence
            << '\n';

        std::cout
            << "right_primer\t"
            << sites.right.primer_sequence
            << '\n';

        std::cout
            << "right_binding\t"
            << sites.right.binding_sequence
            << '\n';

        std::cout
            << "mismatch_5prime_position_1based\t"
            << (
                sites.right
                    .mismatches
                    .at(0)
                    .primer_position
                +
                1
            )
            << '\n';

        std::cout
            << "distance_from_3prime\t"
            << sites.right
                .mismatches
                .at(0)
                .distance_from_3prime
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
