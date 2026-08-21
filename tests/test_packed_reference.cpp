#include "primerpair/packed_reference.hpp"

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

}  // namespace

int main() {
    try {
        /*
         * Positions:
         *
         * 0  A
         * 1  C
         * 2  G
         * 3  T
         * 4  N
         *
         * R/Y/W/... should normalize to N.
         */
        const std::string sequence =
            "ACGTNRYWSKMBDHVACGT";

        const primerpair::PackedReference
            reference(sequence);

        expect(
            reference.size() ==
                sequence.size(),
            "Packed reference length"
        );

        expect(
            reference.base_at(0) == 'A',
            "Decode A"
        );

        expect(
            reference.base_at(1) == 'C',
            "Decode C"
        );

        expect(
            reference.base_at(2) == 'G',
            "Decode G"
        );

        expect(
            reference.base_at(3) == 'T',
            "Decode T"
        );

        expect(
            reference.base_at(4) == 'N',
            "Decode N"
        );

        /*
         * IUPAC ambiguity normalization.
         */
        expect(
            reference.base_at(5) == 'N',
            "R normalized to N"
        );

        expect(
            reference.base_at(6) == 'N',
            "Y normalized to N"
        );

        expect(
            reference.base_at(7) == 'N',
            "W normalized to N"
        );

        /*
         * Last four bases = ACGT.
         */
        expect(
            reference.base_at(15) == 'A',
            "Trailing A"
        );

        expect(
            reference.base_at(18) == 'T',
            "Trailing T"
        );

        expect(
            reference.bounded_hamming_distance(
                15,
                "ACGT",
                0
            ) == 0,
            "Exact Hamming verification"
        );

        expect(
            reference.bounded_hamming_distance(
                15,
                "ACGA",
                1
            ) == 1,
            "One mismatch verification"
        );

        expect(
            reference.bounded_hamming_distance(
                15,
                "TCGA",
                1
            ) > 1,
            "Bounded Hamming early rejection"
        );

        expect(
            reference.bounded_hamming_distance(
                15,
                "acgt",
                0
            ) == 0,
            "Lowercase query normalization"
        );

        /*
         * Query crosses N.
         *
         * Candidate must be rejected regardless
         * of ordinary mismatch count.
         */
        expect(
            reference.bounded_hamming_distance(
                3,
                "TA",
                3
            ) > 3,
            "Reference N rejects candidate"
        );

        bool invalid_query_rejected = false;

        try {
            static_cast<void>(
                reference
                    .bounded_hamming_distance(
                        15,
                        "ACGN",
                        1
                    )
            );

        } catch (
            const std::invalid_argument&
        ) {
            invalid_query_rejected = true;
        }

        expect(
            invalid_query_rejected,
            "Invalid query nucleotide rejection"
        );

        bool position_rejected = false;

        try {
            static_cast<void>(
                reference.base_at(
                    reference.size()
                )
            );

        } catch (
            const std::out_of_range&
        ) {
            position_rejected = true;
        }

        expect(
            position_rejected,
            "Out-of-range base rejection"
        );

        bool span_rejected = false;

        try {
            static_cast<void>(
                reference
                    .bounded_hamming_distance(
                        reference.size() - 1,
                        "ACGT",
                        1
                    )
            );

        } catch (
            const std::out_of_range&
        ) {
            span_rejected = true;
        }

        expect(
            span_rejected,
            "Out-of-range span rejection"
        );

        expect(
            reference.memory_bytes() > 0,
            "Packed reference memory accounting"
        );

        std::cout
            << "packed_reference_bytes\t"
            << reference.memory_bytes()
            << '\n';

        std::cout
            << "All packed reference tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
