#include "primerpair/packed_bwt.hpp"

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

}  // namespace

int main() {
    using primerpair::PackedBWT;

    try {
        const std::string original =
            "ACTGA$TANCGTN";

        const PackedBWT packed(
            original
        );

        expect(
            packed.size() ==
                original.size(),
            "Packed BWT length"
        );

        for (std::size_t i = 0;
             i < original.size();
             ++i) {

            expect(
                packed.at(i) ==
                    original.at(i),
                "Packed BWT decode position " +
                    std::to_string(i)
            );
        }

        expect(
            packed.sentinel_position() == 5,
            "Packed BWT sentinel position"
        );

        /*
         * 13 sembol tek 64-bit word'e sığar.
         *
         * low plane  = 8 byte
         * high plane = 8 byte
         */
        expect(
            packed.sequence_memory_bytes() == 16,
            "Packed bit-plane memory"
        );

        expect(
            packed.n_mask_memory_bytes() == 8,
            "Packed N-mask memory"
        );

        expect(
            packed.count(
                'A',
                0,
                original.size()
            ) == 3,
            "Whole-range A count"
        );

        expect(
            packed.count(
                'C',
                0,
                original.size()
            ) == 2,
            "Whole-range C count"
        );

        expect(
            packed.count(
                'G',
                0,
                original.size()
            ) == 2,
            "Whole-range G count"
        );

        expect(
            packed.count(
                'T',
                0,
                original.size()
            ) == 3,
            "Whole-range T count"
        );

        expect(
            packed.count(
                'N',
                0,
                original.size()
            ) == 2,
            "Whole-range N count"
        );

        expect(
            packed.count(
                '$',
                0,
                original.size()
            ) == 1,
            "Whole-range sentinel count"
        );

        /*
         * [4,12):
         * A $ T A N C G T
         */
        expect(
            packed.count('A', 4, 12) == 2,
            "Partial-range A count"
        );

        expect(
            packed.count('C', 4, 12) == 1,
            "Partial-range C count"
        );

        expect(
            packed.count('G', 4, 12) == 1,
            "Partial-range G count"
        );

        expect(
            packed.count('T', 4, 12) == 2,
            "Partial-range T count"
        );

        expect(
            packed.count('N', 4, 12) == 1,
            "Partial-range N count"
        );

        expect(
            packed.count('$', 4, 12) == 1,
            "Partial-range sentinel count"
        );

        /*
         * 64-bit word sınırını geçen test.
         */
        std::string long_bwt(
            140,
            'A'
        );

        long_bwt.at(31) = 'C';
        long_bwt.at(63) = 'G';
        long_bwt.at(64) = 'T';
        long_bwt.at(65) = 'N';
        long_bwt.at(127) = 'C';
        long_bwt.at(128) = '$';

        const PackedBWT long_packed(
            long_bwt
        );

        expect(
            long_packed.at(63) == 'G',
            "Word-boundary decode 63"
        );

        expect(
            long_packed.at(64) == 'T',
            "Word-boundary decode 64"
        );

        expect(
            long_packed.at(65) == 'N',
            "Word-boundary decode 65"
        );

        expect(
            long_packed.count(
                'G',
                60,
                68
            ) == 1,
            "Cross-word G count"
        );

        expect(
            long_packed.count(
                'T',
                60,
                68
            ) == 1,
            "Cross-word T count"
        );

        expect(
            long_packed.count(
                'N',
                60,
                68
            ) == 1,
            "Cross-word N count"
        );

        expect(
            long_packed.count(
                '$',
                128,
                129
            ) == 1,
            "Cross-word sentinel count"
        );

        bool missing_sentinel_rejected = false;

        try {
            const PackedBWT invalid(
                "ACGTACGT"
            );
        } catch (
            const std::invalid_argument&
        ) {
            missing_sentinel_rejected = true;
        }

        expect(
            missing_sentinel_rejected,
            "Missing sentinel rejection"
        );

        bool duplicate_sentinel_rejected = false;

        try {
            const PackedBWT invalid(
                "AC$GT$AC"
            );
        } catch (
            const std::invalid_argument&
        ) {
            duplicate_sentinel_rejected = true;
        }

        expect(
            duplicate_sentinel_rejected,
            "Duplicate sentinel rejection"
        );

        bool invalid_symbol_rejected = false;

        try {
            const PackedBWT invalid(
                "ACGTX$"
            );
        } catch (
            const std::invalid_argument&
        ) {
            invalid_symbol_rejected = true;
        }

        expect(
            invalid_symbol_rejected,
            "Invalid symbol rejection"
        );

        bool invalid_range_rejected = false;

        try {
            static_cast<void>(
                packed.count(
                    'A',
                    5,
                    4
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            invalid_range_rejected = true;
        }

        expect(
            invalid_range_rejected,
            "Invalid count range rejection"
        );

        bool out_of_range_rejected = false;

        try {
            static_cast<void>(
                packed.at(
                    original.size()
                )
            );
        } catch (
            const std::out_of_range&
        ) {
            out_of_range_rejected = true;
        }

        expect(
            out_of_range_rejected,
            "Packed BWT out-of-range rejection"
        );

        std::cout
            << "All packed-BWT popcount tests passed.\n";

        return 0;

    } catch (const std::exception& error) {
        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
