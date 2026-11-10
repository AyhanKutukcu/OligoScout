#include "primerpair/packed_bwt.hpp"
#include "primerpair/rank_support.hpp"

#include <array>
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
    using primerpair::CheckpointRank;
    using primerpair::PackedBWT;

    try {
        /*
         * Birden fazla 64-bit word ve birden fazla
         * 128-bp rank checkpoint sınırı içerir.
         */
        std::string bwt;
        bwt.reserve(333);

        const std::string motif =
            "ACGTN";

        for (std::size_t i = 0;
             i < 333;
             ++i) {

            bwt.push_back(
                motif.at(
                    i % motif.size()
                )
            );
        }

        /*
         * Tam olarak bir sentinel.
         *
         * 128 sınırından sonra ve 256 sınırından
         * önce olacak şekilde yerleştiriyoruz.
         */
        bwt.at(173) = '$';

        CheckpointRank rank;
        rank.build(
            bwt,
            128
        );

        const PackedBWT packed(
            bwt
        );

        const std::array<char, 6> symbols{
            '$',
            'A',
            'C',
            'G',
            'N',
            'T'
        };

        std::uint64_t comparisons = 0;

        /*
         * position exclusive-prefix semantiğine sahip:
         *
         * rank(symbol, position)
         * -> [0, position)
         */
        for (const char symbol : symbols) {
            for (
                std::uint64_t position = 0;
                position <=
                    static_cast<std::uint64_t>(
                        bwt.size()
                    );
                ++position
            ) {
                const std::uint64_t scalar =
                    rank.rank(
                        bwt,
                        symbol,
                        position
                    );

                const std::uint64_t packed_value =
                    rank.rank(
                        packed,
                        symbol,
                        position
                    );

                if (scalar != packed_value) {
                    throw std::runtime_error(
                        "Packed/scalar rank mismatch for symbol " +
                        std::string(1, symbol) +
                        " at position " +
                        std::to_string(position) +
                        ": scalar=" +
                        std::to_string(scalar) +
                        " packed=" +
                        std::to_string(packed_value)
                    );
                }

                ++comparisons;
            }
        }

        expect(
            comparisons ==
                static_cast<std::uint64_t>(
                    symbols.size()
                ) *
                (
                    static_cast<std::uint64_t>(
                        bwt.size()
                    ) +
                    1
                ),
            "All scalar/packed rank comparisons executed"
        );

        expect(
            rank.rank(
                packed,
                '$',
                174
            ) == 1,
            "Packed rank includes sentinel after position 173"
        );

        expect(
            rank.rank(
                packed,
                '$',
                173
            ) == 0,
            "Packed rank excludes sentinel at exclusive boundary"
        );

        bool invalid_symbol_rejected = false;

        try {
            static_cast<void>(
                rank.rank(
                    packed,
                    'X',
                    100
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            invalid_symbol_rejected = true;
        }

        expect(
            invalid_symbol_rejected,
            "Packed rank invalid-symbol rejection"
        );

        bool out_of_range_rejected = false;

        try {
            static_cast<void>(
                rank.rank(
                    packed,
                    'A',
                    static_cast<std::uint64_t>(
                        bwt.size()
                    ) +
                    1
                )
            );
        } catch (
            const std::out_of_range&
        ) {
            out_of_range_rejected = true;
        }

        expect(
            out_of_range_rejected,
            "Packed rank out-of-range rejection"
        );

        std::cout
            << "scalar_packed_rank_comparisons\t"
            << comparisons
            << '\n';

        std::cout
            << "All packed-rank equivalence tests passed.\n";

        return 0;

    } catch (const std::exception& error) {
        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
