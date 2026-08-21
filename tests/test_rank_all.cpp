#include "primerpair/packed_bwt.hpp"
#include "primerpair/rank_support.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

int main() {
    try {
        using namespace primerpair;

        const std::string bwt =
            "ACGTNACGTACGTNACGTACGTACGTNACGT"
            "ACGTACGTNACGTACGTACGTNACGTACGT"
            "$"
            "ACGTNACGTACGTACGTNACGTACGTACGT"
            "NACGTACGTACGTNACGTACGTACGT";

        constexpr std::array<char, 6>
            symbols{
                '$',
                'A',
                'C',
                'G',
                'N',
                'T'
            };

        PackedBWT packed(
            bwt
        );

        CheckpointRank rank;

        rank.build(
            bwt,
            32
        );

        std::size_t checks = 0;

        /*
         * count_all() vs six independent count() calls.
         */
        for (
            std::uint64_t begin = 0;
            begin <= packed.size();
            ++begin
        ) {
            for (
                std::uint64_t end = begin;
                end <= packed.size();
                ++end
            ) {
                const auto all =
                    packed.count_all(
                        begin,
                        end
                    );

                for (
                    std::size_t i = 0;
                    i < symbols.size();
                    ++i
                ) {
                    const auto scalar =
                        packed.count(
                            symbols[i],
                            begin,
                            end
                        );

                    if (
                        all[i] != scalar
                    ) {
                        std::cerr
                            << "count_all mismatch "
                            << "begin=" << begin
                            << " end=" << end
                            << " symbol="
                            << symbols[i]
                            << " all="
                            << all[i]
                            << " scalar="
                            << scalar
                            << '\n';

                        return 1;
                    }

                    ++checks;
                }
            }
        }

        /*
         * rank_all() vs six independent rank() calls.
         */
        for (
            std::uint64_t position = 0;
            position <= packed.size();
            ++position
        ) {
            const auto all =
                rank.rank_all(
                    packed,
                    position
                );

            for (
                std::size_t i = 0;
                i < symbols.size();
                ++i
            ) {
                const auto scalar =
                    rank.rank(
                        packed,
                        symbols[i],
                        position
                    );

                if (
                    all[i] != scalar
                ) {
                    std::cerr
                        << "rank_all mismatch "
                        << "position="
                        << position
                        << " symbol="
                        << symbols[i]
                        << " all="
                        << all[i]
                        << " scalar="
                        << scalar
                        << '\n';

                    return 1;
                }

                ++checks;
            }
        }

        std::cout
            << "bwt_length\t"
            << packed.size()
            << '\n';

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
