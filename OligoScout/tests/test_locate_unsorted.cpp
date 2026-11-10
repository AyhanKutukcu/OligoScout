#include "primerpair/bidirectional_fm_index.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

int main() {
    try {
        using namespace primerpair;

        constexpr std::array<char, 4>
            alphabet{
                'A',
                'C',
                'G',
                'T'
            };

        std::string reference;
        reference.reserve(
            30000
        );

        std::uint64_t rng =
            0xD1B54A32D192ED03ULL;

        for (
            std::size_t i = 0;
            i < 24000;
            ++i
        ) {
            rng =
                rng *
                6364136223846793005ULL +
                1442695040888963407ULL;

            reference.push_back(
                alphabet.at(
                    static_cast<std::size_t>(
                        (rng >> 32U) &
                        3ULL
                    )
                )
            );

            if (
                i != 0 &&
                i % 701 == 0
            ) {
                reference +=
                    "ACGTACGTACGT"
                    "AAAAAAAAAAAA"
                    "ACGTACGTACGT";
            }
        }

        BidirectionalFMIndex index(
            reference,
            8
        );

        std::size_t checks = 0;

        for (
            std::size_t trial = 0;
            trial < 3000;
            ++trial
        ) {
            rng =
                rng *
                6364136223846793005ULL +
                1442695040888963407ULL;

            const std::size_t length =
                4 +
                static_cast<std::size_t>(
                    rng % 21
                );

            rng =
                rng *
                6364136223846793005ULL +
                1442695040888963407ULL;

            const std::size_t start =
                static_cast<std::size_t>(
                    rng %
                    (
                        reference.size() -
                        length +
                        1
                    )
                );

            const std::string pattern =
                reference.substr(
                    start,
                    length
                );

            const auto state =
                index.search(
                    pattern
                );

            const auto sorted =
                index.locate(
                    state
                );

            auto unsorted =
                index.locate_unsorted(
                    state
                );

            std::sort(
                unsorted.begin(),
                unsorted.end()
            );

            if (
                sorted !=
                unsorted
            ) {
                std::cerr
                    << "locate_unsorted mismatch "
                    << "trial="
                    << trial
                    << '\n';

                return 1;
            }

            ++checks;
        }

        std::cout
            << "trials\t3000\n";

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
