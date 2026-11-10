#include "primerpair/bidirectional_fm_index.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

bool same_interval(
    const primerpair::Interval& a,
    const primerpair::Interval& b
) {
    return
        a.begin == b.begin &&
        a.end == b.end;
}

bool same_state(
    const primerpair::BidirectionalInterval& a,
    const primerpair::BidirectionalInterval& b
) {
    return
        same_interval(
            a.forward,
            b.forward
        )
        &&
        same_interval(
            a.reverse,
            b.reverse
        )
        &&
        a.length == b.length;
}

}  // namespace


int main() {
    try {
        using namespace primerpair;

        std::string reference;
        reference.reserve(
            24000
        );

        static constexpr
        std::array<char, 4>
            bases{
                'A',
                'C',
                'G',
                'T'
            };

        /*
         * Tekrarlı + daha karmaşık bölgeleri birlikte
         * içeren deterministik synthetic reference.
         */
        for (
            std::size_t i = 0;
            i < 16000;
            ++i
        ) {
            if (
                i != 0 &&
                i % 257 == 0
            ) {
                reference.push_back(
                    'N'
                );

            } else {
                const std::size_t code =
                    (
                        i * 17 +
                        i / 7 +
                        i / 31
                    ) %
                    bases.size();

                reference.push_back(
                    bases.at(code)
                );
            }

            if (
                i != 0 &&
                i % 503 == 0
            ) {
                reference +=
                    "ACGTACGTACGT"
                    "GATTACA"
                    "ACGTACGT";
            }
        }

        BidirectionalFMIndex index(
            reference,
            8
        );

        constexpr std::array<char, 4>
            extension_bases{
                'A',
                'C',
                'G',
                'T'
            };

        std::uint64_t rng =
            0x9E3779B97F4A7C15ULL;

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

            const std::size_t limit =
                reference.size() -
                length;

            const std::size_t start =
                static_cast<std::size_t>(
                    rng %
                    (
                        static_cast<
                            std::uint64_t
                        >(limit) +
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

            if (state.empty()) {
                std::cerr
                    << "Unexpected empty seed state "
                    << "trial="
                    << trial
                    << '\n';

                return 1;
            }

            const auto batch =
                index.extend_left_all(
                    state
                );

            for (
                std::size_t i = 0;
                i <
                    extension_bases.size();
                ++i
            ) {
                const auto scalar =
                    index.extend_left(
                        state,
                        extension_bases.at(i)
                    );

                if (
                    !same_state(
                        batch.at(i),
                        scalar
                    )
                ) {
                    std::cerr
                        << "extend_left_all mismatch "
                        << "trial="
                        << trial
                        << " base="
                        << extension_bases.at(i)
                        << '\n';

                    std::cerr
                        << "batch forward="
                        << batch.at(i)
                               .forward.begin
                        << ","
                        << batch.at(i)
                               .forward.end
                        << " reverse="
                        << batch.at(i)
                               .reverse.begin
                        << ","
                        << batch.at(i)
                               .reverse.end
                        << '\n';

                    std::cerr
                        << "scalar forward="
                        << scalar.forward.begin
                        << ","
                        << scalar.forward.end
                        << " reverse="
                        << scalar.reverse.begin
                        << ","
                        << scalar.reverse.end
                        << '\n';

                    return 1;
                }

                ++checks;
            }
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
