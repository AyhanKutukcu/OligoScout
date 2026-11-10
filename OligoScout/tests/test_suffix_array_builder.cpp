#include "primerpair/suffix_array_builder.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

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


std::vector<std::uint32_t>
naive_suffix_array(
    const std::string& text
) {
    std::vector<std::uint32_t>
        sa(
            text.size()
        );

    std::iota(
        sa.begin(),
        sa.end(),
        std::uint32_t{0}
    );


    std::sort(
        sa.begin(),
        sa.end(),
        [&text](
            const std::uint32_t lhs,
            const std::uint32_t rhs
        ) {
            return
                std::lexicographical_compare(
                    text.begin() +
                        static_cast<
                            std::ptrdiff_t
                        >(lhs),
                    text.end(),
                    text.begin() +
                        static_cast<
                            std::ptrdiff_t
                        >(rhs),
                    text.end()
                );
        }
    );

    return sa;
}


void compare(
    const std::string& text
) {
    const auto expected =
        naive_suffix_array(
            text
        );

    const auto observed =
        primerpair::
            build_suffix_array_prefix_doubling(
                text
            );

    if (
        expected !=
        observed
    ) {
        std::cerr
            << "SA mismatch for text: "
            << text
            << '\n';

        throw std::runtime_error(
            "Prefix-doubling SA differs "
            "from naive suffix ordering."
        );
    }
}

}  // namespace


int main() {
    try {
        compare(
            "CATTATTAGGA$"
        );

        compare(
            "ACGTACGTACGT$"
        );

        compare(
            "AAAAAAAAAAAA$"
        );

        compare(
            "TGCATGCATGCA$"
        );

        compare(
            "ACGTNNNNACGT$"
        );

        compare(
            "GATTACAGATTACA$"
        );


        /*
         * Larger deterministic reference.
         *
         * FM and IP-BWT now share this SA builder,
         * so validate it independently against
         * naive lexicographic suffix ordering.
         */
        {
            constexpr char alphabet[] = {
                'A',
                'C',
                'G',
                'T'
            };

            std::uint64_t state =
                0x9E3779B97F4A7C15ULL;

            std::string long_text;

            long_text.reserve(
                5000
            );

            for (
                std::size_t i = 0;
                i < 4096;
                ++i
            ) {
                state =
                    state *
                    6364136223846793005ULL
                    +
                    1442695040888963407ULL;

                const std::size_t index =
                    static_cast<std::size_t>(
                        (
                            state >>
                            32
                        )
                        &
                        3ULL
                    );

                long_text.push_back(
                    alphabet[index]
                );
            }


            const std::string motif =
                "ACGTACGTACGTACGTACGTACGT"
                "TTTTAAAACCCCGGGG"
                "ACGTACGTACGTACGTACGTACGT";

            for (
                std::size_t repeat = 0;
                repeat < 12;
                ++repeat
            ) {
                long_text +=
                    motif;
            }


            long_text.push_back(
                '$'
            );


            compare(
                long_text
            );


            expect(
                true,
                "Large deterministic suffix array "
                "matches naive ground truth"
            );
        }


        expect(
            true,
            "Shared suffix-array builder matches "
            "naive ordering"
        );


        bool missing_sentinel_rejected =
            false;

        try {
            (void)primerpair::
                build_suffix_array_prefix_doubling(
                    "ACGT"
                );

        } catch (
            const std::invalid_argument&
        ) {
            missing_sentinel_rejected =
                true;
        }

        expect(
            missing_sentinel_rejected,
            "Missing sentinel rejected"
        );


        bool duplicate_sentinel_rejected =
            false;

        try {
            (void)primerpair::
                build_suffix_array_prefix_doubling(
                    "AC$GT$"
                );

        } catch (
            const std::invalid_argument&
        ) {
            duplicate_sentinel_rejected =
                true;
        }

        expect(
            duplicate_sentinel_rejected,
            "Duplicate sentinel rejected"
        );


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
