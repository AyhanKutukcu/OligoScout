#include "primerpair/ip_bwt_index.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

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

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string reference =
            "ACGTACGTACGTACGTACGT"
            "GATTACAGATTACAGATTACA"
            "TTTTAAAACCCCGGGG"
            "ACACACACGTGTGTGT"
            "TGCATGCATGCATGCA";


        IPBWTIndex index(
            reference,
            3
        );


        constexpr char alphabet[] = {
            'A',
            'C',
            'G',
            'T'
        };


        std::size_t checks = 0;


        /*
         * Exhaustively generate groups of four
         * 3-mers.
         */
        for (
            std::size_t a = 0;
            a < 4;
            ++a
        ) {
            for (
                std::size_t b = 0;
                b < 4;
                ++b
            ) {
                for (
                    std::size_t c = 0;
                    c < 4;
                    ++c
                ) {
                    std::string first;

                    first.push_back(
                        alphabet[a]
                    );

                    first.push_back(
                        alphabet[b]
                    );

                    first.push_back(
                        alphabet[c]
                    );


                    std::array<std::string, 4>
                        storage{
                            first,
                            "AAA",
                            "CGT",
                            "TTT"
                        };


                    const std::array<
                        std::string_view,
                        4
                    >
                    chunks{
                        storage[0],
                        storage[1],
                        storage[2],
                        storage[3]
                    };


                    const std::array<
                        std::uint64_t,
                        4
                    >
                    rows{
                        0,
                        index.row_count() / 4,
                        index.row_count() / 2,
                        index.row_count()
                    };


                    const auto batch =
                        index.lower_bound_batch4(
                            chunks,
                            rows
                        );


                    for (
                        std::size_t lane = 0;
                        lane < 4;
                        ++lane
                    ) {
                        const auto scalar =
                            index.lower_bound(
                                chunks.at(lane),
                                rows.at(lane)
                            );


                        if (
                            batch.at(lane) !=
                            scalar
                        ) {
                            std::cerr
                                << "BATCH_MISMATCH\t"
                                << lane
                                << '\t'
                                << chunks.at(lane)
                                << '\t'
                                << scalar
                                << '\t'
                                << batch.at(lane)
                                << '\n';

                            throw std::runtime_error(
                                "Batch4 lower_bound differs "
                                "from scalar lower_bound."
                            );
                        }


                        ++checks;
                    }
                }
            }
        }


        expect(
            checks == 256,
            "Batch4/scalar lower-bound matrix completed"
        );


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
