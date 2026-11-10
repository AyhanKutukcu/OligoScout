#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/neural_window_predictor.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

std::string make_reference(
    const std::size_t length
) {
    constexpr char alphabet[] = {
        'A',
        'C',
        'G',
        'T'
    };


    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;


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
            6364136223846793005ULL
            +
            1442695040888963407ULL;


        reference.push_back(
            alphabet[
                static_cast<std::size_t>(
                    (
                        state >>
                        32
                    )
                    &
                    3ULL
                )
            ]
        );
    }


    return reference;
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string reference =
            make_reference(
                200000
            );


        IPBWTIndex index(
            reference,
            21
        );


        NeuralWindowPredictor predictor;


        std::size_t checks = 0;


        for (
            std::size_t base = 0;
            base + 8 * 101 + 21 <
                reference.size();
            base += 997
        ) {
            std::array<std::string, 8>
                owned{};

            std::array<std::string_view, 8>
                queries{};


            for (
                std::size_t lane = 0;
                lane < 8;
                ++lane
            ) {
                owned[lane] =
                    reference.substr(
                        base +
                        lane *
                        101,
                        21
                    );

                queries[lane] =
                    owned[lane];
            }


            const auto batch_binary =
                index.exact_prefix_search_batch8(
                    queries
                );


            const auto predictions =
                predictor.predict_rows_batch8(
                    queries,
                    index.row_count()
                );


            const auto batch_neural =
                index
                .exact_prefix_search_certified_window_batch8(
                    queries,
                    predictions,
                    64
                );


            for (
                std::size_t lane = 0;
                lane < 8;
                ++lane
            ) {
                const Interval expected =
                    index.exact_search(
                        queries[lane]
                    );


                if (
                    batch_binary[lane].begin !=
                        expected.begin
                    ||
                    batch_binary[lane].end !=
                        expected.end
                ) {
                    throw std::runtime_error(
                        "Batch8 binary interval mismatch."
                    );
                }


                if (
                    batch_neural[lane].interval.begin !=
                        expected.begin
                    ||
                    batch_neural[lane].interval.end !=
                        expected.end
                ) {
                    throw std::runtime_error(
                        "Batch8 certified neural "
                        "interval mismatch."
                    );
                }


                ++checks;
            }
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
