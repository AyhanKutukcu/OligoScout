#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/neural_window_predictor.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

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
                100000
            );


        IPBWTIndex index(
            reference,
            21
        );


        NeuralWindowPredictor predictor;


        std::size_t checks = 0;


        for (
            std::size_t start = 0;
            start + 21 <= reference.size();
            start += 97
        ) {
            const std::string query =
                reference.substr(
                    start,
                    21
                );


            const Interval expected =
                index.exact_search(
                    query
                );


            const std::uint64_t prediction =
                predictor.predict_row(
                    query,
                    index.row_count()
                );


            /*
             * This model was trained on the 1 Mb
             * synthetic distribution, not this 100 kb
             * index. Use a deliberately tiny window:
             * failed predictions MUST therefore exercise
             * the certified global fallback.
             */
            const auto observed =
                index.exact_prefix_search_certified_window(
                    query,
                    prediction,
                    32
                );


            if (
                expected.begin !=
                    observed.interval.begin
                ||
                expected.end !=
                    observed.interval.end
            ) {
                std::cerr
                    << "INTERVAL_MISMATCH\t"
                    << query
                    << '\t'
                    << expected.begin
                    << '\t'
                    << expected.end
                    << '\t'
                    << observed.interval.begin
                    << '\t'
                    << observed.interval.end
                    << '\n';

                throw std::runtime_error(
                    "Certified neural-window search "
                    "differs from exact IP-BWT."
                );
            }


            ++checks;
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
