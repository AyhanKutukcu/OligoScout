#include "primerpair/neural_window_predictor.hpp"
#include "primerpair/generated/neural_window_tiny64x32_v3_golden.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace primerpair;

    try {
        NeuralWindowPredictor predictor;


        for (
            std::size_t i = 0;
            i <
            neural_window_golden::kmers.size();
            ++i
        ) {
            const float observed =
                predictor.predict_residual_rows(
                    neural_window_golden::kmers.at(
                        i
                    )
                );


            const float expected =
                neural_window_golden::residual_rows.at(
                    i
                );


            const float difference =
                std::abs(
                    observed -
                    expected
                );


            if (
                difference >
                1.0f
            ) {
                std::cerr
                    << "PREDICTION_MISMATCH\t"
                    << i
                    << '\t'
                    << observed
                    << '\t'
                    << expected
                    << '\t'
                    << difference
                    << '\n';

                throw std::runtime_error(
                    "C++ neural inference differs "
                    "from PyTorch."
                );
            }
        }


        std::cout
            << "golden_checks\t"
            << neural_window_golden::kmers.size()
            << '\n';

        std::cout
            << "parameters\t"
            << NeuralWindowPredictor::parameter_count
            << '\n';

        std::cout
            << "fp32_weight_bytes\t"
            << NeuralWindowPredictor::fp32_weight_bytes
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
