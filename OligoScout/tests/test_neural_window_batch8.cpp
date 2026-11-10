#include "primerpair/neural_window_predictor.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string_view>

int main() {
    using namespace primerpair;

    try {
        NeuralWindowPredictor predictor;


        constexpr
        std::array<std::string_view, 8>
        kmers{
            "AAAAAAAAAAAAAAAAAAAAA",
            "ACGTACGTACGTACGTACGTA",
            "CGCGGAAGCAAAGTGACTTCC",
            "GAAATATAGGTATCAACGGAG",
            "CTGAATGGAATTCCTCCGATC",
            "CAAATGACGATGTCCTTGGGT",
            "GGGTTTTTTTTACACACACGT",
            "TGCATGCATGCATGCATGCAT"
        };


        constexpr
        std::uint64_t row_count =
            1000001;


        const auto batch =
            predictor.predict_rows_batch8(
                kmers,
                row_count
            );


        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            const std::uint64_t scalar =
                predictor.predict_row(
                    kmers[lane],
                    row_count
                );


            const std::uint64_t difference =
                (
                    scalar >
                    batch[lane]
                )
                ?
                scalar -
                    batch[lane]
                :
                batch[lane] -
                    scalar;


            /*
             * Floating-point vectorization may alter
             * the last bit of accumulation.
             *
             * At most one SA row is acceptable.
             */
            if (
                difference >
                1
            ) {
                std::cerr
                    << "BATCH8_MISMATCH"
                    << '\t'
                    << lane
                    << '\t'
                    << scalar
                    << '\t'
                    << batch[lane]
                    << '\t'
                    << difference
                    << '\n';

                throw std::runtime_error(
                    "Batch8 predictor differs "
                    "from scalar predictor."
                );
            }
        }


        std::cout
            << "batch_size\t8\n";

#if defined(PRIMERPAIR_HAVE_AVX2_BACKEND) && PRIMERPAIR_HAVE_AVX2_BACKEND

        std::cout
            << "backend\tAVX2\n";

#else

        std::cout
            << "backend\tSCALAR_FALLBACK\n";

#endif

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
