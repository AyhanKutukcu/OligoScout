// BIOLOGICAL_MISMATCH_PROFILE_V1
#include "primerpair/biological_mismatch_profile.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {

bool approximately_equal(
    const double lhs,
    const double rhs,
    const double tolerance = 1e-12
) {
    return
        std::abs(
            lhs - rhs
        )
        <=
        tolerance;
}


void require(
    const bool condition,
    const std::string_view message
) {
    if (!condition) {
        throw std::runtime_error(
            std::string(
                message
            )
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
        constexpr std::size_t length = 20;

        /*
         * -------------------------------------------------
         * Exact match
         * -------------------------------------------------
         */
        const auto exact =
            build_biological_mismatch_profile(
                0,
                length
            );

        require(
            exact.primer_length == length,
            "Exact profile primer length"
        );

        require(
            exact.mismatch_count == 0,
            "Exact profile mismatch count"
        );

        require(
            approximately_equal(
                exact.mismatch_fraction,
                0.0
            ),
            "Exact profile mismatch fraction"
        );

        require(
            exact.last_1_count == 0 &&
            exact.last_2_count == 0 &&
            exact.last_3_count == 0 &&
            exact.last_5_count == 0 &&
            exact.last_8_count == 0 &&
            exact.last_12_count == 0,
            "Exact profile 3-prime windows"
        );

        require(
            !exact.terminal_3prime_mismatch,
            "Exact profile terminal 3-prime state"
        );

        require(
            !exact
                .nearest_mismatch_to_3prime
                .has_value(),
            "Exact profile nearest mismatch absent"
        );

        require(
            exact.exact_3prime_run_length ==
                length,
            "Exact profile full 3-prime exact run"
        );

        require(
            approximately_equal(
                exact
                    .normalized_3prime_positional_burden,
                0.0
            ),
            "Exact profile positional burden zero"
        );


        /*
         * -------------------------------------------------
         * Single 5-prime-most mismatch
         * -------------------------------------------------
         */
        const auto five_prime =
            build_biological_mismatch_profile(
                std::uint64_t{1},
                length
            );

        require(
            five_prime.mismatch_count == 1,
            "5-prime mismatch count"
        );

        require(
            five_prime
                .nearest_mismatch_to_3prime
                .has_value() &&
            five_prime
                .nearest_mismatch_to_3prime
                .value() == 19,
            "5-prime mismatch distance from 3-prime"
        );

        require(
            five_prime.exact_3prime_run_length ==
                19,
            "5-prime mismatch leaves 19-base exact 3-prime run"
        );

        require(
            five_prime.last_12_count == 0,
            "5-prime mismatch outside terminal 12 bases"
        );


        /*
         * -------------------------------------------------
         * Terminal 3-prime mismatch
         * -------------------------------------------------
         */
        const std::uint64_t terminal_mask =
            std::uint64_t{1}
            <<
            19;

        const auto terminal =
            build_biological_mismatch_profile(
                terminal_mask,
                length
            );

        require(
            terminal.mismatch_count == 1,
            "Terminal mismatch count"
        );

        require(
            terminal.terminal_3prime_mismatch,
            "Terminal 3-prime mismatch detected"
        );

        require(
            terminal.last_1_count == 1,
            "Terminal mismatch represented in last-1 window"
        );

        require(
            terminal
                .nearest_mismatch_to_3prime
                .has_value() &&
            terminal
                .nearest_mismatch_to_3prime
                .value() == 0,
            "Terminal mismatch distance equals zero"
        );

        require(
            terminal.exact_3prime_run_length == 0,
            "Terminal mismatch gives zero-length exact 3-prime run"
        );

        require(
            terminal
                .normalized_3prime_positional_burden
            >
            five_prime
                .normalized_3prime_positional_burden,
            "Terminal mismatch burden exceeds 5-prime mismatch"
        );


        /*
         * -------------------------------------------------
         * Monotonic 3-prime positional weighting
         *
         * Move one mismatch one nucleotide at a time
         * from 5-prime toward 3-prime.
         * Burden must strictly increase.
         * -------------------------------------------------
         */
        double previous_burden = -1.0;

        for (
            std::size_t position = 0;
            position < length;
            ++position
        ) {
            const std::uint64_t mask =
                std::uint64_t{1}
                <<
                position;

            const auto profile =
                build_biological_mismatch_profile(
                    mask,
                    length
                );

            require(
                profile
                    .normalized_3prime_positional_burden
                >
                previous_burden,
                "Single-mismatch burden increases toward 3-prime"
            );

            previous_burden =
                profile
                    .normalized_3prime_positional_burden;
        }


        /*
         * -------------------------------------------------
         * Multi-mismatch window semantics
         *
         * length = 20
         *
         * positions:
         *
         * 0  = first 5-prime base
         * 12 = distance 7 from 3-prime
         * 17 = distance 2 from 3-prime
         * 19 = terminal 3-prime
         * -------------------------------------------------
         */
        const std::uint64_t mixed_mask =
            (
                std::uint64_t{1}
                <<
                0
            )
            |
            (
                std::uint64_t{1}
                <<
                12
            )
            |
            (
                std::uint64_t{1}
                <<
                17
            )
            |
            (
                std::uint64_t{1}
                <<
                19
            );

        const auto mixed =
            build_biological_mismatch_profile(
                mixed_mask,
                length
            );

        require(
            mixed.mismatch_count == 4,
            "Mixed profile total mismatch count"
        );

        require(
            approximately_equal(
                mixed.mismatch_fraction,
                0.2
            ),
            "Mixed profile mismatch fraction"
        );

        require(
            mixed.last_1_count == 1,
            "Mixed profile last-1 count"
        );

        require(
            mixed.last_2_count == 1,
            "Mixed profile last-2 count"
        );

        require(
            mixed.last_3_count == 2,
            "Mixed profile last-3 count"
        );

        require(
            mixed.last_5_count == 2,
            "Mixed profile last-5 count"
        );

        require(
            mixed.last_8_count == 3,
            "Mixed profile last-8 count"
        );

        require(
            mixed.last_12_count == 3,
            "Mixed profile last-12 count"
        );

        require(
            mixed.exact_3prime_run_length == 0,
            "Mixed profile terminal mismatch breaks exact 3-prime run"
        );


        /*
         * -------------------------------------------------
         * Full mismatch mask normalization
         * -------------------------------------------------
         */
        const std::uint64_t full_mask =
            (
                std::uint64_t{1}
                <<
                length
            )
            -
            1;

        const auto full =
            build_biological_mismatch_profile(
                full_mask,
                length
            );

        require(
            approximately_equal(
                full
                    .normalized_3prime_positional_burden,
                1.0
            ),
            "Fully mismatched primer has normalized burden one"
        );


        /*
         * -------------------------------------------------
         * Pair-level aggregation
         * -------------------------------------------------
         */
        PrimerPairHit pair_hit{};

        pair_hit.left_mismatch_mask =
            std::uint64_t{1}
            <<
            19;

        pair_hit.right_mismatch_mask =
            (
                std::uint64_t{1}
                <<
                18
            )
            |
            (
                std::uint64_t{1}
                <<
                19
            );

        pair_hit.left_mismatches = 1;
        pair_hit.right_mismatches = 2;

        const auto pair_profile =
            build_pair_biological_mismatch_profile(
                pair_hit,
                20,
                20
            );

        require(
            pair_profile.total_mismatches == 3,
            "Pair profile total mismatch count"
        );

        require(
            pair_profile
                .left
                .mismatch_count == 1 &&
            pair_profile
                .right
                .mismatch_count == 2,
            "Pair profile side-specific mismatch counts"
        );

        require(
            pair_profile
                .max_normalized_3prime_positional_burden
            >=
            pair_profile
                .mean_normalized_3prime_positional_burden,
            "Pair profile maximum burden >= mean burden"
        );


        /*
         * -------------------------------------------------
         * Existing PrimerPairHit invariant remains active.
         * -------------------------------------------------
         */
        bool count_disagreement_rejected = false;

        try {
            PrimerPairHit invalid_pair{};

            invalid_pair.left_mismatch_mask =
                std::uint64_t{1};

            invalid_pair.right_mismatch_mask = 0;

            invalid_pair.left_mismatches = 0;
            invalid_pair.right_mismatches = 0;

            static_cast<void>(
                build_pair_biological_mismatch_profile(
                    invalid_pair,
                    20,
                    20
                )
            );

        } catch (const std::logic_error&) {
            count_disagreement_rejected = true;
        }

        require(
            count_disagreement_rejected,
            "Pair mismatch count/mask disagreement rejected"
        );


        /*
         * -------------------------------------------------
         * Invalid primer lengths / masks rejected by
         * the validated mismatch feature layer.
         * -------------------------------------------------
         */
        bool zero_length_rejected = false;

        try {
            static_cast<void>(
                build_biological_mismatch_profile(
                    0,
                    0
                )
            );
        } catch (const std::invalid_argument&) {
            zero_length_rejected = true;
        }

        require(
            zero_length_rejected,
            "Zero primer length rejected"
        );


        bool outside_mask_rejected = false;

        try {
            static_cast<void>(
                build_biological_mismatch_profile(
                    std::uint64_t{1}
                        <<
                        20,
                    20
                )
            );
        } catch (const std::invalid_argument&) {
            outside_mask_rejected = true;
        }

        require(
            outside_mask_rejected,
            "Mismatch bits outside primer length rejected"
        );


        std::cout
            << "profile_primer_length\t"
            << length
            << '\n';

        std::cout
            << "terminal_burden\t"
            << terminal
                .normalized_3prime_positional_burden
            << '\n';

        std::cout
            << "five_prime_burden\t"
            << five_prime
                .normalized_3prime_positional_burden
            << '\n';

        std::cout
            << "mixed_mismatch_fraction\t"
            << mixed.mismatch_fraction
            << '\n';

        std::cout
            << "BIOLOGICAL_PROFILE_3PRIME_MONOTONIC\tYES\n";

        std::cout
            << "BIOLOGICAL_PROFILE_WINDOWS_VALID\tYES\n";

        std::cout
            << "BIOLOGICAL_PROFILE_PAIR_AGGREGATION\tYES\n";

        std::cout
            << "BIOLOGICAL_PROFILE_CALIBRATION_READY\tYES\n";

        std::cout
            << "ALL_CHECKS\tYES\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';

        return 1;
    }
}
