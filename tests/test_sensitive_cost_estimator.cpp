#include "primerpair/sensitive_cost_estimator.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "primerpair/strand_aware_primer_search.hpp"

namespace {

void expect(
    const bool condition,
    const std::string& name
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " + name
        );
    }

    std::cout
        << "[PASS] "
        << name
        << '\n';
}

std::size_t hamming(
    const std::string_view lhs,
    const std::string_view rhs
) {
    if (
        lhs.size() !=
        rhs.size()
    ) {
        throw std::runtime_error(
            "Hamming length mismatch."
        );
    }

    std::size_t distance = 0;

    for (
        std::size_t i = 0;
        i < lhs.size();
        ++i
    ) {
        if (
            lhs.at(i) !=
            rhs.at(i)
        ) {
            ++distance;
        }
    }

    return distance;
}

std::uint64_t brute_seed_count(
    const std::string& reference,
    const std::string_view seed
) {
    std::uint64_t total = 0;

    if (
        reference.size() <
        seed.size()
    ) {
        return 0;
    }

    for (
        std::size_t start = 0;
        start + seed.size() <=
            reference.size();
        ++start
    ) {
        const std::string_view candidate(
            reference.data() + start,
            seed.size()
        );

        bool canonical = true;

        for (const char base : candidate) {

            if (
                base != 'A' &&
                base != 'C' &&
                base != 'G' &&
                base != 'T'
            ) {
                canonical = false;
                break;
            }
        }

        if (!canonical) {
            continue;
        }

        if (
            hamming(
                candidate,
                seed
            ) <= 1
        ) {
            ++total;
        }
    }

    return total;
}

}  // namespace

int main() {
    try {
        const std::string primer =
            "ACGTTGCAAGTCCTGAACGA";

        const std::string reverse =
            primerpair::
                reverse_complement(
                    primer
                );

        /*
         * A deliberately repetitive synthetic reference.
         */
        std::string reference =
            "NNNN";

        reference += primer;
        reference += "ACGTTGCAAA";
        reference += "ACGTTGCAAG";
        reference += "TCCTGAACGA";
        reference += "TCCTGAACGG";

        reference += reverse;

        reference +=
            "NNNN";

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::
            SensitiveCandidateCostEstimator
                estimator(
                    index
                );

        const auto estimate =
            estimator.estimate_k3(
                primer
            );

        const std::size_t split =
            primer.size() / 2;

        const std::string
            forward_left =
                primer.substr(
                    0,
                    split
                );

        const std::string
            forward_right =
                primer.substr(
                    split
                );

        const std::string
            reverse_left =
                reverse.substr(
                    0,
                    split
                );

        const std::string
            reverse_right =
                reverse.substr(
                    split
                );

        const auto brute_fl =
            brute_seed_count(
                reference,
                forward_left
            );

        const auto brute_fr =
            brute_seed_count(
                reference,
                forward_right
            );

        const auto brute_rl =
            brute_seed_count(
                reference,
                reverse_left
            );

        const auto brute_rr =
            brute_seed_count(
                reference,
                reverse_right
            );

        expect(
            estimate.forward_left_occurrences ==
                brute_fl,
            "Forward-left Hamming<=1 count"
        );

        expect(
            estimate.forward_right_occurrences ==
                brute_fr,
            "Forward-right Hamming<=1 count"
        );

        expect(
            estimate.reverse_left_occurrences ==
                brute_rl,
            "Reverse-left Hamming<=1 count"
        );

        expect(
            estimate.reverse_right_occurrences ==
                brute_rr,
            "Reverse-right Hamming<=1 count"
        );

        expect(
            estimate.total_seed_occurrences ==
                brute_fl +
                brute_fr +
                brute_rl +
                brute_rr,
            "Total seed occurrence estimate"
        );

        expect(
            estimate.max_seed_occurrences ==
                std::max({
                    brute_fl,
                    brute_fr,
                    brute_rl,
                    brute_rr
                }),
            "Maximum seed occurrence estimate"
        );

        std::cout
            << "total_seed_occurrences\t"
            << estimate.total_seed_occurrences
            << '\n';

        std::cout
            << "max_seed_occurrences\t"
            << estimate.max_seed_occurrences
            << '\n';

        std::cout
            << "All sensitive cost-estimator tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
