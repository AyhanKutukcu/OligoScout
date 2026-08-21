#include "primerpair/sensitive_cost_estimator.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {

namespace {

constexpr std::array<char, 4>
    kBases{
        'A',
        'C',
        'G',
        'T'
    };

std::string normalize(
    const std::string_view primer
) {
    if (
        primer.size() < 2 ||
        primer.size() > 64
    ) {
        throw std::invalid_argument(
            "Primer length for sensitive-cost "
            "estimation must be in range 2..64."
        );
    }

    std::string output;

    output.reserve(
        primer.size()
    );

    for (const char raw : primer) {

        const char base =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        raw
                    )
                )
            );

        if (
            base != 'A' &&
            base != 'C' &&
            base != 'G' &&
            base != 'T'
        ) {
            throw std::invalid_argument(
                "Primer must contain only A/C/G/T."
            );
        }

        output.push_back(
            base
        );
    }

    return output;
}

}  // namespace

SensitiveCandidateCostEstimator::
SensitiveCandidateCostEstimator(
    const BidirectionalFMIndex& index
)
    : index_(
          index
      ) {
}

std::uint64_t
SensitiveCandidateCostEstimator::
count_seed_hamming1(
    const std::string_view seed
) const {
    if (seed.empty()) {
        throw std::invalid_argument(
            "Cost-estimator seed cannot be empty."
        );
    }

    std::uint64_t total = 0;

    /*
     * Exact seed.
     */
    total +=
        index_.search(
            seed
        ).size();

    /*
     * Every sequence exactly one mismatch away.
     *
     * For an L-mer there are 3*L such sequences.
     * They are mutually distinct, so interval sizes
     * can be summed directly.
     */
    std::string mutated(
        seed
    );

    for (
        std::size_t position = 0;
        position < seed.size();
        ++position
    ) {
        const char original =
            seed.at(
                position
            );

        for (const char base : kBases) {

            if (base == original) {
                continue;
            }

            mutated.at(
                position
            ) = base;

            total +=
                index_.search(
                    mutated
                ).size();
        }

        mutated.at(
            position
        ) = original;
    }

    return total;
}

SensitiveCandidateCostEstimate
SensitiveCandidateCostEstimator::
estimate_k3(
    const std::string_view primer
) const {
    const std::string normalized =
        normalize(
            primer
        );

    const std::size_t split =
        normalized.size() / 2;

    if (
        split == 0 ||
        split ==
            normalized.size()
    ) {
        throw std::invalid_argument(
            "Primer cannot be split into two seeds."
        );
    }

    const std::string_view
        forward_left(
            normalized.data(),
            split
        );

    const std::string_view
        forward_right(
            normalized.data() + split,
            normalized.size() - split
        );

    const std::string reverse =
        reverse_complement(
            normalized
        );

    const std::string_view
        reverse_left(
            reverse.data(),
            split
        );

    const std::string_view
        reverse_right(
            reverse.data() + split,
            reverse.size() - split
        );

    SensitiveCandidateCostEstimate
        estimate;

    estimate.forward_left_occurrences =
        count_seed_hamming1(
            forward_left
        );

    estimate.forward_right_occurrences =
        count_seed_hamming1(
            forward_right
        );

    estimate.reverse_left_occurrences =
        count_seed_hamming1(
            reverse_left
        );

    estimate.reverse_right_occurrences =
        count_seed_hamming1(
            reverse_right
        );

    estimate.total_seed_occurrences =
        estimate.forward_left_occurrences +
        estimate.forward_right_occurrences +
        estimate.reverse_left_occurrences +
        estimate.reverse_right_occurrences;

    estimate.max_seed_occurrences =
        std::max({
            estimate.forward_left_occurrences,
            estimate.forward_right_occurrences,
            estimate.reverse_left_occurrences,
            estimate.reverse_right_occurrences
        });

    return estimate;
}

}  // namespace primerpair
