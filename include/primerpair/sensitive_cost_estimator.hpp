#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "primerpair/bidirectional_fm_index.hpp"

namespace primerpair {

struct SensitiveCandidateCostEstimate {
    std::uint64_t forward_left_occurrences{0};
    std::uint64_t forward_right_occurrences{0};

    std::uint64_t reverse_left_occurrences{0};
    std::uint64_t reverse_right_occurrences{0};

    /*
     * Sum before candidate-start deduplication.
     *
     * For the k=3 half-seed backend this closely
     * describes how much locate/verification work
     * may be generated.
     */
    std::uint64_t total_seed_occurrences{0};

    std::uint64_t max_seed_occurrences{0};

    bool operator==(
        const SensitiveCandidateCostEstimate&
    ) const = default;
};

class SensitiveCandidateCostEstimator {
public:
    explicit SensitiveCandidateCostEstimator(
        const BidirectionalFMIndex& index
    );

    /*
     * Estimate the current 20-mer / k=3 half-seed
     * strategy.
     *
     * Primer is split into two contiguous halves.
     * For both forward and reverse-complement query,
     * every seed sequence within Hamming distance <=1
     * is counted by FM interval size WITHOUT locate().
     */
    [[nodiscard]]
    SensitiveCandidateCostEstimate estimate_k3(
        std::string_view primer
    ) const;

private:
    [[nodiscard]]
    std::uint64_t count_seed_hamming1(
        std::string_view seed
    ) const;

    const BidirectionalFMIndex& index_;
};

}  // namespace primerpair
