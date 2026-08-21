#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/search_strategy.hpp"

namespace primerpair {

enum class AnchorPlacement {
    Suffix,
    Prefix
};

struct BatchedAnchorRequest {
    std::string_view query{};

    AnchorPlacement placement{
        AnchorPlacement::Suffix
    };

    std::size_t max_mismatches{3};
};

struct BatchedAnchorDecision {
    Interval interval{};

    std::uint64_t occurrences{0};

    SearchDifficulty difficulty{
        SearchDifficulty::Easy
    };

    SearchStrategy strategy{
        SearchStrategy::DirectBranching
    };

    std::size_t query_length{0};
    std::size_t anchor_length{0};
};

class BatchedAnchorLookup {
public:
    BatchedAnchorLookup(
        const IPBWTIndex& index,
        const SearchDifficultyEstimator& estimator
    ) noexcept;

    [[nodiscard]]
    std::vector<BatchedAnchorDecision>
    lookup(
        const std::vector<BatchedAnchorRequest>& requests,
        std::size_t anchor_length = 12
    ) const;

private:
    const IPBWTIndex& index_;

    const SearchDifficultyEstimator&
        estimator_;
};

}  // namespace primerpair
