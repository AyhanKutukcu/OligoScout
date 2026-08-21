#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"

namespace primerpair {

struct AnchorSearchResult {
    BidirectionalInterval state{};

    std::size_t primer_length{0};
    std::size_t anchor_length{0};
    std::size_t extension_steps{0};

    [[nodiscard]]
    bool empty() const noexcept {
        return state.empty();
    }

    [[nodiscard]]
    std::uint64_t match_count() const noexcept {
        return state.size();
    }
};

class AnchorSearcher {
public:
    explicit AnchorSearcher(
        const BidirectionalFMIndex& index
    ) noexcept;

    [[nodiscard]]
    AnchorSearchResult search_exact(
        std::string_view primer,
        std::size_t anchor_length = 12
    ) const;

    [[nodiscard]]
    std::vector<std::uint64_t> locate(
        const AnchorSearchResult& result
    ) const;

private:
    const BidirectionalFMIndex& index_;
};

}  // namespace primerpair
