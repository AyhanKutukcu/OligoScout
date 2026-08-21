#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/fm_index.hpp"

namespace primerpair {

struct BidirectionalInterval {
    Interval forward{};
    Interval reverse{};

    std::size_t length{0};

    [[nodiscard]]
    bool empty() const noexcept {
        return forward.empty();
    }

    [[nodiscard]]
    std::uint64_t size() const noexcept {
        return forward.size();
    }

    [[nodiscard]]
    std::size_t matched_length() const noexcept {
        return length;
    }
};

class BidirectionalFMIndex {
    friend class PpfmIO;

public:
    explicit BidirectionalFMIndex(
        std::string text,
        std::size_t suffix_array_sample_rate =
            FMIndex::kDefaultSuffixArraySampleRate
    );

    [[nodiscard]]
    BidirectionalInterval search(
        std::string_view pattern
    ) const;

    [[nodiscard]]
    BidirectionalInterval extend_left(
        const BidirectionalInterval& state,
        char base
    ) const;


    /*
     * A, C, G, T left-extension child'larını
     * aynı parent BiFM state üzerinden toplu üretir.
     *
     * Dönüş sırası:
     *   0=A, 1=C, 2=G, 3=T
     */
    [[nodiscard]]
    std::array<BidirectionalInterval, 4>
    extend_left_all(
        const BidirectionalInterval& state
    ) const;

    [[nodiscard]]
    BidirectionalInterval extend_right(
        const BidirectionalInterval& state,
        char base
    ) const;

    [[nodiscard]]
    std::vector<std::uint64_t> locate(
        const BidirectionalInterval& state
    ) const;


    /*
     * locate() ile aynı koordinat kümesini döndürür
     * ancak sonuçları sıralamaz.
     */
    [[nodiscard]]
    std::vector<std::uint64_t> locate_unsorted(
        const BidirectionalInterval& state
    ) const;

    [[nodiscard]]
    const FMIndex& forward_index() const noexcept {
        return forward_index_;
    }

    [[nodiscard]]
    const FMIndex& reverse_index() const noexcept {
        return reverse_index_;
    }

private:
    static constexpr std::array<char, 6>
        kAlphabet{
            '$',
            'A',
            'C',
            'G',
            'N',
            'T'
        };

    FMIndex forward_index_;
    FMIndex reverse_index_;

    [[nodiscard]]
    static char normalize_base(
        char base
    );

    [[nodiscard]]
    static std::string normalize_pattern(
        std::string_view pattern
    );

    [[nodiscard]]
    static std::string reverse_string(
        std::string_view text
    );

    [[nodiscard]]
    static std::size_t alphabet_index(
        char base
    );

    [[nodiscard]]
    static std::uint64_t prefix_less_than(
        const std::array<std::uint64_t, 6>& counts,
        char base
    );

    static void validate_pair(
        const Interval& forward,
        const Interval& reverse
    );
};

}  // namespace primerpair
