#include "primerpair/bidirectional_fm_index.hpp"

#include <cctype>
#include <stdexcept>
#include <utility>

namespace primerpair {

BidirectionalFMIndex::BidirectionalFMIndex(
    std::string text,
    const std::size_t suffix_array_sample_rate
)
    : forward_index_(
          text,
          suffix_array_sample_rate
      ),
      reverse_index_(
          reverse_string(text),
          suffix_array_sample_rate
      ) {

    if (
        forward_index_.bwt_size() !=
        reverse_index_.bwt_size()
    ) {
        throw std::logic_error(
            "Forward and reverse FM-index lengths differ."
        );
    }
}

char BidirectionalFMIndex::normalize_base(
    const char base
) {
    const char normalized =
        static_cast<char>(
            std::toupper(
                static_cast<unsigned char>(
                    base
                )
            )
        );

    switch (normalized) {
        case 'A':
        case 'C':
        case 'G':
        case 'N':
        case 'T':
            return normalized;

        default:
            throw std::invalid_argument(
                "Invalid bidirectional-search nucleotide."
            );
    }
}

std::string BidirectionalFMIndex::normalize_pattern(
    const std::string_view pattern
) {
    if (pattern.empty()) {
        throw std::invalid_argument(
            "Bidirectional search pattern cannot be empty."
        );
    }

    std::string normalized;
    normalized.reserve(
        pattern.size()
    );

    for (const char base : pattern) {
        normalized.push_back(
            normalize_base(base)
        );
    }

    return normalized;
}

std::string BidirectionalFMIndex::reverse_string(
    const std::string_view text
) {
    return std::string(
        text.rbegin(),
        text.rend()
    );
}

std::size_t BidirectionalFMIndex::alphabet_index(
    const char base
) {
    switch (base) {
        case '$':
            return 0;

        case 'A':
            return 1;

        case 'C':
            return 2;

        case 'G':
            return 3;

        case 'N':
            return 4;

        case 'T':
            return 5;

        default:
            throw std::invalid_argument(
                "Base is not part of the BiFM alphabet."
            );
    }
}

std::uint64_t BidirectionalFMIndex::prefix_less_than(
    const std::array<std::uint64_t, 6>& counts,
    const char base
) {
    const std::size_t target =
        alphabet_index(
            base
        );

    std::uint64_t total = 0;

    for (std::size_t i = 0;
         i < target;
         ++i) {

        total +=
            counts.at(i);
    }

    return total;
}

void BidirectionalFMIndex::validate_pair(
    const Interval& forward,
    const Interval& reverse
) {
    if (
        forward.size() !=
        reverse.size()
    ) {
        throw std::logic_error(
            "Forward and reverse BiFM intervals "
            "represent different occurrence counts."
        );
    }
}

BidirectionalInterval BidirectionalFMIndex::search(
    const std::string_view pattern
) const {
    const std::string normalized =
        normalize_pattern(
            pattern
        );

    const std::string reversed =
        reverse_string(
            normalized
        );

    const Interval forward =
        forward_index_.backward_search(
            normalized
        );

    const Interval reverse =
        reverse_index_.backward_search(
            reversed
        );

    validate_pair(
        forward,
        reverse
    );

    return BidirectionalInterval{
        forward,
        reverse,
        normalized.size()
    };
}

BidirectionalInterval BidirectionalFMIndex::extend_left(
    const BidirectionalInterval& state,
    const char base
) const {
    if (state.length == 0) {
        throw std::invalid_argument(
            "Cannot extend an empty BiFM state."
        );
    }

    validate_pair(
        state.forward,
        state.reverse
    );

    const char normalized =
        normalize_base(
            base
        );

    /*
     * Önce native forward extension.
     *
     * Eğer sonuç boşsa karşı interval için rank/count
     * hesaplamaya gerek yoktur. İki tarafı da kanonik
     * empty interval {0,0} olarak döndürürüz.
     */
    const Interval forward =
        forward_index_.backward_extend(
            state.forward,
            normalized
        );

    if (forward.empty()) {
        return BidirectionalInterval{
            Interval{},
            Interval{},
            state.length + 1
        };
    }

    /*
     * I(P) içindeki predecessor dağılımı:
     *
     * $P, AP, CP, GP, NP, TP
     */
    const auto counts =
        forward_index_
            .interval_symbol_counts(
                state.forward
            );

    const std::uint64_t offset =
        prefix_less_than(
            counts,
            normalized
        );

    const std::uint64_t reverse_begin =
        state.reverse.begin +
        offset;

    const Interval reverse{
        reverse_begin,
        reverse_begin +
            forward.size()
    };

    validate_pair(
        forward,
        reverse
    );

    return BidirectionalInterval{
        forward,
        reverse,
        state.length + 1
    };
}


std::array<BidirectionalInterval, 4>
BidirectionalFMIndex::extend_left_all(
    const BidirectionalInterval& state
) const {
    if (state.length == 0) {
        throw std::invalid_argument(
            "Cannot extend an empty BiFM state."
        );
    }

    validate_pair(
        state.forward,
        state.reverse
    );

    /*
     * Lexicographic FM alphabet:
     *
     *   0=$
     *   1=A
     *   2=C
     *   3=G
     *   4=N
     *   5=T
     *
     * Tek rank_all(begin/end) çiftiyle bütün
     * forward child interval'ları üretilir.
     */
    const auto forward_children =
        forward_index_
            .backward_extend_all(
                state.forward
            );

    std::array<
        BidirectionalInterval,
        4
    > result{};

    std::uint64_t reverse_offset = 0;

    std::size_t output_index = 0;

    for (
        std::size_t alphabet_index = 0;
        alphabet_index <
            forward_children.size();
        ++alphabet_index
    ) {
        const Interval& forward =
            forward_children.at(
                alphabet_index
            );

        /*
         * DirectBranching yalnız A/C/G/T üretir.
         *
         * $ ve N dönüş array'ine girmez; ancak
         * lexicographic reverse offset'e katkıları
         * korunmalıdır.
         */
        const bool requested =
            alphabet_index == 1 ||
            alphabet_index == 2 ||
            alphabet_index == 3 ||
            alphabet_index == 5;

        if (requested) {

            if (forward.empty()) {
                result.at(
                    output_index
                ) = BidirectionalInterval{
                    Interval{},
                    Interval{},
                    state.length + 1
                };

            } else {
                const std::uint64_t
                    reverse_begin =
                        state.reverse.begin +
                        reverse_offset;

                const Interval reverse{
                    reverse_begin,
                    reverse_begin +
                        forward.size()
                };

                validate_pair(
                    forward,
                    reverse
                );

                result.at(
                    output_index
                ) = BidirectionalInterval{
                    forward,
                    reverse,
                    state.length + 1
                };
            }

            ++output_index;
        }

        /*
         * Burada $ ve N dahil bütün sembollerin
         * occurrence sayısı offset'e eklenir.
         *
         * T child için N'nin de prefix'te bulunması
         * özellikle önemlidir.
         */
        reverse_offset +=
            forward.size();
    }

    if (
        output_index !=
        result.size()
    ) {
        throw std::logic_error(
            "extend_left_all did not produce "
            "four A/C/G/T children."
        );
    }

    return result;
}


BidirectionalInterval BidirectionalFMIndex::extend_right(
    const BidirectionalInterval& state,
    const char base
) const {
    if (state.length == 0) {
        throw std::invalid_argument(
            "Cannot extend an empty BiFM state."
        );
    }

    validate_pair(
        state.forward,
        state.reverse
    );

    const char normalized =
        normalize_base(
            base
        );

    /*
     * reverse(Pc) = c reverse(P)
     *
     * Önce reverse index üzerinde native backward
     * extension yapılır.
     */
    const Interval reverse =
        reverse_index_.backward_extend(
            state.reverse,
            normalized
        );

    /*
     * Başarısız branch:
     * iki intervali de kanonik {0,0} yap.
     *
     * Özellikle ileride mismatch branching sırasında
     * ölen dallarda gereksiz rank hesaplarını önler.
     */
    if (reverse.empty()) {
        return BidirectionalInterval{
            Interval{},
            Interval{},
            state.length + 1
        };
    }

    /*
     * reverse(P) intervalindeki predecessor dağılımı,
     * forward metinde P'nin right-extension
     * dağılımına karşılık gelir.
     */
    const auto counts =
        reverse_index_
            .interval_symbol_counts(
                state.reverse
            );

    const std::uint64_t offset =
        prefix_less_than(
            counts,
            normalized
        );

    const std::uint64_t forward_begin =
        state.forward.begin +
        offset;

    const Interval forward{
        forward_begin,
        forward_begin +
            reverse.size()
    };

    validate_pair(
        forward,
        reverse
    );

    return BidirectionalInterval{
        forward,
        reverse,
        state.length + 1
    };
}

std::vector<std::uint64_t>
BidirectionalFMIndex::locate_unsorted(
    const BidirectionalInterval& state
) const {
    validate_pair(
        state.forward,
        state.reverse
    );

    return
        forward_index_.locate_unsorted(
            state.forward
        );
}


std::vector<std::uint64_t>
BidirectionalFMIndex::locate(
    const BidirectionalInterval& state
) const {
    validate_pair(
        state.forward,
        state.reverse
    );

    return forward_index_.locate(
        state.forward
    );
}

}  // namespace primerpair
