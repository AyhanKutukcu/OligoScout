#include "primerpair/packed_bwt.hpp"

#include <algorithm>
#include <bit>
#include <cctype>
#include <limits>
#include <stdexcept>
#include <string>

namespace primerpair {

PackedBWT::PackedBWT(
    const std::string_view bwt
) {
    build(bwt);
}

std::uint8_t PackedBWT::encode_acgt(
    const char symbol
) {
    const char normalized =
        static_cast<char>(
            std::toupper(
                static_cast<unsigned char>(
                    symbol
                )
            )
        );

    switch (normalized) {
        case 'A':
            return 0;

        case 'C':
            return 1;

        case 'G':
            return 2;

        case 'T':
            return 3;

        default:
            throw std::invalid_argument(
                std::string(
                    "Unsupported ACGT symbol: "
                ) +
                normalized
            );
    }
}

void PackedBWT::build(
    const std::string_view bwt
) {
    if (bwt.empty()) {
        throw std::invalid_argument(
            "BWT cannot be empty."
        );
    }

    if (
        bwt.size() >
        static_cast<std::size_t>(
            std::numeric_limits<
                std::uint32_t
            >::max()
        )
    ) {
        throw std::length_error(
            "BWT is too large for 32-bit sentinel position."
        );
    }

    length_ =
        static_cast<std::uint64_t>(
            bwt.size()
        );

    sentinel_position_ = 0;

    const std::size_t word_count =
        (bwt.size() + 63) / 64;

    low_bits_.assign(
        word_count,
        std::uint64_t{0}
    );

    high_bits_.assign(
        word_count,
        std::uint64_t{0}
    );

    n_mask_words_.assign(
        word_count,
        std::uint64_t{0}
    );

    bool sentinel_found = false;

    for (std::size_t position = 0;
         position < bwt.size();
         ++position) {

        const char normalized =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        bwt.at(position)
                    )
                )
            );

        const std::size_t word =
            position / 64;

        const std::size_t offset =
            position % 64;

        const std::uint64_t bit =
            std::uint64_t{1}
            << offset;

        if (normalized == '$') {
            if (sentinel_found) {
                throw std::invalid_argument(
                    "BWT contains more than one sentinel."
                );
            }

            sentinel_found = true;

            sentinel_position_ =
                static_cast<std::uint32_t>(
                    position
                );

            continue;
        }

        if (normalized == 'N') {
            n_mask_words_.at(word) |= bit;
            continue;
        }

        const std::uint8_t code =
            encode_acgt(
                normalized
            );

        if ((code & 0x1U) != 0U) {
            low_bits_.at(word) |= bit;
        }

        if ((code & 0x2U) != 0U) {
            high_bits_.at(word) |= bit;
        }
    }

    if (!sentinel_found) {
        throw std::invalid_argument(
            "BWT does not contain a sentinel."
        );
    }
}

std::uint8_t PackedBWT::acgt_code_at(
    const std::uint64_t position
) const {
    if (position >= length_) {
        throw std::out_of_range(
            "Packed BWT position exceeds length."
        );
    }

    const std::size_t word =
        static_cast<std::size_t>(
            position / 64
        );

    const std::size_t offset =
        static_cast<std::size_t>(
            position % 64
        );

    const std::uint64_t bit =
        std::uint64_t{1}
        << offset;

    const std::uint8_t low =
        (low_bits_.at(word) & bit) != 0
            ? 1U
            : 0U;

    const std::uint8_t high =
        (high_bits_.at(word) & bit) != 0
            ? 2U
            : 0U;

    return static_cast<std::uint8_t>(
        low | high
    );
}

bool PackedBWT::is_n(
    const std::uint64_t position
) const {
    if (position >= length_) {
        throw std::out_of_range(
            "Packed BWT position exceeds length."
        );
    }

    const std::size_t word =
        static_cast<std::size_t>(
            position / 64
        );

    const std::size_t offset =
        static_cast<std::size_t>(
            position % 64
        );

    return (
        (
            n_mask_words_.at(word) >>
            offset
        ) &
        std::uint64_t{1}
    ) != 0;
}

std::size_t PackedBWT::symbol_id_at(
    const std::uint64_t position
) const {
    if (position >= length_) {
        throw std::out_of_range(
            "Packed BWT position exceeds length."
        );
    }

    if (
        position ==
        static_cast<std::uint64_t>(
            sentinel_position_
        )
    ) {
        return 0;
    }

    const std::size_t word =
        static_cast<std::size_t>(
            position >> 6U
        );

    const std::size_t offset =
        static_cast<std::size_t>(
            position & 63ULL
        );

    const std::uint64_t bit =
        std::uint64_t{1}
        << offset;

    if (
        (
            n_mask_words_.at(word) &
            bit
        ) != 0
    ) {
        return 4;
    }

    const std::size_t low =
        (
            low_bits_.at(word) &
            bit
        ) != 0
            ? 1U
            : 0U;

    const std::size_t high =
        (
            high_bits_.at(word) &
            bit
        ) != 0
            ? 2U
            : 0U;

    /*
     * Packed ACGT code:
     *
     * 00 -> A -> FM id 1
     * 01 -> C -> FM id 2
     * 10 -> G -> FM id 3
     * 11 -> T -> FM id 5
     */
    switch (low | high) {
        case 0:
            return 1;

        case 1:
            return 2;

        case 2:
            return 3;

        case 3:
            return 5;

        default:
            throw std::logic_error(
                "Invalid packed BWT symbol code."
            );
    }
}


char PackedBWT::at(
    const std::uint64_t position
) const {
    if (position >= length_) {
        throw std::out_of_range(
            "Packed BWT position exceeds length."
        );
    }

    if (
        position ==
        static_cast<std::uint64_t>(
            sentinel_position_
        )
    ) {
        return '$';
    }

    if (is_n(position)) {
        return 'N';
    }

    switch (
        acgt_code_at(position)
    ) {
        case 0:
            return 'A';

        case 1:
            return 'C';

        case 2:
            return 'G';

        case 3:
            return 'T';

        default:
            throw std::logic_error(
                "Invalid packed BWT code."
            );
    }
}


std::array<std::uint64_t, 6>
PackedBWT::count_all(
    const std::uint64_t begin,
    const std::uint64_t end
) const {
    if (begin > end) {
        throw std::invalid_argument(
            "Packed BWT count_all begin exceeds end."
        );
    }

    if (end > length_) {
        throw std::out_of_range(
            "Packed BWT count_all range exceeds length."
        );
    }

    std::array<std::uint64_t, 6> counts{};

    if (begin == end) {
        return counts;
    }

    const std::uint64_t sentinel =
        static_cast<std::uint64_t>(
            sentinel_position_
        );

    if (
        sentinel >= begin &&
        sentinel < end
    ) {
        counts[0] = 1;
    }

    std::uint64_t current = begin;

    while (current < end) {
        const std::size_t word =
            static_cast<std::size_t>(
                current / 64
            );

        const std::size_t offset =
            static_cast<std::size_t>(
                current % 64
            );

        const std::uint64_t word_end =
            std::min(
                end,
                (
                    static_cast<std::uint64_t>(
                        word
                    ) +
                    1
                ) *
                64
            );

        const std::size_t bit_count =
            static_cast<std::size_t>(
                word_end - current
            );

        std::uint64_t range_mask = 0;

        if (
            offset == 0 &&
            bit_count == 64
        ) {
            range_mask =
                std::numeric_limits<
                    std::uint64_t
                >::max();

        } else {
            range_mask =
                (
                    (
                        std::uint64_t{1}
                        << bit_count
                    ) -
                    1
                )
                << offset;
        }

        const std::uint64_t n_mask =
            n_mask_words_.at(word);

        const std::uint64_t n_bits =
            n_mask &
            range_mask;

        counts[4] +=
            static_cast<std::uint64_t>(
                std::popcount(
                    n_bits
                )
            );

        std::uint64_t valid_mask =
            range_mask &
            ~n_mask;

        /*
         * Sentinel bit-plane kodu 00'dır.
         * Bu nedenle A maskesinden çıkarılmalıdır.
         */
        if (
            sentinel >= current &&
            sentinel < word_end
        ) {
            valid_mask &=
                ~(
                    std::uint64_t{1}
                    <<
                    static_cast<std::size_t>(
                        sentinel % 64
                    )
                );
        }

        const std::uint64_t low =
            low_bits_.at(word);

        const std::uint64_t high =
            high_bits_.at(word);

        const std::uint64_t a_bits =
            (~low) &
            (~high) &
            valid_mask;

        const std::uint64_t c_bits =
            low &
            (~high) &
            valid_mask;

        const std::uint64_t g_bits =
            (~low) &
            high &
            valid_mask;

        const std::uint64_t t_bits =
            low &
            high &
            valid_mask;

        counts[1] +=
            static_cast<std::uint64_t>(
                std::popcount(a_bits)
            );

        counts[2] +=
            static_cast<std::uint64_t>(
                std::popcount(c_bits)
            );

        counts[3] +=
            static_cast<std::uint64_t>(
                std::popcount(g_bits)
            );

        counts[5] +=
            static_cast<std::uint64_t>(
                std::popcount(t_bits)
            );

        current = word_end;
    }

    return counts;
}


std::uint64_t PackedBWT::count_by_symbol_id(
    const std::size_t symbol_id,
    const std::uint64_t begin,
    const std::uint64_t end
) const {
    if (begin > end) {
        throw std::invalid_argument(
            "Packed BWT count begin exceeds end."
        );
    }

    if (end > length_) {
        throw std::out_of_range(
            "Packed BWT count range exceeds length."
        );
    }

    if (symbol_id >= 6) {
        throw std::invalid_argument(
            "Packed BWT symbol ID exceeds alphabet."
        );
    }

    if (begin == end) {
        return 0;
    }

    if (symbol_id == 0) {
        const std::uint64_t sentinel =
            static_cast<std::uint64_t>(
                sentinel_position_
            );

        return (
            sentinel >= begin &&
            sentinel < end
        )
            ? 1
            : 0;
    }

    std::uint64_t total = 0;
    std::uint64_t current = begin;

    const std::uint64_t sentinel =
        static_cast<std::uint64_t>(
            sentinel_position_
        );

    while (current < end) {
        const std::size_t word =
            static_cast<std::size_t>(
                current >> 6U
            );

        const std::size_t offset =
            static_cast<std::size_t>(
                current & 63ULL
            );

        const std::uint64_t word_end =
            std::min(
                end,
                (
                    static_cast<std::uint64_t>(
                        word
                    ) +
                    1
                ) *
                64
            );

        const std::size_t bit_count =
            static_cast<std::size_t>(
                word_end -
                current
            );

        std::uint64_t range_mask = 0;

        if (
            offset == 0 &&
            bit_count == 64
        ) {
            range_mask =
                std::numeric_limits<
                    std::uint64_t
                >::max();
        } else {
            range_mask =
                (
                    (
                        std::uint64_t{1}
                        << bit_count
                    ) -
                    1
                )
                << offset;
        }

        const std::uint64_t n_mask =
            n_mask_words_.at(
                word
            );

        if (symbol_id == 4) {
            total +=
                static_cast<std::uint64_t>(
                    std::popcount(
                        n_mask &
                        range_mask
                    )
                );

            current = word_end;
            continue;
        }

        std::uint64_t valid_mask =
            range_mask &
            ~n_mask;

        /*
         * Sentinel bit-plane kodu 00 olduğundan A
         * maskesinden açık biçimde çıkarılır.
         */
        if (
            sentinel >= current &&
            sentinel < word_end
        ) {
            valid_mask &=
                ~(
                    std::uint64_t{1}
                    <<
                    static_cast<std::size_t>(
                        sentinel & 63ULL
                    )
                );
        }

        const std::uint64_t low =
            low_bits_.at(
                word
            );

        const std::uint64_t high =
            high_bits_.at(
                word
            );

        std::uint64_t matches = 0;

        switch (symbol_id) {
            case 1:  // A
                matches =
                    (~low) &
                    (~high);
                break;

            case 2:  // C
                matches =
                    low &
                    (~high);
                break;

            case 3:  // G
                matches =
                    (~low) &
                    high;
                break;

            case 5:  // T
                matches =
                    low &
                    high;
                break;

            default:
                throw std::logic_error(
                    "Invalid non-special packed BWT symbol ID."
                );
        }

        total +=
            static_cast<std::uint64_t>(
                std::popcount(
                    matches &
                    valid_mask
                )
            );

        current = word_end;
    }

    return total;
}


std::uint64_t PackedBWT::count(
    const char symbol,
    const std::uint64_t begin,
    const std::uint64_t end
) const {
    const char normalized =
        static_cast<char>(
            std::toupper(
                static_cast<unsigned char>(
                    symbol
                )
            )
        );

    std::size_t symbol_id = 0;

    switch (normalized) {
        case '$':
            symbol_id = 0;
            break;

        case 'A':
            symbol_id = 1;
            break;

        case 'C':
            symbol_id = 2;
            break;

        case 'G':
            symbol_id = 3;
            break;

        case 'N':
            symbol_id = 4;
            break;

        case 'T':
            symbol_id = 5;
            break;

        default:
            throw std::invalid_argument(
                std::string(
                    "Unsupported packed BWT count symbol: "
                ) +
                normalized
            );
    }

    return
        count_by_symbol_id(
            symbol_id,
            begin,
            end
        );
}


}  // namespace primerpair
