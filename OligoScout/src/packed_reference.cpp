#include "primerpair/packed_reference.hpp"

#include <cctype>
#include <limits>
#include <stdexcept>

namespace primerpair {

namespace {

char normalize_reference_base(
    const char raw
) {
    const char base =
        static_cast<char>(
            std::toupper(
                static_cast<unsigned char>(
                    raw
                )
            )
        );

    switch (base) {

        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'N':
            return base;

        /*
         * Standard IUPAC ambiguity codes.
         *
         * Exact candidate verification açısından
         * bunları N/unknown olarak değerlendiriyoruz.
         */
        case 'R':
        case 'Y':
        case 'S':
        case 'W':
        case 'K':
        case 'M':
        case 'B':
        case 'D':
        case 'H':
        case 'V':
            return 'N';

        default:
            throw std::invalid_argument(
                "Unsupported reference nucleotide."
            );
    }
}

char normalize_query_base(
    const char raw
) {
    const char base =
        static_cast<char>(
            std::toupper(
                static_cast<unsigned char>(
                    raw
                )
            )
        );

    switch (base) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return base;

        default:
            throw std::invalid_argument(
                "Hamming query must contain "
                "only A/C/G/T."
            );
    }
}

}  // namespace

PackedReference::PackedReference(
    const std::string_view sequence
) {
    if (sequence.empty()) {
        throw std::invalid_argument(
            "Packed reference cannot be empty."
        );
    }

    if (
        sequence.size() >
        std::numeric_limits<
            std::uint64_t
        >::max()
    ) {
        throw std::length_error(
            "Reference sequence is too large."
        );
    }

    length_ =
        static_cast<std::uint64_t>(
            sequence.size()
        );

    const std::size_t word_count =
        static_cast<std::size_t>(
            (length_ + 63ULL) /
            64ULL
        );

    low_bits_.assign(
        word_count,
        0
    );

    high_bits_.assign(
        word_count,
        0
    );

    n_mask_.assign(
        word_count,
        0
    );

    for (std::uint64_t i = 0;
         i < length_;
         ++i) {

        const char base =
            normalize_reference_base(
                sequence.at(
                    static_cast<std::size_t>(
                        i
                    )
                )
            );

        /*
         * Encoding:
         *
         * A 00
         * C 01
         * G 10
         * T 11
         */
        switch (base) {

            case 'A':
                break;

            case 'C':
                set_bit(
                    low_bits_,
                    i
                );
                break;

            case 'G':
                set_bit(
                    high_bits_,
                    i
                );
                break;

            case 'T':
                set_bit(
                    low_bits_,
                    i
                );

                set_bit(
                    high_bits_,
                    i
                );
                break;

            case 'N':
                set_bit(
                    n_mask_,
                    i
                );
                break;

            default:
                throw std::logic_error(
                    "Unexpected normalized base."
                );
        }
    }
}

void PackedReference::set_bit(
    std::vector<std::uint64_t>& words,
    const std::uint64_t position
) {
    const std::size_t word =
        static_cast<std::size_t>(
            position >> 6
        );

    const std::uint64_t bit =
        position & 63ULL;

    words.at(word) |=
        (
            std::uint64_t{1} <<
            bit
        );
}

bool PackedReference::get_bit(
    const std::vector<std::uint64_t>& words,
    const std::uint64_t position
) noexcept {
    const std::size_t word =
        static_cast<std::size_t>(
            position >> 6
        );

    const std::uint64_t bit =
        position & 63ULL;

    return
        (
            words[word] >>
            bit
        ) &
        std::uint64_t{1};
}

char PackedReference::base_at(
    const std::uint64_t position
) const {
    if (position >= length_) {
        throw std::out_of_range(
            "Packed reference position "
            "out of range."
        );
    }

    if (
        get_bit(
            n_mask_,
            position
        )
    ) {
        return 'N';
    }

    const bool low =
        get_bit(
            low_bits_,
            position
        );

    const bool high =
        get_bit(
            high_bits_,
            position
        );

    if (!high && !low) {
        return 'A';
    }

    if (!high && low) {
        return 'C';
    }

    if (high && !low) {
        return 'G';
    }

    return 'T';
}

std::size_t
PackedReference::bounded_hamming_distance(
    const std::uint64_t start,
    const std::string_view query,
    const std::size_t max_mismatches
) const {
    if (query.empty()) {
        throw std::invalid_argument(
            "Hamming query cannot be empty."
        );
    }

    /*
     * Overflow-safe bounds test.
     */
    if (
        start > length_ ||
        query.size() >
            length_ - start
    ) {
        throw std::out_of_range(
            "Hamming query exceeds "
            "reference bounds."
        );
    }

    std::size_t mismatches = 0;

    for (std::size_t i = 0;
         i < query.size();
         ++i) {

        const char expected =
            normalize_query_base(
                query.at(i)
            );

        const char observed =
            base_at(
                start +
                static_cast<std::uint64_t>(
                    i
                )
            );

        /*
         * Unknown reference sequence is not accepted
         * as a biological primer-binding candidate.
         */
        if (observed == 'N') {

            if (
                max_mismatches ==
                std::numeric_limits<
                    std::size_t
                >::max()
            ) {
                return max_mismatches;
            }

            return
                max_mismatches + 1;
        }

        if (observed != expected) {
            ++mismatches;

            if (
                mismatches >
                max_mismatches
            ) {
                return mismatches;
            }
        }
    }

    return mismatches;
}

std::size_t
PackedReference::memory_bytes() const noexcept {
    return
        low_bits_.size() *
            sizeof(std::uint64_t) +
        high_bits_.size() *
            sizeof(std::uint64_t) +
        n_mask_.size() *
            sizeof(std::uint64_t);
}

}  // namespace primerpair
