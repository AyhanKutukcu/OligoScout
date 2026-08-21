#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace primerpair {

/*
 * Compact reference sequence representation.
 *
 * A = 00
 * C = 01
 * G = 10
 * T = 11
 *
 * Ambiguous/N bases are represented by an
 * independent N-mask.
 *
 * Approximate persistent storage:
 *
 * 2 bits/bp sequence
 * 1 bit/bp N mask
 * ----------------
 * 3 bits/bp = 0.375 bytes/bp
 */
class PackedReference {
    friend class PpfmIO;

public:
    explicit PackedReference(
        std::string_view sequence
    );

    [[nodiscard]]
    std::uint64_t size() const noexcept {
        return length_;
    }

    [[nodiscard]]
    bool empty() const noexcept {
        return length_ == 0;
    }

    /*
     * Returns A/C/G/T/N.
     */
    [[nodiscard]]
    char base_at(
        std::uint64_t position
    ) const;

    /*
     * Bounded Hamming verification.
     *
     * - Query must contain only A/C/G/T.
     * - Any N/ambiguous reference base rejects
     *   the candidate immediately.
     * - Stops once mismatches > max_mismatches.
     *
     * Return:
     *
     * 0..max_mismatches:
     *     valid Hamming distance
     *
     * max_mismatches + 1:
     *     candidate failed verification
     */
    [[nodiscard]]
    std::size_t bounded_hamming_distance(
        std::uint64_t start,
        std::string_view query,
        std::size_t max_mismatches
    ) const;

    [[nodiscard]]
    std::size_t memory_bytes() const noexcept;

private:
    std::uint64_t length_{0};

    std::vector<std::uint64_t>
        low_bits_;

    std::vector<std::uint64_t>
        high_bits_;

    std::vector<std::uint64_t>
        n_mask_;

    static void set_bit(
        std::vector<std::uint64_t>& words,
        std::uint64_t position
    );

    [[nodiscard]]
    static bool get_bit(
        const std::vector<std::uint64_t>& words,
        std::uint64_t position
    ) noexcept;
};

}  // namespace primerpair
