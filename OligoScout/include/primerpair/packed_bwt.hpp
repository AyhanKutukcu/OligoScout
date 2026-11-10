#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace primerpair {

class PackedBWT {
    friend class PpfmIO;

public:
    PackedBWT() = default;

    explicit PackedBWT(
        std::string_view bwt
    );

    void build(
        std::string_view bwt
    );

    [[nodiscard]]
    char at(
        std::uint64_t position
    ) const;


    /*
     * BWT sembolünü doğrudan FM-index alfabetik
     * indeksine decode eder.
     *
     * 0=$, 1=A, 2=C, 3=G, 4=N, 5=T
     *
     * LF hot-path'te char dönüşümünü ortadan
     * kaldırmak için kullanılır.
     */
    [[nodiscard]]
    std::size_t symbol_id_at(
        std::uint64_t position
    ) const;

    /*
     * [begin, end) aralığında symbol sayısını döndürür.
     *
     * Bu fonksiyon ileride checkpoint rank tarafından
     * en fazla ~128 bp'lik aralıkları popcount ile
     * saymak için kullanılacaktır.
     */
    [[nodiscard]]
    std::uint64_t count(
        char symbol,
        std::uint64_t begin,
        std::uint64_t end
    ) const;


    /*
     * Numeric alphabet-ID variant of count().
     *
     * 0=$, 1=A, 2=C, 3=G, 4=N, 5=T
     *
     * char normalization/switch maliyetini LF
     * hot-path'ten çıkarır.
     */
    [[nodiscard]]
    std::uint64_t count_by_symbol_id(
        std::size_t symbol_id,
        std::uint64_t begin,
        std::uint64_t end
    ) const;


    /*
     * [begin,end) aralığındaki tüm BWT sembollerini
     * tek bit-plane geçişinde sayar.
     *
     * Sıra:
     *   0=$, 1=A, 2=C, 3=G, 4=N, 5=T
     */
    [[nodiscard]]
    std::array<std::uint64_t, 6>
    count_all(
        std::uint64_t begin,
        std::uint64_t end
    ) const;

    [[nodiscard]]
    std::uint64_t size() const noexcept {
        return length_;
    }

    [[nodiscard]]
    std::uint32_t sentinel_position() const noexcept {
        return sentinel_position_;
    }

    [[nodiscard]]
    std::size_t sequence_memory_bytes() const noexcept {
        return
            (
                low_bits_.size() +
                high_bits_.size()
            ) *
            sizeof(std::uint64_t);
    }

    [[nodiscard]]
    std::size_t n_mask_memory_bytes() const noexcept {
        return
            n_mask_words_.size() *
            sizeof(std::uint64_t);
    }

    [[nodiscard]]
    std::size_t memory_bytes() const noexcept {
        return
            sequence_memory_bytes() +
            n_mask_memory_bytes() +
            sizeof(sentinel_position_) +
            sizeof(length_);
    }

private:
    /*
     * 64 BWT sembolü / uint64_t.
     *
     * A = 00
     * C = 01
     * G = 10
     * T = 11
     *
     * low_bits_  -> düşük bit
     * high_bits_ -> yüksek bit
     */
    std::vector<std::uint64_t> low_bits_;
    std::vector<std::uint64_t> high_bits_;

    /*
     * Bit = 1 ise pozisyon N'dir.
     */
    std::vector<std::uint64_t> n_mask_words_;

    std::uint64_t length_{0};

    std::uint32_t sentinel_position_{0};

    [[nodiscard]]
    static std::uint8_t encode_acgt(
        char symbol
    );

    [[nodiscard]]
    std::uint8_t acgt_code_at(
        std::uint64_t position
    ) const;

    [[nodiscard]]
    bool is_n(
        std::uint64_t position
    ) const;
};

}  // namespace primerpair
