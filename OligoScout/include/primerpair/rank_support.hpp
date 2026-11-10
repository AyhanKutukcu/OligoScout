#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace primerpair {

class PackedBWT;

class CheckpointRank {
    friend class PpfmIO;

public:
    static constexpr std::size_t kAlphabetSize = 6;
    static constexpr std::size_t kDefaultCheckpointRate = 128;

    void build(
        std::string_view bwt,
        std::size_t checkpoint_rate =
            kDefaultCheckpointRate
    );

    [[nodiscard]]
    std::uint64_t rank(
        std::string_view bwt,
        char symbol,
        std::uint64_t position
    ) const;

    /*
     * Packed-BWT query yolu.
     *
     * Checkpoint'a kadar olan kümülatif sayaç aynı,
     * checkpoint sonrası kısa aralık PackedBWT::count()
     * ile popcount kullanılarak hesaplanır.
     */
    [[nodiscard]]
    std::uint64_t rank(
        const PackedBWT& bwt,
        char symbol,
        std::uint64_t position
    ) const;


    /*
     * Numeric alphabet-ID hot-path.
     *
     * 0=$, 1=A, 2=C, 3=G, 4=N, 5=T
     */
    [[nodiscard]]
    std::uint64_t rank_by_symbol_id(
        const PackedBWT& bwt,
        std::size_t symbol_id,
        std::uint64_t position
    ) const;


    /*
     * $, A, C, G, N, T rank değerlerini
     * tek checkpoint + tek PackedBWT geçişiyle döndürür.
     */
    [[nodiscard]]
    std::array<std::uint64_t, kAlphabetSize>
    rank_all(
        const PackedBWT& bwt,
        std::uint64_t position
    ) const;

    [[nodiscard]]
    std::size_t checkpoint_rate() const noexcept {
        return checkpoint_rate_;
    }

    [[nodiscard]]
    std::size_t checkpoint_count() const noexcept {
        return checkpoints_.size();
    }

    [[nodiscard]]
    std::size_t memory_bytes() const noexcept {
        return checkpoints_.size() *
               sizeof(
                   std::array<
                       std::uint32_t,
                       kAlphabetSize
                   >
               );
    }

private:
    std::size_t checkpoint_rate_{
        kDefaultCheckpointRate
    };

    std::vector<
        std::array<
            std::uint32_t,
            kAlphabetSize
        >
    > checkpoints_;

    [[nodiscard]]
    static std::size_t symbol_index(char symbol);
};

}  // namespace primerpair
