#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/packed_bwt.hpp"
#include "primerpair/rank_support.hpp"

namespace primerpair {

struct Interval {
    std::uint64_t begin{0};
    std::uint64_t end{0};

    [[nodiscard]] bool empty() const noexcept {
        return begin >= end;
    }

    [[nodiscard]] std::uint64_t size() const noexcept {
        return end > begin ? end - begin : 0;
    }
};

class FMIndex {
    friend class PpfmIO;

public:
    static constexpr std::size_t
        kDefaultSuffixArraySampleRate = 8;

    explicit FMIndex(
        std::string text,
        std::size_t suffix_array_sample_rate =
            kDefaultSuffixArraySampleRate
    );

    [[nodiscard]] Interval backward_extend(
        const Interval& interval,
        char base
    ) const;


    /*
     * Aynı parent FM intervalinden bütün alfabetik
     * child interval'larını tek rank_all(begin/end)
     * çiftiyle üretir.
     *
     * Sıra:
     *   0=$, 1=A, 2=C, 3=G, 4=N, 5=T
     */
    [[nodiscard]]
    std::array<Interval, 6>
    backward_extend_all(
        const Interval& interval
    ) const;

    [[nodiscard]] Interval backward_search(
        std::string_view pattern
    ) const;

    /*
     * BWT intervali [begin,end) içindeki
     * $, A, C, G, N, T sembol sayılarını döndürür.
     *
     * Lexicographic sıra:
     * 0=$, 1=A, 2=C, 3=G, 4=N, 5=T
     */
    [[nodiscard]]
    std::array<std::uint64_t, 6>
    interval_symbol_counts(
        const Interval& interval
    ) const;

    [[nodiscard]] std::vector<std::uint64_t> locate(
        const Interval& interval
    ) const;


    /*
     * locate() ile aynı koordinatları döndürür ancak
     * sonuçları sıralamaz.
     *
     * Hot-path kullanımı içindir; çağıran taraf sonuç
     * sırasına ihtiyaç duymuyorsa branch başına sort
     * maliyetini ortadan kaldırır.
     */
    [[nodiscard]]
    std::vector<std::uint64_t> locate_unsorted(
        const Interval& interval
    ) const;

    [[nodiscard]] std::uint64_t lf(
        std::uint64_t row
    ) const;

    [[nodiscard]] std::uint64_t indexed_text_length() const noexcept {
        return indexed_text_length_;
    }

    [[nodiscard]] std::size_t text_memory_bytes() const noexcept {
        return text_.size();
    }

    /*
     * Debug/test amacıyla PackedBWT'den BWT stringini
     * yeniden oluşturur. Hot query path'te kullanılmaz.
     */
    [[nodiscard]]
    std::string bwt_string() const;

    [[nodiscard]]
    std::uint64_t bwt_size() const noexcept {
        return packed_bwt_.size();
    }

    /*
     * Build sonrasında tutulmaya devam eden klasik
     * BWT string belleği. release_bwt() sonrası 0 olmalı.
     */
    [[nodiscard]]
    std::size_t bwt_memory_bytes() const noexcept {
        return bwt_.size();
    }

    [[nodiscard]]
    std::size_t packed_bwt_memory_bytes() const noexcept {
        return packed_bwt_.memory_bytes();
    }

    [[nodiscard]] std::size_t rank_checkpoint_count() const noexcept {
        return rank_support_.checkpoint_count();
    }

    [[nodiscard]] std::size_t rank_memory_bytes() const noexcept {
        return rank_support_.memory_bytes();
    }

    [[nodiscard]] std::size_t suffix_array_sample_rate() const noexcept {
        return suffix_array_sample_rate_;
    }

    [[nodiscard]] std::size_t sampled_sa_count() const noexcept {
        return sampled_sa_values_.size();
    }

    [[nodiscard]]
    std::size_t sampled_sa_marker_memory_bytes() const noexcept {
        return
            sampled_sa_marker_words_.size() *
            sizeof(std::uint64_t);
    }

    [[nodiscard]]
    std::size_t sampled_sa_prefix_memory_bytes() const noexcept {
        return
            sampled_sa_block_prefix_.size() *
            sizeof(std::uint32_t);
    }

    [[nodiscard]]
    std::size_t sampled_sa_values_memory_bytes() const noexcept {
        return
            sampled_sa_values_.size() *
            sizeof(std::uint32_t);
    }

    [[nodiscard]]
    std::size_t sampled_sa_memory_bytes() const noexcept {
        return
            sampled_sa_marker_memory_bytes() +
            sampled_sa_prefix_memory_bytes() +
            sampled_sa_values_memory_bytes();
    }

private:
    static constexpr std::size_t kAlphabetSize = 6;

    /*
     * Sample-marker rank için bir prefix checkpoint
     * her 8 x 64 = 512 FM-index satırında tutulur.
     */
    static constexpr std::size_t
        kSampleMarkerWordsPerBlock = 8;

    /*
     * Yalnızca indeks oluşturma sırasında tutulur.
     * Constructor sonunda release_text() ile bırakılır.
     */
    std::string text_;

    /*
     * Sentinel dahil FM-index satır sayısı.
     * Query sırasında text_.size() yerine kullanılır.
     */
    std::uint64_t indexed_text_length_{0};

    /*
     * Yalnızca indeks oluşturma sırasında kullanılır.
     * Sampled SA oluşturulduktan sonra bellekten silinir.
     */
    std::vector<std::uint32_t> suffix_array_;

    /*
     * Build sırasında klasik BWT stringi oluşturulur.
     * Query yolu PackedBWT üzerinden çalışacaktır.
     *
     * Bu aşamada bwt_ henüz release edilmiyor;
     * önce packed-query performansını izole ediyoruz.
     */
    std::string bwt_;

    PackedBWT packed_bwt_;

    std::array<
        std::uint64_t,
        kAlphabetSize
    > c_table_{};

    CheckpointRank rank_support_;

    std::size_t suffix_array_sample_rate_;

    /*
     * sampled_sa_marker_words_:
     *   Bit 1 ise ilgili FM-index satırında sampled SA vardır.
     *
     * sampled_sa_block_prefix_[b]:
     *   b numaralı 512-row bloktan önce bulunan toplam
     *   sampled-SA sayısıdır.
     *
     * sampled_sa_values_:
     *   Marker bitleri satır sırasıyla tarandığında karşılık
     *   gelen gerçek SA değerlerini tutar.
     */
    std::vector<std::uint64_t>
        sampled_sa_marker_words_;

    std::vector<std::uint32_t>
        sampled_sa_block_prefix_;

    std::vector<std::uint32_t>
        sampled_sa_values_;

    [[nodiscard]] static char normalize_base(
        char base
    );

    [[nodiscard]] static std::size_t symbol_index(
        char symbol
    );

    void build_suffix_array();
    void build_bwt();
    void build_packed_bwt();
    void build_c_table();
    void build_rank_support();
    void build_sampled_suffix_array();
    void release_full_suffix_array();
    void release_bwt();
    void release_text();

    [[nodiscard]] bool find_sampled_sa(
        std::uint64_t row,
        std::uint64_t& value
    ) const;

    [[nodiscard]] std::uint64_t resolve_sa(
        std::uint64_t row
    ) const;
};

}  // namespace primerpair
