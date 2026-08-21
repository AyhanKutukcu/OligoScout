#include "primerpair/fm_index.hpp"
#include "primerpair/suffix_array_builder.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace primerpair {

FMIndex::FMIndex(
    std::string text,
    const std::size_t suffix_array_sample_rate
)
    : text_(std::move(text)),
      suffix_array_sample_rate_(
          suffix_array_sample_rate
      ) {

    if (text_.empty()) {
        throw std::invalid_argument(
            "Reference sequence cannot be empty."
        );
    }

    if (suffix_array_sample_rate_ == 0) {
        throw std::invalid_argument(
            "Suffix-array sample rate must be greater than zero."
        );
    }

    if (text_.find('$') != std::string::npos) {
        throw std::invalid_argument(
            "Reference must not contain the '$' sentinel."
        );
    }

    for (char& base : text_) {
        base = normalize_base(base);
    }

    /*
     * Build-time suffix-array konumları uint32_t
     * kullanır.
     *
     * İnsan kromozom shard'ları bu sınırın çok
     * altındadır. Sentinel için de bir pozisyon
     * ayırıyoruz.
     */
    if (
        text_.size() >=
        static_cast<std::size_t>(
            std::numeric_limits<std::uint32_t>::max()
        )
    ) {
        throw std::length_error(
            "Reference shard is too large for 32-bit suffix-array positions."
        );
    }

    /*
     * FM-index için benzersiz ve leksikografik olarak
     * en küçük sonlandırıcı.
     */
    text_.push_back('$');

    indexed_text_length_ =
        static_cast<std::uint64_t>(
            text_.size()
        );

    build_suffix_array();
    build_bwt();
    build_packed_bwt();
    build_c_table();
    build_rank_support();

    /*
     * Query yolu PackedBWT kullanır.
     * Checkpoint rank oluşturulduktan sonra klasik
     * byte-BWT artık kalıcı olarak gerekli değildir.
     */
    release_bwt();

    build_sampled_suffix_array();
    release_full_suffix_array();
    release_text();
}

char FMIndex::normalize_base(char base) {
    const char normalized =
        static_cast<char>(
            std::toupper(
                static_cast<unsigned char>(base)
            )
        );

    switch (normalized) {
        case '$':
        case 'A':
        case 'C':
        case 'G':
        case 'N':
        case 'T':
            return normalized;

        default:
            throw std::invalid_argument(
                std::string("Unsupported base: ") +
                normalized
            );
    }
}

std::size_t FMIndex::symbol_index(char symbol) {
    switch (normalize_base(symbol)) {
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
                "Unsupported FM-index symbol."
            );
    }
}

void FMIndex::build_suffix_array() {
    suffix_array_ =
        build_suffix_array_prefix_doubling(
            text_
        );
}

void FMIndex::build_bwt() {
    bwt_.clear();
    bwt_.reserve(text_.size());

    for (const std::uint64_t suffix_position :
         suffix_array_) {

        if (suffix_position == 0) {
            bwt_.push_back('$');
        } else {
            bwt_.push_back(
                text_.at(
                    static_cast<std::size_t>(
                        suffix_position - 1
                    )
                )
            );
        }
    }
}

void FMIndex::build_packed_bwt() {
    packed_bwt_.build(
        bwt_
    );

    if (
        packed_bwt_.size() !=
        static_cast<std::uint64_t>(
            bwt_.size()
        )
    ) {
        throw std::logic_error(
            "Packed BWT length does not match string BWT length."
        );
    }
}

void FMIndex::build_c_table() {
    std::array<
        std::uint64_t,
        kAlphabetSize
    > counts{};

    for (const char symbol : text_) {
        ++counts.at(
            symbol_index(symbol)
        );
    }

    std::uint64_t cumulative = 0;

    for (std::size_t index = 0;
         index < kAlphabetSize;
         ++index) {

        c_table_.at(index) = cumulative;
        cumulative += counts.at(index);
    }
}

void FMIndex::build_rank_support() {
    rank_support_.build(
        bwt_,
        CheckpointRank::kDefaultCheckpointRate
    );
}

void FMIndex::build_sampled_suffix_array() {
    sampled_sa_marker_words_.clear();
    sampled_sa_block_prefix_.clear();
    sampled_sa_values_.clear();

    const std::size_t row_count =
        suffix_array_.size();

    const std::size_t marker_word_count =
        (
            row_count + 63
        ) /
        64;

    sampled_sa_marker_words_.assign(
        marker_word_count,
        std::uint64_t{0}
    );

    const std::size_t estimated_samples =
        (
            row_count +
            suffix_array_sample_rate_ - 1
        ) /
        suffix_array_sample_rate_;

    sampled_sa_values_.reserve(
        estimated_samples
    );

    /*
     * Marker bitvector + sampled values oluştur.
     */
    for (std::size_t row = 0;
         row < row_count;
         ++row) {

        const std::uint32_t sa_value =
            suffix_array_.at(row);

        if (
            static_cast<std::size_t>(
                sa_value
            ) %
                suffix_array_sample_rate_
            !=
            0
        ) {
            continue;
        }

        const std::size_t word =
            row / 64;

        const std::size_t offset =
            row % 64;

        sampled_sa_marker_words_
            .at(word) |=
            (
                std::uint64_t{1}
                << offset
            );

        sampled_sa_values_.push_back(
            sa_value
        );
    }

    if (sampled_sa_values_.empty()) {
        throw std::logic_error(
            "Sampled suffix array is empty."
        );
    }

    /*
     * Her 8 marker word = 512 row için
     * bir prefix checkpoint.
     */
    const std::size_t block_count =
        (
            marker_word_count +
            kSampleMarkerWordsPerBlock - 1
        ) /
        kSampleMarkerWordsPerBlock;

    sampled_sa_block_prefix_.assign(
        block_count + 1,
        std::uint32_t{0}
    );

    std::uint32_t cumulative = 0;

    for (std::size_t block = 0;
         block < block_count;
         ++block) {

        sampled_sa_block_prefix_.at(
            block
        ) = cumulative;

        const std::size_t first_word =
            block *
            kSampleMarkerWordsPerBlock;

        const std::size_t last_word =
            std::min(
                marker_word_count,
                first_word +
                kSampleMarkerWordsPerBlock
            );

        for (std::size_t word = first_word;
             word < last_word;
             ++word) {

            cumulative +=
                static_cast<std::uint32_t>(
                    std::popcount(
                        sampled_sa_marker_words_
                            .at(word)
                    )
                );
        }
    }

    sampled_sa_block_prefix_.at(
        block_count
    ) = cumulative;

    if (
        static_cast<std::size_t>(
            cumulative
        ) !=
        sampled_sa_values_.size()
    ) {
        throw std::logic_error(
            "Sample-marker count does not match sampled SA values."
        );
    }
}

void FMIndex::release_full_suffix_array() {
    std::vector<std::uint32_t> empty;
    suffix_array_.swap(empty);
}

void FMIndex::release_bwt() {
    std::string empty;
    bwt_.swap(empty);
}

std::string FMIndex::bwt_string() const {
    std::string decoded;

    decoded.reserve(
        static_cast<std::size_t>(
            packed_bwt_.size()
        )
    );

    for (
        std::uint64_t position = 0;
        position < packed_bwt_.size();
        ++position
    ) {
        decoded.push_back(
            packed_bwt_.at(position)
        );
    }

    return decoded;
}

void FMIndex::release_text() {
    std::string empty;
    text_.swap(empty);
}

bool FMIndex::find_sampled_sa(
    const std::uint64_t row,
    std::uint64_t& value
) const {
    if (row >= indexed_text_length_) {
        return false;
    }

    const std::size_t word =
        static_cast<std::size_t>(
            row / 64
        );

    const std::size_t offset =
        static_cast<std::size_t>(
            row % 64
        );

    if (
        word >=
        sampled_sa_marker_words_.size()
    ) {
        return false;
    }

    const std::uint64_t marker_word =
        sampled_sa_marker_words_.at(
            word
        );

    const std::uint64_t marker_bit =
        std::uint64_t{1}
        << offset;

    /*
     * İlgili FM-index satırı sample değilse
     * hiçbir ek rank hesabı yapmadan çık.
     */
    if (
        (marker_word & marker_bit) == 0
    ) {
        return false;
    }

    const std::size_t block =
        word /
        kSampleMarkerWordsPerBlock;

    std::uint32_t sample_index =
        sampled_sa_block_prefix_.at(
            block
        );

    const std::size_t first_word =
        block *
        kSampleMarkerWordsPerBlock;

    /*
     * Blok başlangıcından mevcut word'e kadar
     * en fazla 7 tam 64-bit word sayılır.
     */
    for (std::size_t current = first_word;
         current < word;
         ++current) {

        sample_index +=
            static_cast<std::uint32_t>(
                std::popcount(
                    sampled_sa_marker_words_
                        .at(current)
                )
            );
    }

    /*
     * Mevcut word içinde row'dan önceki marker'lar.
     */
    if (offset != 0) {
        const std::uint64_t before_mask =
            (
                std::uint64_t{1}
                << offset
            ) -
            1;

        sample_index +=
            static_cast<std::uint32_t>(
                std::popcount(
                    marker_word &
                    before_mask
                )
            );
    }

    if (
        static_cast<std::size_t>(
            sample_index
        ) >=
        sampled_sa_values_.size()
    ) {
        throw std::logic_error(
            "Sample-marker index exceeds sampled SA values."
        );
    }

    value =
        static_cast<std::uint64_t>(
            sampled_sa_values_.at(
                static_cast<std::size_t>(
                    sample_index
                )
            )
        );

    return true;
}

std::uint64_t FMIndex::lf(
    const std::uint64_t row
) const {
    if (row >= packed_bwt_.size()) {
        throw std::out_of_range(
            "LF row exceeds BWT length."
        );
    }

    /*
     * Numeric LF hot-path.
     *
     * Alfabe:
     *   0=$, 1=A, 2=C, 3=G, 4=N, 5=T
     *
     * Eski char decode + iki ayrı symbol_index()
     * zincirini kaldırır.
     */
    const std::size_t symbol_id =
        packed_bwt_.symbol_id_at(
            row
        );

    return
        c_table_.at(
            symbol_id
        ) +
        rank_support_.rank_by_symbol_id(
            packed_bwt_,
            symbol_id,
            row
        );
}

std::uint64_t FMIndex::resolve_sa(
    const std::uint64_t row
) const {
    if (row >= packed_bwt_.size()) {
        throw std::out_of_range(
            "Suffix-array row exceeds BWT length."
        );
    }

    std::uint64_t current_row = row;
    std::uint64_t sampled_value = 0;

    for (std::size_t steps = 0;
         steps < suffix_array_sample_rate_;
         ++steps) {

        if (
            find_sampled_sa(
                current_row,
                sampled_value
            )
        ) {
            return (
                sampled_value +
                static_cast<std::uint64_t>(
                    steps
                )
            ) %
            indexed_text_length_;
        }

        current_row = lf(current_row);
    }

    throw std::logic_error(
        "Sampled suffix-array lookup exceeded sample rate."
    );
}

Interval FMIndex::backward_extend(
    const Interval& interval,
    char base
) const {
    if (interval.begin > interval.end) {
        throw std::invalid_argument(
            "FM-index interval begin exceeds interval end."
        );
    }

    if (interval.end > packed_bwt_.size()) {
        throw std::out_of_range(
            "FM-index interval exceeds BWT length."
        );
    }

    const std::size_t index =
        symbol_index(base);

    const std::uint64_t new_begin =
        c_table_.at(index) +
        rank_support_.rank(
            packed_bwt_,
            base,
            interval.begin
        );

    const std::uint64_t new_end =
        c_table_.at(index) +
        rank_support_.rank(
            packed_bwt_,
            base,
            interval.end
        );

    return Interval{
        .begin = new_begin,
        .end = new_end
    };
}

std::array<Interval, 6>
FMIndex::backward_extend_all(
    const Interval& interval
) const {
    if (interval.begin > interval.end) {
        throw std::invalid_argument(
            "FM-index interval begin exceeds interval end."
        );
    }

    if (
        interval.end >
        packed_bwt_.size()
    ) {
        throw std::out_of_range(
            "FM-index interval exceeds BWT length."
        );
    }

    /*
     * Parent intervalin iki ucunda tüm alfabetik
     * rank değerlerini yalnızca birer kez hesapla.
     *
     * Alphabet:
     *   0=$
     *   1=A
     *   2=C
     *   3=G
     *   4=N
     *   5=T
     */
    const auto before =
        rank_support_.rank_all(
            packed_bwt_,
            interval.begin
        );

    const auto after =
        rank_support_.rank_all(
            packed_bwt_,
            interval.end
        );

    std::array<Interval, 6>
        children{};

    for (
        std::size_t i = 0;
        i < children.size();
        ++i
    ) {
        children.at(i) = Interval{
            .begin =
                c_table_.at(i) +
                before.at(i),

            .end =
                c_table_.at(i) +
                after.at(i)
        };
    }

    return children;
}


std::array<std::uint64_t, 6>
FMIndex::interval_symbol_counts(
    const Interval& interval
) const {
    if (interval.begin > interval.end) {
        throw std::invalid_argument(
            "FM-index interval begin exceeds end."
        );
    }

    if (
        interval.end >
        packed_bwt_.size()
    ) {
        throw std::out_of_range(
            "FM-index interval exceeds BWT length."
        );
    }

    const auto before =
        rank_support_.rank_all(
            packed_bwt_,
            interval.begin
        );

    const auto after =
        rank_support_.rank_all(
            packed_bwt_,
            interval.end
        );

    std::array<std::uint64_t, 6>
        counts{};

    for (
        std::size_t i = 0;
        i < counts.size();
        ++i
    ) {
        counts[i] =
            after[i] -
            before[i];
    }

    return counts;
}


Interval FMIndex::backward_search(
    std::string_view pattern
) const {
    Interval interval{
        .begin = 0,
        .end = static_cast<std::uint64_t>(
            packed_bwt_.size()
        )
    };

    if (pattern.empty()) {
        return interval;
    }

    for (auto iterator = pattern.rbegin();
         iterator != pattern.rend();
         ++iterator) {

        interval = backward_extend(
            interval,
            *iterator
        );

        if (interval.empty()) {
            break;
        }
    }

    return interval;
}

std::vector<std::uint64_t>
FMIndex::locate_unsorted(
    const Interval& interval
) const {
    if (interval.begin > interval.end) {
        throw std::invalid_argument(
            "FM-index locate interval begin exceeds end."
        );
    }

    if (
        interval.end >
        indexed_text_length_
    ) {
        throw std::out_of_range(
            "FM-index locate interval exceeds row count."
        );
    }

    std::vector<std::uint64_t>
        positions;

    positions.reserve(
        static_cast<std::size_t>(
            interval.size()
        )
    );

    const std::uint64_t sentinel_position =
        indexed_text_length_ - 1;

    for (
        std::uint64_t row = interval.begin;
        row < interval.end;
        ++row
    ) {
        const std::uint64_t position =
            resolve_sa(
                row
            );

        if (
            position !=
            sentinel_position
        ) {
            positions.push_back(
                position
            );
        }
    }

    return positions;
}


std::vector<std::uint64_t>
FMIndex::locate(
    const Interval& interval
) const {
    auto positions =
        locate_unsorted(
            interval
        );

    std::sort(
        positions.begin(),
        positions.end()
    );

    return positions;
}


}  // namespace primerpair
