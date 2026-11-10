#include "primerpair/ppfm_io.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace primerpair {

namespace {

constexpr std::array<char, 8>
    kPpfmMagic{
        'P', 'P', 'F', 'M',
        'S', 'H', 'R', 'D'
    };

constexpr std::uint64_t
    kMaximumChromosomeNameLength = 1024;


void require_write(
    std::ostream& output
) {
    if (!output) {
        throw std::runtime_error(
            "PPFM binary write failed."
        );
    }
}


void require_read(
    std::istream& input
) {
    if (!input) {
        throw std::runtime_error(
            "PPFM binary read failed or file is truncated."
        );
    }
}


void write_u32(
    std::ostream& output,
    const std::uint32_t value
) {
    const std::array<unsigned char, 4> bytes{
        static_cast<unsigned char>(
            value & 0xffU
        ),
        static_cast<unsigned char>(
            (value >> 8U) & 0xffU
        ),
        static_cast<unsigned char>(
            (value >> 16U) & 0xffU
        ),
        static_cast<unsigned char>(
            (value >> 24U) & 0xffU
        )
    };

    output.write(
        reinterpret_cast<const char*>(
            bytes.data()
        ),
        static_cast<std::streamsize>(
            bytes.size()
        )
    );

    require_write(output);
}


void write_u64(
    std::ostream& output,
    const std::uint64_t value
) {
    std::array<unsigned char, 8> bytes{};

    for (std::size_t i = 0;
         i < bytes.size();
         ++i) {

        bytes.at(i) =
            static_cast<unsigned char>(
                (
                    value >>
                    static_cast<unsigned>(
                        i * 8
                    )
                ) &
                0xffULL
            );
    }

    output.write(
        reinterpret_cast<const char*>(
            bytes.data()
        ),
        static_cast<std::streamsize>(
            bytes.size()
        )
    );

    require_write(output);
}


std::uint32_t read_u32(
    std::istream& input
) {
    std::array<unsigned char, 4> bytes{};

    input.read(
        reinterpret_cast<char*>(
            bytes.data()
        ),
        static_cast<std::streamsize>(
            bytes.size()
        )
    );

    require_read(input);

    return
        static_cast<std::uint32_t>(
            bytes.at(0)
        )
        |
        (
            static_cast<std::uint32_t>(
                bytes.at(1)
            ) << 8U
        )
        |
        (
            static_cast<std::uint32_t>(
                bytes.at(2)
            ) << 16U
        )
        |
        (
            static_cast<std::uint32_t>(
                bytes.at(3)
            ) << 24U
        );
}


std::uint64_t read_u64(
    std::istream& input
) {
    std::array<unsigned char, 8> bytes{};

    input.read(
        reinterpret_cast<char*>(
            bytes.data()
        ),
        static_cast<std::streamsize>(
            bytes.size()
        )
    );

    require_read(input);

    std::uint64_t value = 0;

    for (std::size_t i = 0;
         i < bytes.size();
         ++i) {

        value |=
            (
                static_cast<std::uint64_t>(
                    bytes.at(i)
                )
                <<
                static_cast<unsigned>(
                    i * 8
                )
            );
    }

    return value;
}


std::size_t checked_size_t(
    const std::uint64_t value,
    const char* label
) {
    if (
        value >
        static_cast<std::uint64_t>(
            std::numeric_limits<
                std::size_t
            >::max()
        )
    ) {
        throw std::runtime_error(
            std::string(
                "PPFM "
            ) +
            label +
            " exceeds platform size_t."
        );
    }

    return static_cast<std::size_t>(
        value
    );
}


void write_u64_vector(
    std::ostream& output,
    const std::vector<std::uint64_t>& values
) {
    write_u64(
        output,
        static_cast<std::uint64_t>(
            values.size()
        )
    );

    if (values.empty()) {
        return;
    }

    if constexpr (
        std::endian::native ==
        std::endian::little
    ) {
        output.write(
            reinterpret_cast<const char*>(
                values.data()
            ),
            static_cast<std::streamsize>(
                values.size() *
                sizeof(std::uint64_t)
            )
        );

        require_write(output);

    } else {

        for (const auto value : values) {
            write_u64(
                output,
                value
            );
        }
    }
}


std::vector<std::uint64_t>
read_u64_vector(
    std::istream& input
) {
    const std::uint64_t raw_count =
        read_u64(input);

    const std::size_t count =
        checked_size_t(
            raw_count,
            "uint64 vector length"
        );

    if (
        count >
        std::numeric_limits<
            std::size_t
        >::max() /
        sizeof(std::uint64_t)
    ) {
        throw std::runtime_error(
            "PPFM uint64 vector byte size overflow."
        );
    }

    std::vector<std::uint64_t>
        values(
            count
        );

    if (values.empty()) {
        return values;
    }

    if constexpr (
        std::endian::native ==
        std::endian::little
    ) {
        input.read(
            reinterpret_cast<char*>(
                values.data()
            ),
            static_cast<std::streamsize>(
                values.size() *
                sizeof(std::uint64_t)
            )
        );

        require_read(input);

    } else {

        for (auto& value : values) {
            value =
                read_u64(input);
        }
    }

    return values;
}


void write_u32_vector(
    std::ostream& output,
    const std::vector<std::uint32_t>& values
) {
    write_u64(
        output,
        static_cast<std::uint64_t>(
            values.size()
        )
    );

    if (values.empty()) {
        return;
    }

    if constexpr (
        std::endian::native ==
        std::endian::little
    ) {
        output.write(
            reinterpret_cast<const char*>(
                values.data()
            ),
            static_cast<std::streamsize>(
                values.size() *
                sizeof(std::uint32_t)
            )
        );

        require_write(output);

    } else {

        for (const auto value : values) {
            write_u32(
                output,
                value
            );
        }
    }
}


std::vector<std::uint32_t>
read_u32_vector(
    std::istream& input
) {
    const std::uint64_t raw_count =
        read_u64(input);

    const std::size_t count =
        checked_size_t(
            raw_count,
            "uint32 vector length"
        );

    if (
        count >
        std::numeric_limits<
            std::size_t
        >::max() /
        sizeof(std::uint32_t)
    ) {
        throw std::runtime_error(
            "PPFM uint32 vector byte size overflow."
        );
    }

    std::vector<std::uint32_t>
        values(
            count
        );

    if (values.empty()) {
        return values;
    }

    if constexpr (
        std::endian::native ==
        std::endian::little
    ) {
        input.read(
            reinterpret_cast<char*>(
                values.data()
            ),
            static_cast<std::streamsize>(
                values.size() *
                sizeof(std::uint32_t)
            )
        );

        require_read(input);

    } else {

        for (auto& value : values) {
            value =
                read_u32(input);
        }
    }

    return values;
}


void write_string(
    std::ostream& output,
    const std::string& value
) {
    if (
        value.size() >
        std::numeric_limits<
            std::uint32_t
        >::max()
    ) {
        throw std::length_error(
            "PPFM string is too large."
        );
    }

    write_u32(
        output,
        static_cast<std::uint32_t>(
            value.size()
        )
    );

    if (!value.empty()) {
        output.write(
            value.data(),
            static_cast<std::streamsize>(
                value.size()
            )
        );

        require_write(output);
    }
}


std::string read_string(
    std::istream& input
) {
    const std::uint32_t length =
        read_u32(input);

    if (
        static_cast<std::uint64_t>(
            length
        ) >
        kMaximumChromosomeNameLength
    ) {
        throw std::runtime_error(
            "PPFM chromosome name is unreasonably long."
        );
    }

    std::string value(
        static_cast<std::size_t>(
            length
        ),
        '\0'
    );

    if (!value.empty()) {
        input.read(
            value.data(),
            static_cast<std::streamsize>(
                value.size()
            )
        );

        require_read(input);
    }

    return value;
}


std::uint64_t word_count_for_length(
    const std::uint64_t length
) {
    return (
        length + 63ULL
    ) / 64ULL;
}

}  // namespace


void PpfmIO::write_reference(
    std::ostream& output,
    const PackedReference& reference
) {
    write_u64(
        output,
        reference.length_
    );

    write_u64_vector(
        output,
        reference.low_bits_
    );

    write_u64_vector(
        output,
        reference.high_bits_
    );

    write_u64_vector(
        output,
        reference.n_mask_
    );
}


void PpfmIO::read_reference(
    std::istream& input,
    PackedReference& reference
) {
    reference.length_ =
        read_u64(input);

    reference.low_bits_ =
        read_u64_vector(input);

    reference.high_bits_ =
        read_u64_vector(input);

    reference.n_mask_ =
        read_u64_vector(input);
}


void PpfmIO::write_packed_bwt(
    std::ostream& output,
    const PackedBWT& bwt
) {
    write_u64(
        output,
        bwt.length_
    );

    write_u32(
        output,
        bwt.sentinel_position_
    );

    write_u64_vector(
        output,
        bwt.low_bits_
    );

    write_u64_vector(
        output,
        bwt.high_bits_
    );

    write_u64_vector(
        output,
        bwt.n_mask_words_
    );
}


void PpfmIO::read_packed_bwt(
    std::istream& input,
    PackedBWT& bwt
) {
    bwt.length_ =
        read_u64(input);

    bwt.sentinel_position_ =
        read_u32(input);

    bwt.low_bits_ =
        read_u64_vector(input);

    bwt.high_bits_ =
        read_u64_vector(input);

    bwt.n_mask_words_ =
        read_u64_vector(input);
}


void PpfmIO::write_rank(
    std::ostream& output,
    const CheckpointRank& rank
) {
    write_u64(
        output,
        static_cast<std::uint64_t>(
            rank.checkpoint_rate_
        )
    );

    write_u64(
        output,
        static_cast<std::uint64_t>(
            rank.checkpoints_.size()
        )
    );

    for (const auto& checkpoint :
         rank.checkpoints_) {

        for (const auto value :
             checkpoint) {

            write_u32(
                output,
                value
            );
        }
    }
}


void PpfmIO::read_rank(
    std::istream& input,
    CheckpointRank& rank
) {
    rank.checkpoint_rate_ =
        checked_size_t(
            read_u64(input),
            "checkpoint rate"
        );

    const std::size_t count =
        checked_size_t(
            read_u64(input),
            "checkpoint count"
        );

    rank.checkpoints_.assign(
        count,
        {}
    );

    for (auto& checkpoint :
         rank.checkpoints_) {

        for (auto& value :
             checkpoint) {

            value =
                read_u32(input);
        }
    }
}


void PpfmIO::write_fm_index(
    std::ostream& output,
    const FMIndex& index
) {
    write_u64(
        output,
        index.indexed_text_length_
    );

    write_u64(
        output,
        static_cast<std::uint64_t>(
            index.suffix_array_sample_rate_
        )
    );

    write_packed_bwt(
        output,
        index.packed_bwt_
    );

    for (const auto value :
         index.c_table_) {

        write_u64(
            output,
            value
        );
    }

    write_rank(
        output,
        index.rank_support_
    );

    write_u64_vector(
        output,
        index.sampled_sa_marker_words_
    );

    write_u32_vector(
        output,
        index.sampled_sa_block_prefix_
    );

    write_u32_vector(
        output,
        index.sampled_sa_values_
    );
}


void PpfmIO::read_fm_index(
    std::istream& input,
    FMIndex& index
) {
    index.indexed_text_length_ =
        read_u64(input);

    index.suffix_array_sample_rate_ =
        checked_size_t(
            read_u64(input),
            "suffix-array sample rate"
        );

    read_packed_bwt(
        input,
        index.packed_bwt_
    );

    for (auto& value :
         index.c_table_) {

        value =
            read_u64(input);
    }

    read_rank(
        input,
        index.rank_support_
    );

    index.sampled_sa_marker_words_ =
        read_u64_vector(input);

    index.sampled_sa_block_prefix_ =
        read_u32_vector(input);

    index.sampled_sa_values_ =
        read_u32_vector(input);

    /*
     * Build-only state must remain absent after
     * persistent load.
     */
    index.text_.clear();
    index.text_.shrink_to_fit();

    index.suffix_array_.clear();
    index.suffix_array_.shrink_to_fit();

    index.bwt_.clear();
    index.bwt_.shrink_to_fit();
}


void PpfmIO::validate_reference(
    const PackedReference& reference
) {
    if (reference.length_ == 0) {
        throw std::runtime_error(
            "PPFM reference length is zero."
        );
    }

    const std::uint64_t expected_words =
        word_count_for_length(
            reference.length_
        );

    if (
        reference.low_bits_.size() !=
            expected_words
        ||
        reference.high_bits_.size() !=
            expected_words
        ||
        reference.n_mask_.size() !=
            expected_words
    ) {
        throw std::runtime_error(
            "PPFM packed-reference word count mismatch."
        );
    }
}


void PpfmIO::validate_fm_index(
    const FMIndex& index,
    const std::uint64_t reference_length
) {
    if (
        reference_length ==
        std::numeric_limits<
            std::uint64_t
        >::max()
    ) {
        throw std::runtime_error(
            "PPFM reference length cannot accommodate sentinel."
        );
    }

    const std::uint64_t expected_length =
        reference_length + 1ULL;

    if (
        index.indexed_text_length_ !=
        expected_length
    ) {
        throw std::runtime_error(
            "PPFM FM-index text length does not match reference + sentinel."
        );
    }

    if (
        index.packed_bwt_.length_ !=
        expected_length
    ) {
        throw std::runtime_error(
            "PPFM PackedBWT length mismatch."
        );
    }

    if (
        static_cast<std::uint64_t>(
            index.packed_bwt_
                .sentinel_position_
        ) >=
        expected_length
    ) {
        throw std::runtime_error(
            "PPFM sentinel position is out of range."
        );
    }

    const std::uint64_t expected_words =
        word_count_for_length(
            expected_length
        );

    if (
        index.packed_bwt_
            .low_bits_.size() !=
            expected_words
        ||
        index.packed_bwt_
            .high_bits_.size() !=
            expected_words
        ||
        index.packed_bwt_
            .n_mask_words_.size() !=
            expected_words
    ) {
        throw std::runtime_error(
            "PPFM PackedBWT word count mismatch."
        );
    }

    if (
        index.suffix_array_sample_rate_ ==
        0
    ) {
        throw std::runtime_error(
            "PPFM suffix-array sample rate is zero."
        );
    }

    if (
        index.rank_support_
            .checkpoint_rate_ ==
        0
    ) {
        throw std::runtime_error(
            "PPFM checkpoint rate is zero."
        );
    }

    if (
        index.rank_support_
            .checkpoints_.empty()
    ) {
        throw std::runtime_error(
            "PPFM rank checkpoint table is empty."
        );
    }

    if (
        index.c_table_.front() != 0
    ) {
        throw std::runtime_error(
            "PPFM C-table must begin at zero."
        );
    }

    for (std::size_t i = 1;
         i < index.c_table_.size();
         ++i) {

        if (
            index.c_table_.at(i) <
            index.c_table_.at(i - 1)
        ) {
            throw std::runtime_error(
                "PPFM C-table is not monotonic."
            );
        }
    }

    const std::size_t expected_marker_words =
        static_cast<std::size_t>(
            (
                expected_length +
                63ULL
            ) /
            64ULL
        );

    if (
        index.sampled_sa_marker_words_
            .size() !=
        expected_marker_words
    ) {
        throw std::runtime_error(
            "PPFM sampled-SA marker length mismatch."
        );
    }

    const std::size_t block_count =
        (
            expected_marker_words +
            FMIndex::
                kSampleMarkerWordsPerBlock -
            1
        ) /
        FMIndex::
            kSampleMarkerWordsPerBlock;

    if (
        index.sampled_sa_block_prefix_
            .size() !=
        block_count + 1
    ) {
        throw std::runtime_error(
            "PPFM sampled-SA block-prefix length mismatch."
        );
    }

    std::uint64_t marker_total = 0;

    for (const auto word :
         index.sampled_sa_marker_words_) {

        marker_total +=
            static_cast<std::uint64_t>(
                std::popcount(word)
            );
    }

    if (
        marker_total !=
        static_cast<std::uint64_t>(
            index.sampled_sa_values_
                .size()
        )
    ) {
        throw std::runtime_error(
            "PPFM sampled-SA marker/value count mismatch."
        );
    }

    if (
        index.sampled_sa_block_prefix_
            .empty()
        ||
        static_cast<std::uint64_t>(
            index.sampled_sa_block_prefix_
                .back()
        ) !=
            marker_total
    ) {
        throw std::runtime_error(
            "PPFM sampled-SA final prefix mismatch."
        );
    }

    if (
        index.sampled_sa_values_
            .empty()
    ) {
        throw std::runtime_error(
            "PPFM sampled-SA values are empty."
        );
    }

    for (const auto value :
         index.sampled_sa_values_) {

        if (
            static_cast<std::uint64_t>(
                value
            ) >=
            expected_length
        ) {
            throw std::runtime_error(
                "PPFM sampled-SA value is out of range."
            );
        }

        if (
            static_cast<std::size_t>(
                value
            ) %
                index
                    .suffix_array_sample_rate_
            !=
            0
        ) {
            throw std::runtime_error(
                "PPFM sampled-SA value violates sampling rate."
            );
        }
    }
}


void PpfmIO::validate_bidirectional_index(
    const BidirectionalFMIndex& index,
    const std::uint64_t reference_length
) {
    validate_fm_index(
        index.forward_index_,
        reference_length
    );

    validate_fm_index(
        index.reverse_index_,
        reference_length
    );

    if (
        index.forward_index_
            .suffix_array_sample_rate_
        !=
        index.reverse_index_
            .suffix_array_sample_rate_
    ) {
        throw std::runtime_error(
            "PPFM forward/reverse SA sample rates differ."
        );
    }
}


void PpfmIO::save_shard(
    const std::filesystem::path& path,
    const std::string& chromosome,
    const PackedReference& reference,
    const BidirectionalFMIndex& index
) {
    if (chromosome.empty()) {
        throw std::invalid_argument(
            "PPFM chromosome name cannot be empty."
        );
    }

    validate_reference(
        reference
    );

    validate_bidirectional_index(
        index,
        reference.length_
    );

    std::ofstream output(
        path,
        std::ios::binary |
        std::ios::trunc
    );

    if (!output) {
        throw std::runtime_error(
            "Cannot open PPFM file for writing: " +
            path.string()
        );
    }

    output.write(
        kPpfmMagic.data(),
        static_cast<std::streamsize>(
            kPpfmMagic.size()
        )
    );

    require_write(output);

    write_u32(
        output,
        kFormatVersion
    );

    write_string(
        output,
        chromosome
    );

    write_reference(
        output,
        reference
    );

    write_fm_index(
        output,
        index.forward_index_
    );

    write_fm_index(
        output,
        index.reverse_index_
    );

    output.flush();

    require_write(output);
}


PpfmShardData PpfmIO::load_shard(
    const std::filesystem::path& path
) {
    std::ifstream input(
        path,
        std::ios::binary
    );

    if (!input) {
        throw std::runtime_error(
            "Cannot open PPFM file for reading: " +
            path.string()
        );
    }

    std::array<char, 8> magic{};

    input.read(
        magic.data(),
        static_cast<std::streamsize>(
            magic.size()
        )
    );

    require_read(input);

    if (magic != kPpfmMagic) {
        throw std::runtime_error(
            "Invalid PPFM magic."
        );
    }

    const std::uint32_t version =
        read_u32(input);

    if (
        version !=
        kFormatVersion
    ) {
        throw std::runtime_error(
            "Unsupported PPFM format version."
        );
    }

    std::string chromosome =
        read_string(input);

    if (chromosome.empty()) {
        throw std::runtime_error(
            "PPFM chromosome name is empty."
        );
    }

    /*
     * Tiny valid objects are constructed first.
     *
     * Their build state is then replaced directly
     * with validated persistent state. This avoids
     * adding public raw-state constructors to the
     * search classes.
     */
    PackedReference reference{
        std::string_view{"A"}
    };

    BidirectionalFMIndex index{
        std::string{"A"},
        std::size_t{1}
    };

    read_reference(
        input,
        reference
    );

    read_fm_index(
        input,
        index.forward_index_
    );

    read_fm_index(
        input,
        index.reverse_index_
    );

    validate_reference(
        reference
    );

    validate_bidirectional_index(
        index,
        reference.length_
    );

    char trailing = '\0';

    if (
        input.read(
            &trailing,
            1
        )
    ) {
        throw std::runtime_error(
            "PPFM file contains unexpected trailing data."
        );
    }

    if (!input.eof()) {
        throw std::runtime_error(
            "PPFM stream ended in an invalid state."
        );
    }

    return PpfmShardData{
        std::move(chromosome),
        std::move(reference),
        std::move(index)
    };
}

}  // namespace primerpair
