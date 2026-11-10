#pragma once

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <string>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"

namespace primerpair {

/*
 * Fully restored persistent chromosome/shard state.
 *
 * No FASTA, suffix-array construction or BWT
 * construction is required after load_shard().
 */
struct PpfmShardData {
    std::string chromosome;
    PackedReference reference;
    BidirectionalFMIndex index;
};


class PpfmIO {
public:
    static constexpr std::uint32_t
        kFormatVersion = 1;

    static void save_shard(
        const std::filesystem::path& path,
        const std::string& chromosome,
        const PackedReference& reference,
        const BidirectionalFMIndex& index
    );

    [[nodiscard]]
    static PpfmShardData load_shard(
        const std::filesystem::path& path
    );

private:
    static void write_reference(
        std::ostream& output,
        const PackedReference& reference
    );

    static void read_reference(
        std::istream& input,
        PackedReference& reference
    );

    static void write_packed_bwt(
        std::ostream& output,
        const PackedBWT& bwt
    );

    static void read_packed_bwt(
        std::istream& input,
        PackedBWT& bwt
    );

    static void write_rank(
        std::ostream& output,
        const CheckpointRank& rank
    );

    static void read_rank(
        std::istream& input,
        CheckpointRank& rank
    );

    static void write_fm_index(
        std::ostream& output,
        const FMIndex& index
    );

    static void read_fm_index(
        std::istream& input,
        FMIndex& index
    );

    static void validate_reference(
        const PackedReference& reference
    );

    static void validate_fm_index(
        const FMIndex& index,
        std::uint64_t reference_length
    );

    static void validate_bidirectional_index(
        const BidirectionalFMIndex& index,
        std::uint64_t reference_length
    );
};

}  // namespace primerpair
