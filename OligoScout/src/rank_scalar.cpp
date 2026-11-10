#include "primerpair/rank_support.hpp"
#include "primerpair/packed_bwt.hpp"

#include <cctype>
#include <limits>
#include <stdexcept>
#include <string>

namespace primerpair {

std::size_t CheckpointRank::symbol_index(
    char symbol
) {
    const char normalized = static_cast<char>(
        std::toupper(
            static_cast<unsigned char>(symbol)
        )
    );

    switch (normalized) {
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
                std::string(
                    "Unsupported rank symbol: "
                ) + normalized
            );
    }
}

void CheckpointRank::build(
    std::string_view bwt,
    std::size_t checkpoint_rate
) {
    if (checkpoint_rate == 0) {
        throw std::invalid_argument(
            "Checkpoint rate must be greater than zero."
        );
    }

    if (
        bwt.size() >
        static_cast<std::size_t>(
            std::numeric_limits<std::uint32_t>::max()
        )
    ) {
        throw std::length_error(
            "BWT is too large for 32-bit rank checkpoints."
        );
    }

    checkpoint_rate_ = checkpoint_rate;
    checkpoints_.clear();

    checkpoints_.reserve(
        (bwt.size() / checkpoint_rate_) + 2
    );

    std::array<
        std::uint32_t,
        kAlphabetSize
    > cumulative{};

    // BWT pozisyonu 0 için başlangıç sayacı.
    checkpoints_.push_back(cumulative);

    for (std::size_t position = 0;
         position < bwt.size();
         ++position) {

        const std::size_t index =
            symbol_index(bwt.at(position));

        ++cumulative.at(index);

        /*
         * 128, 256, 384... pozisyonlarında
         * kümülatif sayacı kaydet.
         */
        if ((position + 1) % checkpoint_rate_ == 0) {
            checkpoints_.push_back(cumulative);
        }
    }
}

std::uint64_t CheckpointRank::rank(
    std::string_view bwt,
    char symbol,
    std::uint64_t position
) const {
    if (position > bwt.size()) {
        throw std::out_of_range(
            "Rank position exceeds BWT length."
        );
    }

    if (checkpoints_.empty()) {
        throw std::logic_error(
            "Checkpoint rank structure has not been built."
        );
    }

    const std::size_t symbol_id =
        symbol_index(symbol);

    const std::size_t requested_position =
        static_cast<std::size_t>(position);

    const std::size_t checkpoint_id =
        requested_position / checkpoint_rate_;

    if (checkpoint_id >= checkpoints_.size()) {
        throw std::out_of_range(
            "Rank checkpoint does not exist."
        );
    }

    std::uint64_t count =
        checkpoints_
            .at(checkpoint_id)
            .at(symbol_id);

    const std::size_t checkpoint_position =
        checkpoint_id * checkpoint_rate_;

    /*
     * Checkpoint ile sorgu pozisyonu arasında
     * en fazla 127 sembol scalar olarak taranır.
     */
    for (std::size_t current = checkpoint_position;
         current < requested_position;
         ++current) {

        if (symbol_index(bwt.at(current)) ==
            symbol_id) {

            ++count;
        }
    }

    return count;
}

std::uint64_t CheckpointRank::rank_by_symbol_id(
    const PackedBWT& bwt,
    const std::size_t symbol_id,
    const std::uint64_t position
) const {
    if (position > bwt.size()) {
        throw std::out_of_range(
            "Rank position exceeds packed BWT length."
        );
    }

    if (checkpoints_.empty()) {
        throw std::logic_error(
            "Checkpoint rank structure has not been built."
        );
    }

    if (symbol_id >= kAlphabetSize) {
        throw std::invalid_argument(
            "Rank symbol ID exceeds alphabet."
        );
    }

    const std::size_t requested_position =
        static_cast<std::size_t>(
            position
        );

    const std::size_t checkpoint_id =
        requested_position /
        checkpoint_rate_;

    if (
        checkpoint_id >=
        checkpoints_.size()
    ) {
        throw std::out_of_range(
            "Rank checkpoint does not exist."
        );
    }

    std::uint64_t count =
        static_cast<std::uint64_t>(
            checkpoints_
                .at(checkpoint_id)
                .at(symbol_id)
        );

    const std::uint64_t checkpoint_position =
        static_cast<std::uint64_t>(
            checkpoint_id *
            checkpoint_rate_
        );

    count +=
        bwt.count_by_symbol_id(
            symbol_id,
            checkpoint_position,
            position
        );

    return count;
}


std::uint64_t CheckpointRank::rank(
    const PackedBWT& bwt,
    const char symbol,
    const std::uint64_t position
) const {
    return
        rank_by_symbol_id(
            bwt,
            symbol_index(
                symbol
            ),
            position
        );
}


std::array<
    std::uint64_t,
    CheckpointRank::kAlphabetSize
>
CheckpointRank::rank_all(
    const PackedBWT& bwt,
    const std::uint64_t position
) const {
    if (position > bwt.size()) {
        throw std::out_of_range(
            "Rank-all position exceeds packed BWT length."
        );
    }

    if (checkpoints_.empty()) {
        throw std::logic_error(
            "Checkpoint rank structure has not been built."
        );
    }

    const std::size_t requested_position =
        static_cast<std::size_t>(
            position
        );

    const std::size_t checkpoint_id =
        requested_position /
        checkpoint_rate_;

    if (
        checkpoint_id >=
        checkpoints_.size()
    ) {
        throw std::out_of_range(
            "Rank-all checkpoint does not exist."
        );
    }

    std::array<
        std::uint64_t,
        kAlphabetSize
    > counts{};

    for (
        std::size_t i = 0;
        i < kAlphabetSize;
        ++i
    ) {
        counts[i] =
            static_cast<std::uint64_t>(
                checkpoints_
                    .at(checkpoint_id)
                    .at(i)
            );
    }

    const std::uint64_t checkpoint_position =
        static_cast<std::uint64_t>(
            checkpoint_id *
            checkpoint_rate_
        );

    const auto tail =
        bwt.count_all(
            checkpoint_position,
            position
        );

    for (
        std::size_t i = 0;
        i < kAlphabetSize;
        ++i
    ) {
        counts[i] +=
            tail[i];
    }

    return counts;
}


}  // namespace primerpair
