#include "primerpair/suffix_array_builder.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace primerpair {

namespace {

constexpr std::size_t kAlphabetSize = 6;


std::size_t symbol_index(
    const char symbol
) {
    switch (symbol) {

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
                "Unsupported suffix-array symbol."
            );
    }
}


void validate_text(
    const std::string_view text
) {
    if (text.empty()) {
        throw std::invalid_argument(
            "Suffix-array text cannot be empty."
        );
    }

    if (
        text.size() >
        static_cast<std::size_t>(
            std::numeric_limits<
                std::uint32_t
            >::max()
        )
    ) {
        throw std::length_error(
            "Suffix-array text exceeds uint32_t "
            "coordinate range."
        );
    }

    if (text.back() != '$') {
        throw std::invalid_argument(
            "Suffix-array text must end with '$'."
        );
    }

    std::size_t sentinel_count = 0;

    for (const char symbol : text) {

        if (symbol == '$') {
            ++sentinel_count;
        }

        (void)symbol_index(
            symbol
        );
    }

    if (sentinel_count != 1) {
        throw std::invalid_argument(
            "Suffix-array text must contain "
            "exactly one sentinel."
        );
    }
}

}  // namespace


std::vector<std::uint32_t>
build_suffix_array_prefix_doubling(
    const std::string_view text
) {
    validate_text(
        text
    );

    const std::size_t text_length =
        text.size();


    std::vector<std::uint32_t>
        suffix_array(
            text_length
        );

    std::iota(
        suffix_array.begin(),
        suffix_array.end(),
        std::uint32_t{0}
    );


    /*
     * Prefix-doubling state.
     *
     * All chromosome positions and equivalence
     * classes stay uint32_t.
     */
    std::vector<std::uint32_t>
        ranks(
            text_length
        );

    std::vector<std::uint32_t>
        next_ranks(
            text_length
        );

    std::vector<std::uint32_t>
        buffer(
            text_length
        );


    for (
        std::size_t position = 0;
        position < text_length;
        ++position
    ) {
        ranks.at(position) =
            static_cast<std::uint32_t>(
                symbol_index(
                    text.at(position)
                )
            );
    }


    std::size_t rank_class_count =
        kAlphabetSize;

    std::vector<std::uint32_t>
        counts;


    for (
        std::size_t width = 1;
        width < text_length;
    ) {

        /*
         * rank + 1 is used for real positions.
         *
         * 0 is reserved for an absent second rank.
         */
        const std::size_t key_count =
            rank_class_count + 1;


        /*
         * ====================================================
         * PASS 1
         *
         * Stable counting sort by second rank.
         * ====================================================
         */

        counts.assign(
            key_count,
            std::uint32_t{0}
        );


        for (
            const std::uint32_t position :
            suffix_array
        ) {
            const std::size_t second_position =
                static_cast<std::size_t>(
                    position
                )
                +
                width;

            const std::size_t key =
                (
                    second_position <
                    text_length
                )
                ?
                static_cast<std::size_t>(
                    ranks.at(
                        second_position
                    )
                ) + 1
                :
                std::size_t{0};

            ++counts.at(
                key
            );
        }


        std::size_t cumulative = 0;

        for (
            std::size_t key = 0;
            key < key_count;
            ++key
        ) {
            const std::uint32_t frequency =
                counts.at(
                    key
                );

            counts.at(
                key
            ) =
                static_cast<std::uint32_t>(
                    cumulative
                );

            cumulative +=
                static_cast<std::size_t>(
                    frequency
                );
        }


        for (
            const std::uint32_t position :
            suffix_array
        ) {
            const std::size_t second_position =
                static_cast<std::size_t>(
                    position
                )
                +
                width;

            const std::size_t key =
                (
                    second_position <
                    text_length
                )
                ?
                static_cast<std::size_t>(
                    ranks.at(
                        second_position
                    )
                ) + 1
                :
                std::size_t{0};


            const std::size_t destination =
                static_cast<std::size_t>(
                    counts.at(
                        key
                    )
                );

            buffer.at(
                destination
            ) =
                position;

            ++counts.at(
                key
            );
        }


        /*
         * ====================================================
         * PASS 2
         *
         * Stable counting sort by first rank.
         * ====================================================
         */

        counts.assign(
            key_count,
            std::uint32_t{0}
        );


        for (
            const std::uint32_t position :
            buffer
        ) {
            const std::size_t key =
                static_cast<std::size_t>(
                    ranks.at(
                        static_cast<std::size_t>(
                            position
                        )
                    )
                ) + 1;

            ++counts.at(
                key
            );
        }


        cumulative = 0;

        for (
            std::size_t key = 0;
            key < key_count;
            ++key
        ) {
            const std::uint32_t frequency =
                counts.at(
                    key
                );

            counts.at(
                key
            ) =
                static_cast<std::uint32_t>(
                    cumulative
                );

            cumulative +=
                static_cast<std::size_t>(
                    frequency
                );
        }


        for (
            const std::uint32_t position :
            buffer
        ) {
            const std::size_t key =
                static_cast<std::size_t>(
                    ranks.at(
                        static_cast<std::size_t>(
                            position
                        )
                    )
                ) + 1;


            const std::size_t destination =
                static_cast<std::size_t>(
                    counts.at(
                        key
                    )
                );

            suffix_array.at(
                destination
            ) =
                position;

            ++counts.at(
                key
            );
        }


        /*
         * ====================================================
         * Recompute equivalence classes.
         * ====================================================
         */

        next_ranks.at(
            static_cast<std::size_t>(
                suffix_array.front()
            )
        ) =
            0;


        std::size_t new_class_count =
            1;


        const auto second_key =
            [
                &ranks,
                text_length,
                width
            ](
                const std::uint32_t position
            ) -> std::uint64_t {

                const std::size_t second_position =
                    static_cast<std::size_t>(
                        position
                    )
                    +
                    width;

                if (
                    second_position >=
                    text_length
                ) {
                    return 0;
                }

                return
                    static_cast<std::uint64_t>(
                        ranks.at(
                            second_position
                        )
                    )
                    +
                    1;
            };


        for (
            std::size_t index = 1;
            index < text_length;
            ++index
        ) {
            const std::uint32_t previous =
                suffix_array.at(
                    index - 1
                );

            const std::uint32_t current =
                suffix_array.at(
                    index
                );


            const bool first_rank_differs =
                ranks.at(
                    static_cast<std::size_t>(
                        previous
                    )
                )
                !=
                ranks.at(
                    static_cast<std::size_t>(
                        current
                    )
                );


            const bool second_rank_differs =
                second_key(
                    previous
                )
                !=
                second_key(
                    current
                );


            if (
                first_rank_differs ||
                second_rank_differs
            ) {
                ++new_class_count;
            }


            next_ranks.at(
                static_cast<std::size_t>(
                    current
                )
            ) =
                static_cast<std::uint32_t>(
                    new_class_count - 1
                );
        }


        ranks.swap(
            next_ranks
        );

        rank_class_count =
            new_class_count;


        if (
            rank_class_count ==
            text_length
        ) {
            break;
        }


        /*
         * Overflow-safe doubling.
         */
        if (
            width >
            text_length / 2
        ) {
            break;
        }

        width *= 2;
    }


    return suffix_array;
}

}  // namespace primerpair
