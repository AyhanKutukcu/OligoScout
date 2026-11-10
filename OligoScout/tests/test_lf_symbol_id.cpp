#include "primerpair/fm_index.hpp"
#include "primerpair/packed_bwt.hpp"
#include "primerpair/rank_support.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::size_t symbol_id(
    const char symbol
) {
    switch (symbol) {
        case '$': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'N': return 4;
        case 'T': return 5;

        default:
            throw std::logic_error(
                "Unexpected test symbol."
            );
    }
}

}  // namespace

int main() {
    try {
        using namespace primerpair;

        constexpr std::array<char, 6>
            alphabet{
                '$',
                'A',
                'C',
                'G',
                'N',
                'T'
            };

        std::string text;
        text.reserve(
            12000
        );

        std::uint64_t rng =
            0xA0761D6478BD642FULL;

        constexpr std::array<char, 4>
            dna{
                'A',
                'C',
                'G',
                'T'
            };

        for (
            std::size_t i = 0;
            i < 10000;
            ++i
        ) {
            rng =
                rng *
                6364136223846793005ULL +
                1442695040888963407ULL;

            if (
                i != 0 &&
                i % 257 == 0
            ) {
                text.push_back(
                    'N'
                );
            } else {
                text.push_back(
                    dna.at(
                        static_cast<std::size_t>(
                            (rng >> 32U) &
                            3ULL
                        )
                    )
                );
            }
        }

        FMIndex index(
            text,
            8
        );

        const std::string bwt =
            index.bwt_string();

        PackedBWT packed(
            bwt
        );

        CheckpointRank rank;

        rank.build(
            bwt,
            CheckpointRank::
                kDefaultCheckpointRate
        );

        std::size_t symbol_checks = 0;
        std::size_t count_checks = 0;
        std::size_t rank_checks = 0;
        std::size_t lf_checks = 0;

        /*
         * 1. symbol_id_at() vs decoded BWT.
         */
        for (
            std::size_t position = 0;
            position < bwt.size();
            ++position
        ) {
            const std::size_t expected =
                symbol_id(
                    bwt.at(position)
                );

            const std::size_t observed =
                packed.symbol_id_at(
                    position
                );

            if (
                expected != observed
            ) {
                std::cerr
                    << "symbol_id mismatch position="
                    << position
                    << '\n';

                return 1;
            }

            ++symbol_checks;
        }

        /*
         * 2. Numeric count vs existing char count.
         */
        for (
            std::size_t trial = 0;
            trial < 5000;
            ++trial
        ) {
            rng =
                rng *
                6364136223846793005ULL +
                1442695040888963407ULL;

            std::uint64_t begin =
                rng %
                (
                    bwt.size() +
                    1
                );

            rng =
                rng *
                6364136223846793005ULL +
                1442695040888963407ULL;

            std::uint64_t end =
                rng %
                (
                    bwt.size() +
                    1
                );

            if (begin > end) {
                std::swap(
                    begin,
                    end
                );
            }

            for (
                std::size_t id = 0;
                id < alphabet.size();
                ++id
            ) {
                const auto old_count =
                    packed.count(
                        alphabet.at(id),
                        begin,
                        end
                    );

                const auto new_count =
                    packed.count_by_symbol_id(
                        id,
                        begin,
                        end
                    );

                if (
                    old_count !=
                    new_count
                ) {
                    std::cerr
                        << "count mismatch trial="
                        << trial
                        << " id="
                        << id
                        << '\n';

                    return 1;
                }

                ++count_checks;
            }
        }

        /*
         * 3. Numeric rank vs existing char rank.
         */
        for (
            std::size_t position = 0;
            position <= bwt.size();
            ++position
        ) {
            for (
                std::size_t id = 0;
                id < alphabet.size();
                ++id
            ) {
                const auto old_rank =
                    rank.rank(
                        packed,
                        alphabet.at(id),
                        position
                    );

                const auto new_rank =
                    rank.rank_by_symbol_id(
                        packed,
                        id,
                        position
                    );

                if (
                    old_rank !=
                    new_rank
                ) {
                    std::cerr
                        << "rank mismatch position="
                        << position
                        << " id="
                        << id
                        << '\n';

                    return 1;
                }

                ++rank_checks;
            }
        }

        /*
         * 4. Mevcut FMIndex::lf()yi bağımsız
         *    BWT-prefix referansıyla doğrula.
         */
        std::array<std::uint64_t, 6>
            totals{};

        for (const char symbol : bwt) {
            ++totals.at(
                symbol_id(symbol)
            );
        }

        std::array<std::uint64_t, 6>
            c_table{};

        std::uint64_t cumulative = 0;

        for (
            std::size_t id = 0;
            id < c_table.size();
            ++id
        ) {
            c_table.at(id) =
                cumulative;

            cumulative +=
                totals.at(id);
        }

        std::array<std::uint64_t, 6>
            prefix{};

        for (
            std::size_t row = 0;
            row < bwt.size();
            ++row
        ) {
            const std::size_t id =
                symbol_id(
                    bwt.at(row)
                );

            const std::uint64_t expected_lf =
                c_table.at(id) +
                prefix.at(id);

            const std::uint64_t observed_lf =
                index.lf(
                    row
                );

            if (
                expected_lf !=
                observed_lf
            ) {
                std::cerr
                    << "LF mismatch row="
                    << row
                    << " expected="
                    << expected_lf
                    << " observed="
                    << observed_lf
                    << '\n';

                return 1;
            }

            ++prefix.at(id);
            ++lf_checks;
        }

        std::cout
            << "symbol_checks\t"
            << symbol_checks
            << '\n';

        std::cout
            << "count_checks\t"
            << count_checks
            << '\n';

        std::cout
            << "rank_checks\t"
            << rank_checks
            << '\n';

        std::cout
            << "lf_checks\t"
            << lf_checks
            << '\n';

        std::cout
            << "ALL_CHECKS\tYES\n";

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';

        return 1;
    }
}
