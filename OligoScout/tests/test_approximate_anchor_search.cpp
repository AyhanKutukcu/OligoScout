#include "primerpair/approximate_anchor_search.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect(
    const bool condition,
    const std::string& name
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " + name
        );
    }

    std::cout
        << "[PASS] "
        << name
        << '\n';
}

std::vector<std::uint64_t>
collect_positions(
    const primerpair::ApproximateAnchorSearcher& searcher,
    const primerpair::ApproximateAnchorSearchResult& result
) {
    std::vector<std::uint64_t>
        positions;

    for (const auto& hit : result.hits) {
        auto current =
            searcher.locate(
                hit
            );

        positions.insert(
            positions.end(),
            current.begin(),
            current.end()
        );
    }

    std::sort(
        positions.begin(),
        positions.end()
    );

    return positions;
}

}  // namespace

int main() {
    try {
        /*
         * Primer:
         *
         * ACGT | TGCAACGTACGT
         * ----   ------------
         * 5'     exact 12-nt 3' anchor
         */

        const std::string primer =
            "ACGTTGCAACGTACGT";

        /*
         * 0 mismatch
         */
        const std::string exact =
            "ACGTTGCAACGTACGT";

        /*
         * Primer position 0:
         * A -> C
         */
        const std::string one_mismatch =
            "CCGTTGCAACGTACGT";

        /*
         * Primer positions 0 and 1:
         * AC -> TA
         */
        const std::string two_mismatch =
            "TAGTTGCAACGTACGT";

        /*
         * 3' anchor içinde mismatch.
         *
         * Exact-anchor profilinde bulunmamalı.
         */
        const std::string anchor_mismatch =
            "ACGTTGCAACGTACGA";

        const std::string reference =
            "TTT" +
            exact +
            "GGG" +
            one_mismatch +
            "CCC" +
            two_mismatch +
            "GGG" +
            anchor_mismatch +
            "AAA";

        const primerpair::BidirectionalFMIndex index(
            reference
        );

        const primerpair::ApproximateAnchorSearcher
            searcher(
                index
            );

        /*
         * --------------------------------------------------
         * max mismatch = 0
         * --------------------------------------------------
         */
        const auto zero =
            searcher.search_5prime_mismatches(
                primer,
                12,
                0
            );

        expect(
            zero.total_match_count() == 1,
            "Zero-mismatch occurrence count"
        );

        expect(
            collect_positions(
                searcher,
                zero
            ) ==
                std::vector<std::uint64_t>{
                    3
                },
            "Zero-mismatch genomic position"
        );

        /*
         * --------------------------------------------------
         * max mismatch = 1
         * --------------------------------------------------
         */
        const auto one =
            searcher.search_5prime_mismatches(
                primer,
                12,
                1
            );

        expect(
            one.total_match_count() == 2,
            "One-mismatch occurrence count"
        );

        expect(
            collect_positions(
                searcher,
                one
            ) ==
                std::vector<std::uint64_t>{
                    3,
                    22
                },
            "One-mismatch genomic positions"
        );

        bool found_position_zero_mismatch =
            false;

        for (const auto& hit : one.hits) {
            if (
                hit.mismatches == 1 &&
                hit.mismatch_positions ==
                    std::vector<std::size_t>{
                        0
                    }
            ) {
                found_position_zero_mismatch =
                    true;
            }
        }

        expect(
            found_position_zero_mismatch,
            "Mismatch position tracking"
        );

        /*
         * --------------------------------------------------
         * max mismatch = 2
         * --------------------------------------------------
         */
        const auto two =
            searcher.search_5prime_mismatches(
                primer,
                12,
                2
            );

        expect(
            two.total_match_count() == 3,
            "Two-mismatch occurrence count"
        );

        expect(
            collect_positions(
                searcher,
                two
            ) ==
                std::vector<std::uint64_t>{
                    3,
                    22,
                    41
                },
            "Two-mismatch genomic positions"
        );

        /*
         * Anchor mismatch içeren sekansın genomik
         * başlangıcı:
         *
         * 60
         *
         * max_mismatches=3 olsa bile exact 3' anchor
         * kuralı nedeniyle sonuçta bulunmamalı.
         */
        const auto three =
            searcher.search_5prime_mismatches(
                primer,
                12,
                3
            );

        const auto positions_three =
            collect_positions(
                searcher,
                three
            );

        expect(
            std::find(
                positions_three.begin(),
                positions_three.end(),
                60
            ) ==
                positions_three.end(),
            "Exact 3-prime anchor enforcement"
        );

        /*
         * Lowercase primer de aynı sonucu vermeli.
         */
        const auto lowercase =
            searcher.search_5prime_mismatches(
                "acgttgcaacgtacgt",
                12,
                1
            );

        expect(
            collect_positions(
                searcher,
                lowercase
            ) ==
                std::vector<std::uint64_t>{
                    3,
                    22
                },
            "Lowercase approximate primer"
        );

        bool excessive_mismatch_rejected =
            false;

        try {
            static_cast<void>(
                searcher.search_5prime_mismatches(
                    primer,
                    12,
                    4
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            excessive_mismatch_rejected =
                true;
        }

        expect(
            excessive_mismatch_rejected,
            "More-than-three mismatch rejection"
        );

        std::cout
            << "All 3-prime-aware approximate anchor tests passed.\n";

        return 0;

    } catch (const std::exception& error) {
        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
