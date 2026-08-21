#include "primerpair/anchor_search.hpp"

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

bool same_interval(
    const primerpair::Interval& lhs,
    const primerpair::Interval& rhs
) {
    if (lhs.empty() && rhs.empty()) {
        return true;
    }

    return
        lhs.begin == rhs.begin &&
        lhs.end == rhs.end;
}

}  // namespace

int main() {
    try {
        const std::string reference =
            "GGGACGTTGCAACGTACGTCCC"
            "TTTACGTTGCAACGTACGTAAA";

        const primerpair::BidirectionalFMIndex index(
            reference
        );

        const primerpair::AnchorSearcher searcher(
            index
        );

        const std::string primer =
            "ACGTTGCAACGTACGT";

        const auto anchored =
            searcher.search_exact(
                primer,
                12
            );

        const auto direct =
            index.search(
                primer
            );

        expect(
            anchored.anchor_length == 12,
            "3-prime anchor length"
        );

        expect(
            anchored.primer_length ==
                primer.size(),
            "Primer length"
        );

        expect(
            anchored.extension_steps ==
                primer.size() - 12,
            "5-prime extension count"
        );

        expect(
            same_interval(
                anchored.state.forward,
                direct.forward
            ),
            "Anchor-first forward equivalence"
        );

        expect(
            same_interval(
                anchored.state.reverse,
                direct.reverse
            ),
            "Anchor-first reverse equivalence"
        );

        expect(
            anchored.match_count() ==
                direct.size(),
            "Anchor-first occurrence count"
        );

        auto anchored_positions =
            searcher.locate(
                anchored
            );

        auto direct_positions =
            index.locate(
                direct
            );

        std::sort(
            anchored_positions.begin(),
            anchored_positions.end()
        );

        std::sort(
            direct_positions.begin(),
            direct_positions.end()
        );

        expect(
            anchored_positions ==
                direct_positions,
            "Anchor-first genomic positions"
        );

        /*
         * Primerin 5' tarafını bozuyoruz.
         * Anchor aynı kalabilir ama full primer
         * extension sırasında elenmelidir.
         */
        std::string absent =
            primer;

        absent.at(0) =
            absent.at(0) == 'A'
                ? 'C'
                : 'A';

        const auto absent_result =
            searcher.search_exact(
                absent,
                12
            );

        const auto absent_direct =
            index.search(
                absent
            );

        expect(
            absent_result.match_count() ==
                absent_direct.size(),
            "Absent primer equivalence"
        );

        bool zero_anchor_rejected = false;

        try {
            static_cast<void>(
                searcher.search_exact(
                    primer,
                    0
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            zero_anchor_rejected = true;
        }

        expect(
            zero_anchor_rejected,
            "Zero anchor rejection"
        );

        bool long_anchor_rejected = false;

        try {
            static_cast<void>(
                searcher.search_exact(
                    primer,
                    primer.size() + 1
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            long_anchor_rejected = true;
        }

        expect(
            long_anchor_rejected,
            "Oversized anchor rejection"
        );

        std::cout
            << "All exact 3-prime anchor search tests passed.\n";

        return 0;

    } catch (const std::exception& error) {
        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
