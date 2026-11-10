#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <cstddef>
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

}  // namespace

int main() {
    try {
        const std::string primer =
            "ACGTTGCAACGTACGT";

        const std::string reverse =
            primerpair::reverse_complement(
                primer
            );

        expect(
            reverse ==
                "ACGTACGTTGCAACGT",
            "Reverse complement"
        );

        /*
         * Exact reverse orientation.
         */
        const std::string reverse_exact =
            reverse;

        /*
         * Mismatch in suffix of reverse_query.
         *
         * This corresponds to the original primer's
         * 5-prime region and is allowed.
         */
        std::string reverse_5prime_mismatch =
            reverse;

        reverse_5prime_mismatch.back() =
            reverse_5prime_mismatch.back() == 'A'
                ? 'C'
                : 'A';

        /*
         * Mismatch in prefix of reverse_query.
         *
         * This corresponds to the original primer's
         * biological 3-prime anchor and MUST be
         * rejected even with k=3.
         */
        std::string reverse_3prime_mismatch =
            reverse;

        reverse_3prime_mismatch.front() =
            reverse_3prime_mismatch.front() == 'A'
                ? 'C'
                : 'A';

        const std::string reference =
            "AAA" +
            primer +
            "CCC" +
            reverse_exact +
            "GGG" +
            reverse_5prime_mismatch +
            "TTT" +
            reverse_3prime_mismatch +
            "AAA";

        /*
         * Expected positions:
         *
         * Forward exact:             3
         * Reverse exact:            22
         * Reverse 5-prime mismatch: 41
         * Reverse 3-prime mismatch: 60 -> reject
         */
        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::StrandAwarePrimerSearchEngine
            engine(
                index,
                packed
            );

        const auto exact =
            engine.search(
                primer,
                12,
                0
            );

        expect(
            exact.hits ==
                std::vector<
                    primerpair::OrientedPrimerSearchHit
                >{
                    {
                        3,
                        0,
                        primerpair::PrimerOrientation::Forward
                    },
                    {
                        22,
                        0,
                        primerpair::PrimerOrientation::Reverse
                    }
                },
            "Exact both-strand hits"
        );

        const auto k1 =
            engine.search(
                primer,
                12,
                1
            );

        expect(
            k1.hits ==
                std::vector<
                    primerpair::OrientedPrimerSearchHit
                >{
                    {
                        3,
                        0,
                        primerpair::PrimerOrientation::Forward
                    },
                    {
                        22,
                        0,
                        primerpair::PrimerOrientation::Reverse
                    },
                    {
                        41,
                        1,
                        primerpair::PrimerOrientation::Reverse,
                        std::uint64_t{1}
                    }
                },
            "Reverse 5-prime mismatch accepted"
        );

        const auto reverse_mismatch_it =
            std::find_if(
                k1.hits.begin(),
                k1.hits.end(),
                [](
                    const auto& hit
                ) {
                    return
                        hit.position == 41 &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                Reverse;
                }
            );

        expect(
            reverse_mismatch_it !=
                k1.hits.end(),
            "Reverse mismatch hit found"
        );

        expect(
            reverse_mismatch_it
                ->mismatch_mask ==
                std::uint64_t{1},
            "Reverse mismatch maps to original 5-prime position 0"
        );

        expect(
            reverse_mismatch_it
                ->mismatches == 1,
            "Mismatch mask count agrees with hit count"
        );


        const auto k3 =
            engine.search(
                primer,
                12,
                3
            );

        const bool bad_anchor_present =
            std::any_of(
                k3.hits.begin(),
                k3.hits.end(),
                [](
                    const auto& hit
                ) {
                    return
                        hit.position == 60 &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                Reverse;
                }
            );

        expect(
            !bad_anchor_present,
            "Reverse biological 3-prime mismatch rejected"
        );

        expect(
            k3.forward_decision
                .max_mismatches == 3,
            "Forward decision stored"
        );

        expect(
            k3.reverse_decision
                .max_mismatches == 3,
            "Reverse decision stored"
        );

        expect(
            std::string(
                primerpair::to_string(
                    primerpair::
                        PrimerOrientation::
                        Forward
                )
            ) ==
                "FORWARD",
            "Forward orientation name"
        );

        expect(
            std::string(
                primerpair::to_string(
                    primerpair::
                        PrimerOrientation::
                        Reverse
                )
            ) ==
                "REVERSE",
            "Reverse orientation name"
        );

        const auto lowercase =
            engine.search(
                "acgttgcaacgtacgt",
                12,
                1
            );

        expect(
            lowercase.hits ==
                k1.hits,
            "Lowercase both-strand search"
        );

        std::cout
            << "All strand-aware primer tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
