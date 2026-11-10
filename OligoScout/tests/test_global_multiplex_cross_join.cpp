#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/global_multiplex_cross_join.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {


std::string make_reference(
    const std::size_t length
) {
    constexpr char bases[] = {
        'A', 'C', 'G', 'T'
    };

    std::uint64_t state =
        0xE7037ED1A0B428DBULL;

    std::string reference;

    reference.reserve(
        length
    );


    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        reference.push_back(
            bases[
                static_cast<std::size_t>(
                    (state >> 32U) &
                    3ULL
                )
            ]
        );
    }


    /*
     * Controlled exact repeat.
     *
     * Böylece global sweep yalnız unique single-hit
     * primerlerde değil multi-hit primerlerde de
     * doğrulanır.
     */
    const std::string repeat =
        "AGTCCGATGCTAACGTTGACCTGA";


    const std::vector<std::size_t>
        repeat_positions{
            5000,
            20000,
            35000,
            50000,
            65000,
            80000
        };


    for (
        const std::size_t position :
        repeat_positions
    ) {
        reference.replace(
            position,
            repeat.size(),
            repeat
        );
    }


    return reference;
}


bool product_less(
    const primerpair::
        GlobalMultiplexCrossProduct& lhs,

    const primerpair::
        GlobalMultiplexCrossProduct& rhs
) noexcept {
    return
        std::tie(
            lhs.amplicon_start,
            lhs.amplicon_end_exclusive,

            lhs.forward_primer_id,
            lhs.reverse_primer_id,

            lhs.forward_position,
            lhs.reverse_position,

            lhs.forward_mismatches,
            lhs.reverse_mismatches,

            lhs.forward_mismatch_mask,
            lhs.reverse_mismatch_mask
        )
        <
        std::tie(
            rhs.amplicon_start,
            rhs.amplicon_end_exclusive,

            rhs.forward_primer_id,
            rhs.reverse_primer_id,

            rhs.forward_position,
            rhs.reverse_position,

            rhs.forward_mismatches,
            rhs.reverse_mismatches,

            rhs.forward_mismatch_mask,
            rhs.reverse_mismatch_mask
        );
}


void normalize(
    std::vector<
        primerpair::
            GlobalMultiplexCrossProduct
    >& products
) {
    std::sort(
        products.begin(),
        products.end(),
        product_less
    );

    products.erase(
        std::unique(
            products.begin(),
            products.end()
        ),
        products.end()
    );
}


/*
 * Independent brute-force reference.
 *
 * Global sweep kullanmaz.
 *
 * Every forward hit × every reverse hit.
 *
 * Complexity:
 *
 * O(Hf * Hr)
 */
std::vector<
    primerpair::
        GlobalMultiplexCrossProduct
>
brute_force_join(
    const std::vector<
        primerpair::
            GlobalMultiplexPrimerHits
    >& primers,

    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon
) {
    using namespace primerpair;


    std::vector<
        GlobalMultiplexCrossProduct
    > products;


    for (
        const auto& forward_primer :
        primers
    ) {
        for (
            const auto& forward_hit :
            forward_primer.hits
        ) {
            if (
                forward_hit.orientation !=
                PrimerOrientation::Forward
            ) {
                continue;
            }


            if (
                forward_hit.position >
                std::numeric_limits<
                    std::uint64_t
                >::max() -
                    forward_primer
                        .primer_length
            ) {
                continue;
            }


            const std::uint64_t
                minimum_reverse_position =
                    forward_hit.position +
                    static_cast<
                        std::uint64_t
                    >(
                        forward_primer
                            .primer_length
                    );


            for (
                const auto& reverse_primer :
                primers
            ) {
                for (
                    const auto& reverse_hit :
                    reverse_primer.hits
                ) {
                    if (
                        reverse_hit.orientation !=
                        PrimerOrientation::Reverse
                    ) {
                        continue;
                    }


                    if (
                        reverse_hit.position <
                        minimum_reverse_position
                    ) {
                        continue;
                    }


                    if (
                        reverse_hit.position >
                        std::numeric_limits<
                            std::uint64_t
                        >::max() -
                            reverse_primer
                                .primer_length
                    ) {
                        continue;
                    }


                    const std::uint64_t
                        reverse_end =
                            reverse_hit.position +
                            static_cast<
                                std::uint64_t
                            >(
                                reverse_primer
                                    .primer_length
                            );


                    if (
                        reverse_end <=
                        forward_hit.position
                    ) {
                        continue;
                    }


                    const std::uint64_t
                        amplicon_length =
                            reverse_end -
                            forward_hit.position;


                    if (
                        amplicon_length <
                            min_amplicon
                        ||
                        amplicon_length >
                            max_amplicon
                    ) {
                        continue;
                    }


                    products.push_back(
                        GlobalMultiplexCrossProduct{
                            forward_primer
                                .primer_id,

                            reverse_primer
                                .primer_id,

                            forward_hit
                                .position,

                            reverse_hit
                                .position,

                            forward_hit
                                .position,

                            reverse_end,

                            amplicon_length,

                            forward_hit
                                .mismatches,

                            reverse_hit
                                .mismatches,

                            forward_hit
                                .mismatch_mask,

                            reverse_hit
                                .mismatch_mask
                        }
                    );
                }
            }
        }
    }


    normalize(
        products
    );


    return products;
}


}  // namespace


int main() {
    try {
        using namespace primerpair;


        constexpr std::size_t
            anchor_length = 12;


        const std::string reference =
            make_reference(
                100000
            );


        PackedReference packed(
            reference
        );


        BidirectionalFMIndex bifm(
            reference
        );


        StrandAwarePrimerSearchEngine
            strand_engine(
                bifm,
                packed
            );


        /*
         * 46 normal primers +
         * repeat forward +
         * repeat reverse
         *
         * = 48 unique logical primer IDs.
         */
        constexpr std::size_t
            ordinary_primer_count = 46;


        std::vector<std::string>
            primer_storage;

        primer_storage.reserve(
            ordinary_primer_count + 2
        );


        for (
            std::size_t i = 0;
            i < ordinary_primer_count;
            ++i
        ) {
            /*
             * Different lengths:
             *
             * 18..30 nt
             */
            const std::size_t length =
                18 +
                (
                    i %
                    13
                );


            /*
             * Primers clustered enough that multiple
             * PCR-compatible cross products exist.
             */
            const std::size_t position =
                10000 +
                i * 140;


            const std::string genomic =
                reference.substr(
                    position,
                    length
                );


            /*
             * Alternating orientation.
             *
             * Even IDs:
             *   exact forward genomic primer.
             *
             * Odd IDs:
             *   reverse-complemented genomic primer.
             *
             * Thus the panel contains both forward
             * and reverse binding sites.
             */
            if (
                i %
                2 ==
                0
            ) {
                primer_storage.push_back(
                    genomic
                );

            } else {
                primer_storage.push_back(
                    reverse_complement(
                        genomic
                    )
                );
            }
        }


        const std::string repeat =
            "AGTCCGATGCTAACGTTGACCTGA";


        primer_storage.push_back(
            repeat
        );


        primer_storage.push_back(
            reverse_complement(
                repeat
            )
        );


        std::vector<
            std::vector<
                OrientedPrimerSearchHit
            >
        > hit_storage;


        hit_storage.reserve(
            primer_storage.size()
        );


        for (
            const auto& primer :
            primer_storage
        ) {
            const auto search_result =
                strand_engine.search(
                    primer,
                    anchor_length,
                    0
                );


            hit_storage.push_back(
                search_result.hits
            );
        }


        /*
         * Build spans only after hit_storage is fully
         * populated, so all span targets remain stable.
         */
        std::vector<
            GlobalMultiplexPrimerHits
        > primer_hits;


        primer_hits.reserve(
            primer_storage.size()
        );


        for (
            std::size_t id = 0;
            id < primer_storage.size();
            ++id
        ) {
            primer_hits.push_back(
                GlobalMultiplexPrimerHits{
                    id,
                    primer_storage
                        .at(id)
                        .size(),
                    hit_storage.at(id)
                }
            );
        }


        struct Window {
            std::uint64_t minimum;
            std::uint64_t maximum;
        };


        const std::vector<Window>
            windows{
                {50, 300},
                {100, 1000},
                {500, 2500},
                {50, 5000}
            };


        std::size_t checks = 0;

        std::size_t total_products = 0;

        std::size_t total_window_candidates = 0;


        for (
            const auto& window :
            windows
        ) {
            const auto observed =
                global_multiplex_cross_join(
                    primer_hits,
                    window.minimum,
                    window.maximum
                );


            auto expected =
                brute_force_join(
                    primer_hits,
                    window.minimum,
                    window.maximum
                );


            auto actual =
                observed.products;


            normalize(
                expected
            );

            normalize(
                actual
            );


            if (
                expected !=
                actual
            ) {
                std::cerr
                    << "GLOBAL_SWEEP_MISMATCH"
                    << '\t'
                    << "min="
                    << window.minimum
                    << '\t'
                    << "max="
                    << window.maximum
                    << '\t'
                    << "expected="
                    << expected.size()
                    << '\t'
                    << "actual="
                    << actual.size()
                    << '\n';

                return 1;
            }


            if (
                observed
                    .stats
                    .unique_products !=
                actual.size()
            ) {
                throw std::runtime_error(
                    "Global sweep unique-product "
                    "stat mismatch."
                );
            }


            total_products +=
                actual.size();


            total_window_candidates +=
                observed
                    .stats
                    .window_candidates;


            std::cout
                << "window\t"
                << window.minimum
                << '-'
                << window.maximum
                << '\t'
                << "products="
                << actual.size()
                << '\t'
                << "candidates="
                << observed
                    .stats
                    .window_candidates
                << '\n';


            ++checks;
        }


        std::size_t forward_hits = 0;
        std::size_t reverse_hits = 0;


        for (
            const auto& hits :
            hit_storage
        ) {
            for (
                const auto& hit :
                hits
            ) {
                if (
                    hit.orientation ==
                    PrimerOrientation::Forward
                ) {
                    ++forward_hits;

                } else {
                    ++reverse_hits;
                }
            }
        }


        std::cout
            << "primers\t"
            << primer_storage.size()
            << '\n';


        std::cout
            << "forward_hits\t"
            << forward_hits
            << '\n';


        std::cout
            << "reverse_hits\t"
            << reverse_hits
            << '\n';


        std::cout
            << "windows_checked\t"
            << checks
            << '\n';


        std::cout
            << "total_products\t"
            << total_products
            << '\n';


        std::cout
            << "total_window_candidates\t"
            << total_window_candidates
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
