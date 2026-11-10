#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

std::string
make_reference(
    const std::size_t length
) {
    constexpr char alphabet[] = {
        'A',
        'C',
        'G',
        'T'
    };

    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;

    std::string reference;
    reference.reserve(length);

    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL
            +
            1442695040888963407ULL;

        reference.push_back(
            alphabet[
                static_cast<std::size_t>(
                    (state >> 32U) & 3ULL
                )
            ]
        );
    }

    /*
     * Deliberate high-copy regions.
     */
    const std::string motif =
        "ACGTACGTACGTACGT"
        "AAAAAAAACCCCCCCC"
        "GGGGGGGGTTTTTTTT"
        "ACACACACGTGTGTGT";

    for (
        std::size_t start = 4096;
        start + 512 < reference.size();
        start += 8192
    ) {
        for (
            std::size_t offset = 0;
            offset < 512;
            ++offset
        ) {
            reference[start + offset] =
                motif[
                    offset %
                    motif.size()
                ];
        }
    }

    return reference;
}


char
complement(
    const char base
) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
    }

    throw std::invalid_argument(
        "Invalid DNA base."
    );
}


std::string
reverse_complement(
    const std::string_view sequence
) {
    std::string result;
    result.reserve(
        sequence.size()
    );

    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {
        result.push_back(
            complement(*it)
        );
    }

    return result;
}


void
normalize_positions(
    std::vector<std::uint64_t>& positions
) {
    std::sort(
        positions.begin(),
        positions.end()
    );

    positions.erase(
        std::unique(
            positions.begin(),
            positions.end()
        ),
        positions.end()
    );
}


void
compare_anchor(
    const primerpair::BidirectionalFMIndex& bifm,
    const primerpair::IPBWTIndex& ipbwt,
    const std::string_view anchor,
    const std::size_t check_id,
    std::size_t& checks
) {
    const auto bifm_interval =
        bifm.search(
            anchor
        );

    const auto ip_interval =
        ipbwt.exact_search(
            anchor
        );


    if (
        bifm_interval.size() !=
        ip_interval.size()
    ) {
        std::cerr
            << "SIZE_MISMATCH"
            << '\t'
            << "check="
            << check_id
            << '\t'
            << "anchor="
            << anchor
            << '\t'
            << "bifm="
            << bifm_interval.size()
            << '\t'
            << "ipbwt="
            << ip_interval.size()
            << '\n';

        throw std::runtime_error(
            "BiFM/IP-BWT anchor interval "
            "size mismatch."
        );
    }


    auto bifm_positions =
        bifm.locate(
            bifm_interval
        );

    auto ip_positions =
        ipbwt.locate(
            ip_interval
        );


    normalize_positions(
        bifm_positions
    );

    normalize_positions(
        ip_positions
    );


    if (
        bifm_positions !=
        ip_positions
    ) {
        std::cerr
            << "POSITION_MISMATCH"
            << '\t'
            << "check="
            << check_id
            << '\t'
            << "anchor="
            << anchor
            << '\t'
            << "bifm_hits="
            << bifm_positions.size()
            << '\t'
            << "ipbwt_hits="
            << ip_positions.size()
            << '\n';

        throw std::runtime_error(
            "BiFM/IP-BWT anchor genomic "
            "positions differ."
        );
    }


    ++checks;
}

}  // namespace


int
main() {
    try {
        constexpr std::size_t
            anchor_length = 12;

        const std::string reference =
            make_reference(
                200000
            );


        primerpair::BidirectionalFMIndex
            bifm(
                reference
            );

        primerpair::IPBWTIndex
            ipbwt(
                reference,
                anchor_length
            );


        std::size_t checks = 0;
        std::size_t check_id = 0;


        /*
         * Explicit low-complexity/repeat-rich
         * anchors.
         */
        for (
            const std::string anchor :
            {
                "AAAAAAAAAAAA",
                "CCCCCCCCCCCC",
                "GGGGGGGGGGGG",
                "TTTTTTTTTTTT",
                "ACGTACGTACGT",
                "ACACACACACAC",
                "GTGTGTGTGTGT"
            }
        ) {
            compare_anchor(
                bifm,
                ipbwt,
                anchor,
                check_id++,
                checks
            );
        }


        /*
         * Present genomic anchors plus
         * deterministic one-base mutations.
         */
        for (
            std::size_t i = 0;
            i < 2500;
            ++i
        ) {
            const std::size_t position =
                (
                    i *
                    7919
                    +
                    123
                )
                %
                (
                    reference.size()
                    -
                    anchor_length
                );


            std::string anchor =
                reference.substr(
                    position,
                    anchor_length
                );


            compare_anchor(
                bifm,
                ipbwt,
                anchor,
                check_id++,
                checks
            );


            /*
             * May be present or absent; either is
             * valid. The two indexes must agree.
             */
            anchor[5] =
                (
                    anchor[5] ==
                    'A'
                )
                ?
                'T'
                :
                'A';


            compare_anchor(
                bifm,
                ipbwt,
                anchor,
                check_id++,
                checks
            );
        }


        /*
         * Explicitly exercise the same biological
         * 3-prime semantics used by the primer
         * engines.
         *
         * Forward:
         *     exact anchor = primer suffix
         *
         * Reverse:
         *     reverse_query = RC(primer)
         *     exact anchor = reverse_query prefix
         */
        for (
            std::size_t i = 0;
            i < 1000;
            ++i
        ) {
            constexpr std::size_t
                primer_length = 24;

            const std::size_t position =
                (
                    i *
                    3571
                    +
                    2048
                )
                %
                (
                    reference.size()
                    -
                    primer_length
                );


            const std::string primer =
                reference.substr(
                    position,
                    primer_length
                );


            const std::string_view
                forward_anchor(
                    primer.data()
                        +
                        primer.size()
                        -
                        anchor_length,
                    anchor_length
                );


            compare_anchor(
                bifm,
                ipbwt,
                forward_anchor,
                check_id++,
                checks
            );


            const std::string reverse_query =
                reverse_complement(
                    primer
                );


            const std::string_view
                reverse_anchor(
                    reverse_query.data(),
                    anchor_length
                );


            compare_anchor(
                bifm,
                ipbwt,
                reverse_anchor,
                check_id++,
                checks
            );
        }


        std::cout
            << "anchor_length\t"
            << anchor_length
            << '\n';

        std::cout
            << "checks\t"
            << checks
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
