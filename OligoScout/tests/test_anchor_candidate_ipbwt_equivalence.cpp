#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/packed_reference.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
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
            6364136223846793005ULL
            +
            1442695040888963407ULL;


        reference.push_back(
            alphabet[
                static_cast<std::size_t>(
                    (
                        state >>
                        32U
                    )
                    &
                    3ULL
                )
            ]
        );
    }


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
            reference[
                start +
                offset
            ] =
                motif[
                    offset %
                    motif.size()
                ];
        }
    }


    return reference;
}


void
normalize_hits(
    std::vector<primerpair::AnchorCandidateHit>& hits
) {
    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const auto& lhs,
            const auto& rhs
        ) {
            if (
                lhs.position !=
                rhs.position
            ) {
                return
                    lhs.position <
                    rhs.position;
            }

            return
                lhs.mismatches <
                rhs.mismatches;
        }
    );
}


void
compare_result(
    primerpair::AnchorCandidateSearchResult expected,
    primerpair::AnchorCandidateSearchResult observed,
    const std::size_t check_id
) {
    if (
        expected.primer_length !=
            observed.primer_length
        ||
        expected.anchor_length !=
            observed.anchor_length
        ||
        expected.max_mismatches !=
            observed.max_mismatches
        ||
        expected.anchor_occurrences !=
            observed.anchor_occurrences
        ||
        expected.candidates_verified !=
            observed.candidates_verified
    ) {
        std::cerr
            << "METADATA_MISMATCH"
            << '\t'
            << "check="
            << check_id
            << '\t'
            << "expected_occ="
            << expected.anchor_occurrences
            << '\t'
            << "observed_occ="
            << observed.anchor_occurrences
            << '\t'
            << "expected_verified="
            << expected.candidates_verified
            << '\t'
            << "observed_verified="
            << observed.candidates_verified
            << '\n';


        throw std::runtime_error(
            "Candidate metadata mismatch."
        );
    }


    normalize_hits(
        expected.hits
    );

    normalize_hits(
        observed.hits
    );


    if (
        expected.hits.size() !=
        observed.hits.size()
    ) {
        throw std::runtime_error(
            "Candidate hit-count mismatch."
        );
    }


    for (
        std::size_t i = 0;
        i < expected.hits.size();
        ++i
    ) {
        if (
            expected.hits[i].position !=
                observed.hits[i].position
            ||
            expected.hits[i].mismatches !=
                observed.hits[i].mismatches
        ) {
            std::cerr
                << "HIT_MISMATCH"
                << '\t'
                << "check="
                << check_id
                << '\t'
                << "hit="
                << i
                << '\t'
                << "expected_pos="
                << expected.hits[i].position
                << '\t'
                << "observed_pos="
                << observed.hits[i].position
                << '\t'
                << "expected_mm="
                << expected.hits[i].mismatches
                << '\t'
                << "observed_mm="
                << observed.hits[i].mismatches
                << '\n';


            throw std::runtime_error(
                "Candidate hit mismatch."
            );
        }
    }
}


char
mutated_base(
    const char base
) {
    return
        (
            base ==
            'A'
        )
        ?
        'T'
        :
        'A';
}

}  // namespace


int
main() {
    try {
        using namespace primerpair;


        constexpr std::size_t
            anchor_length = 12;


        constexpr std::size_t
            primer_length = 24;


        const std::string reference =
            make_reference(
                200000
            );


        PackedReference packed(
            reference
        );


        BidirectionalFMIndex bifm(
            reference
        );


        IPBWTIndex ipbwt(
            reference,
            anchor_length
        );


        AnchorCandidateSearcher searcher(
            bifm,
            packed
        );


        std::size_t checks = 0;


        for (
            std::size_t i = 0;
            i < 1500;
            ++i
        ) {
            const std::size_t position =
                (
                    i *
                    3571
                    +
                    123
                )
                %
                (
                    reference.size()
                    -
                    primer_length
                );


            const std::string original =
                reference.substr(
                    position,
                    primer_length
                );


            /*
             * Test four versions of each primer:
             *
             * 0 mismatches in prefix
             * 1 mismatch in prefix
             * 2 mismatches in prefix
             * 3 mismatches in prefix
             *
             * The exact 12-mer 3' anchor remains
             * untouched.
             */
            for (
                std::size_t introduced = 0;
                introduced <= 3;
                ++introduced
            ) {
                std::string primer =
                    original;


                for (
                    std::size_t m = 0;
                    m < introduced;
                    ++m
                ) {
                    const std::size_t p =
                        1 +
                        m *
                        3;


                    primer[p] =
                        mutated_base(
                            primer[p]
                        );
                }


                const auto expected =
                    searcher.search(
                        primer,
                        anchor_length,
                        3
                    );


                const std::string_view anchor(
                    primer.data()
                        +
                        primer.size()
                        -
                        anchor_length,
                    anchor_length
                );


                const Interval interval =
                    ipbwt.exact_search(
                        anchor
                    );


                auto positions =
                    ipbwt.locate(
                        interval
                    );


                const auto observed =
                    searcher
                        .verify_from_anchor_positions(
                            primer,
                            std::move(
                                positions
                            ),
                            interval.size(),
                            anchor_length,
                            3
                        );


                compare_result(
                    expected,
                    observed,
                    checks
                );


                ++checks;
            }
        }


        /*
         * Also exercise an exact anchor that is
         * deliberately changed and may therefore
         * become absent.
         */
        for (
            std::size_t i = 0;
            i < 500;
            ++i
        ) {
            const std::size_t position =
                (
                    i *
                    7919
                    +
                    999
                )
                %
                (
                    reference.size()
                    -
                    primer_length
                );


            std::string primer =
                reference.substr(
                    position,
                    primer_length
                );


            primer[
                primer_length -
                3
            ] =
                mutated_base(
                    primer[
                        primer_length -
                        3
                    ]
                );


            const auto expected =
                searcher.search(
                    primer,
                    anchor_length,
                    3
                );


            const std::string_view anchor(
                primer.data()
                    +
                    primer.size()
                    -
                    anchor_length,
                anchor_length
            );


            const Interval interval =
                ipbwt.exact_search(
                    anchor
                );


            auto positions =
                ipbwt.locate(
                    interval
                );


            const auto observed =
                searcher
                    .verify_from_anchor_positions(
                        primer,
                        std::move(
                            positions
                        ),
                        interval.size(),
                        anchor_length,
                        3
                    );


            compare_result(
                expected,
                observed,
                checks
            );


            ++checks;
        }


        /*
         * Explicit repeat-rich primer.
         */
        {
            const std::string primer =
                "TTTTTTTTTTTT"
                "ACGTACGTACGT";


            const auto expected =
                searcher.search(
                    primer,
                    anchor_length,
                    3
                );


            const std::string_view anchor(
                primer.data()
                    +
                    primer.size()
                    -
                    anchor_length,
                anchor_length
            );


            const Interval interval =
                ipbwt.exact_search(
                    anchor
                );


            auto positions =
                ipbwt.locate(
                    interval
                );


            const auto observed =
                searcher
                    .verify_from_anchor_positions(
                        primer,
                        std::move(
                            positions
                        ),
                        interval.size(),
                        anchor_length,
                        3
                    );


            compare_result(
                expected,
                observed,
                checks
            );


            ++checks;
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
