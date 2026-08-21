#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/search_strategy.hpp"

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


char
complement(
    const char base
) {
    switch (base) {

        case 'A':
            return 'T';

        case 'C':
            return 'G';

        case 'G':
            return 'C';

        case 'T':
            return 'A';
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
            complement(
                *it
            )
        );
    }


    return result;
}


std::string_view
expected_anchor(
    const std::string& query,
    const primerpair::AnchorPlacement placement,
    const std::size_t anchor_length
) {
    if (
        placement ==
        primerpair::AnchorPlacement::Prefix
    ) {
        return std::string_view(
            query.data(),
            anchor_length
        );
    }


    return std::string_view(
        query.data()
            +
            query.size()
            -
            anchor_length,
        anchor_length
    );
}

}  // namespace


int
main() {
    try {
        using namespace primerpair;


        constexpr std::size_t
            anchor_length = 12;


        const std::string reference =
            make_reference(
                200000
            );


        BidirectionalFMIndex bifm(
            reference
        );


        IPBWTIndex ipbwt(
            reference,
            anchor_length
        );


        SearchDifficultyEstimator estimator(
            bifm
        );


        BatchedAnchorLookup lookup(
            ipbwt,
            estimator
        );


        std::vector<std::string>
            owned_queries;


        std::vector<BatchedAnchorRequest>
            requests;


        owned_queries.reserve(
            1024
        );

        requests.reserve(
            1024
        );


        /*
         * Mix forward suffix anchors and the
         * reverse-query prefix semantics used by
         * StrandAwarePrimerSearchEngine.
         */
        for (
            std::size_t i = 0;
            i < 512;
            ++i
        ) {
            constexpr std::size_t
                primer_length = 24;


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


            const std::string primer =
                reference.substr(
                    position,
                    primer_length
                );


            owned_queries.push_back(
                primer
            );


            requests.push_back(
                BatchedAnchorRequest{
                    owned_queries.back(),
                    AnchorPlacement::Suffix,
                    3
                }
            );


            owned_queries.push_back(
                reverse_complement(
                    primer
                )
            );


            requests.push_back(
                BatchedAnchorRequest{
                    owned_queries.back(),
                    AnchorPlacement::Prefix,
                    3
                }
            );
        }


        /*
         * Rebuild request string_views after all
         * strings are inserted, preventing vector
         * relocation from invalidating views.
         */
        requests.clear();

        requests.reserve(
            owned_queries.size()
        );


        for (
            std::size_t i = 0;
            i < owned_queries.size();
            ++i
        ) {
            requests.push_back(
                BatchedAnchorRequest{
                    owned_queries[i],
                    (
                        i %
                        2 ==
                        0
                    )
                    ?
                    AnchorPlacement::Suffix
                    :
                    AnchorPlacement::Prefix,
                    3
                }
            );
        }


        const auto observed =
            lookup.lookup(
                requests,
                anchor_length
            );


        if (
            observed.size() !=
            requests.size()
        ) {
            throw std::runtime_error(
                "Decision count mismatch."
            );
        }


        std::size_t checks = 0;


        for (
            std::size_t i = 0;
            i < requests.size();
            ++i
        ) {
            const auto anchor =
                expected_anchor(
                    owned_queries[i],
                    requests[i].placement,
                    anchor_length
                );


            const auto bifm_interval =
                bifm.search(
                    anchor
                );


            const std::uint64_t
                expected_occurrences =
                    bifm_interval.size();


            const auto
                expected_difficulty =
                    estimator.classify(
                        expected_occurrences
                    );


            const auto
                expected_strategy =
                    SearchStrategyRouter::choose(
                        expected_difficulty,
                        requests[i]
                            .max_mismatches
                    );


            if (
                observed[i].occurrences !=
                    expected_occurrences
            ) {
                throw std::runtime_error(
                    "Occurrence-count mismatch."
                );
            }


            if (
                observed[i].difficulty !=
                    expected_difficulty
            ) {
                throw std::runtime_error(
                    "Difficulty-routing mismatch."
                );
            }


            if (
                observed[i].strategy !=
                    expected_strategy
            ) {
                throw std::runtime_error(
                    "Strategy-routing mismatch."
                );
            }


            auto bifm_positions =
                bifm.locate(
                    bifm_interval
                );


            auto ip_positions =
                ipbwt.locate(
                    observed[i].interval
                );


            std::sort(
                bifm_positions.begin(),
                bifm_positions.end()
            );

            std::sort(
                ip_positions.begin(),
                ip_positions.end()
            );


            bifm_positions.erase(
                std::unique(
                    bifm_positions.begin(),
                    bifm_positions.end()
                ),
                bifm_positions.end()
            );


            ip_positions.erase(
                std::unique(
                    ip_positions.begin(),
                    ip_positions.end()
                ),
                ip_positions.end()
            );


            if (
                bifm_positions !=
                    ip_positions
            ) {
                throw std::runtime_error(
                    "Anchor-position mismatch."
                );
            }


            ++checks;
        }


        /*
         * Empty batch.
         */
        {
            const std::vector<
                BatchedAnchorRequest
            > empty;


            const auto result =
                lookup.lookup(
                    empty,
                    anchor_length
                );


            if (!result.empty()) {
                throw std::runtime_error(
                    "Empty batch produced results."
                );
            }


            ++checks;
        }


        /*
         * Invalid anchor length.
         */
        {
            bool threw = false;


            try {
                static_cast<void>(
                    lookup.lookup(
                        requests,
                        0
                    )
                );

            } catch (
                const std::invalid_argument&
            ) {
                threw = true;
            }


            if (!threw) {
                throw std::runtime_error(
                    "Zero anchor length "
                    "was not rejected."
                );
            }


            ++checks;
        }


        std::cout
            << "anchor_length\t"
            << anchor_length
            << '\n';


        std::cout
            << "requests\t"
            << requests.size()
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
