#include "primerpair/batched_anchor_lookup.hpp"

#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace primerpair {

namespace {

std::string
normalize_query(
    const std::string_view query
) {
    if (query.empty()) {
        throw std::invalid_argument(
            "Anchor query cannot be empty."
        );
    }

    std::string normalized;

    normalized.reserve(
        query.size()
    );


    for (
        const char raw :
        query
    ) {
        const char base =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        raw
                    )
                )
            );


        switch (base) {

            case 'A':
            case 'C':
            case 'G':
            case 'T':

                normalized.push_back(
                    base
                );

                break;

            default:

                throw std::invalid_argument(
                    "Anchor query must contain "
                    "only A/C/G/T."
                );
        }
    }


    return normalized;
}

}  // namespace


BatchedAnchorLookup::
BatchedAnchorLookup(
    const IPBWTIndex& index,
    const SearchDifficultyEstimator& estimator
) noexcept
    : index_(index),
      estimator_(estimator) {
}


std::vector<BatchedAnchorDecision>
BatchedAnchorLookup::lookup(
    const std::vector<BatchedAnchorRequest>& requests,
    const std::size_t anchor_length
) const {
    if (
        anchor_length == 0
    ) {
        throw std::invalid_argument(
            "Anchor length must be > 0."
        );
    }


    /*
     * Own normalized strings so the string_views
     * passed into the batch API remain valid.
     */
    std::vector<std::string>
        normalized_queries;

    normalized_queries.reserve(
        requests.size()
    );


    for (
        const auto& request :
        requests
    ) {
        normalized_queries.push_back(
            normalize_query(
                request.query
            )
        );


        if (
            anchor_length >
            normalized_queries.back().size()
        ) {
            throw std::invalid_argument(
                "Anchor length cannot exceed "
                "query length."
            );
        }


        if (
            request.max_mismatches >
            3
        ) {
            throw std::invalid_argument(
                "BatchedAnchorLookup currently "
                "supports at most 3 mismatches."
            );
        }
    }


    std::vector<std::string_view>
        anchors;

    anchors.reserve(
        requests.size()
    );


    for (
        std::size_t i = 0;
        i < requests.size();
        ++i
    ) {
        const auto& query =
            normalized_queries[i];


        switch (
            requests[i].placement
        ) {

            case AnchorPlacement::Suffix:

                anchors.emplace_back(
                    query.data()
                        +
                        query.size()
                        -
                        anchor_length,
                    anchor_length
                );

                break;


            case AnchorPlacement::Prefix:

                anchors.emplace_back(
                    query.data(),
                    anchor_length
                );

                break;
        }
    }


    const auto intervals =
        index_.exact_prefix_search_many(
            anchors
        );


    std::vector<BatchedAnchorDecision>
        decisions;

    decisions.reserve(
        requests.size()
    );


    for (
        std::size_t i = 0;
        i < requests.size();
        ++i
    ) {
        const std::uint64_t occurrences =
            intervals[i].size();


        const SearchDifficulty difficulty =
            estimator_.classify(
                occurrences
            );


        const SearchStrategy strategy =
            SearchStrategyRouter::choose(
                difficulty,
                requests[i].max_mismatches
            );


        decisions.push_back(
            BatchedAnchorDecision{
                intervals[i],
                occurrences,
                difficulty,
                strategy,
                normalized_queries[i].size(),
                anchor_length
            }
        );
    }


    return decisions;
}

}  // namespace primerpair
