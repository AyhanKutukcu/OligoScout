#include "primerpair/sensitive_primer_search.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace primerpair {

namespace {

constexpr std::array<char, 4>
    kCanonicalBases{
        'A',
        'C',
        'G',
        'T'
    };

std::string normalize_primer(
    const std::string_view primer
) {
    if (
        primer.empty() ||
        primer.size() > 64
    ) {
        throw std::invalid_argument(
            "SENSITIVE primer length must be in "
            "the range 1..64."
        );
    }

    std::string normalized;

    normalized.reserve(
        primer.size()
    );

    for (const char raw : primer) {

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
                    "SENSITIVE primer must contain "
                    "only A/C/G/T."
                );
        }
    }

    return normalized;
}

std::size_t original_position(
    const std::size_t query_position,
    const std::size_t query_length,
    const bool reverse_to_original
) {
    return
        reverse_to_original
            ? query_length -
                1 -
                query_position
            : query_position;
}

bool hit_less(
    const OrientedPrimerSearchHit& lhs,
    const OrientedPrimerSearchHit& rhs
) noexcept {
    if (
        lhs.position !=
        rhs.position
    ) {
        return
            lhs.position <
            rhs.position;
    }

    if (
        lhs.orientation !=
        rhs.orientation
    ) {
        return
            static_cast<int>(
                lhs.orientation
            ) <
            static_cast<int>(
                rhs.orientation
            );
    }

    if (
        lhs.mismatches !=
        rhs.mismatches
    ) {
        return
            lhs.mismatches <
            rhs.mismatches;
    }

    return
        lhs.mismatch_mask <
        rhs.mismatch_mask;
}

}  // namespace

SensitivePrimerSearchEngine::
SensitivePrimerSearchEngine(
    const BidirectionalFMIndex& index
)
    : index_(
          index
      ) {
}

std::vector<OrientedPrimerSearchHit>
SensitivePrimerSearchEngine::
search_reference_oriented_query(
    const std::string_view query,
    const PrimerOrientation orientation,
    const bool reverse_to_original,
    const std::size_t max_mismatches
) const {
    if (query.empty()) {
        throw std::invalid_argument(
            "SENSITIVE query cannot be empty."
        );
    }

    /*
     * --------------------------------------------------------
     * Search direction
     *
     * BiFM extension begins with the RIGHT-most query base
     * and repeatedly extend_left().
     *
     * Unlike STRICT, the initial base is NOT assumed exact.
     * We branch A/C/G/T immediately.
     * --------------------------------------------------------
     */

    const std::size_t
        last_query_position =
            query.size() - 1;

    std::vector<BranchState>
        branches;

    branches.reserve(
        4
    );

    /*
     * Initial one-base states.
     */
    for (
        const char genomic_base :
        kCanonicalBases
    ) {
        const bool mismatch =
            genomic_base !=
            query.at(
                last_query_position
            );

        const std::size_t mismatch_count =
            mismatch
                ? 1
                : 0;

        if (
            mismatch_count >
            max_mismatches
        ) {
            continue;
        }

        const std::string seed(
            1,
            genomic_base
        );

        const auto interval =
            index_.search(
                seed
            );

        if (interval.empty()) {
            continue;
        }

        std::uint64_t mask = 0;

        if (mismatch) {

            const std::size_t position =
                original_position(
                    last_query_position,
                    query.size(),
                    reverse_to_original
                );

            mask |=
                (
                    std::uint64_t{1}
                    <<
                    position
                );
        }

        branches.push_back(
            BranchState{
                interval,
                mismatch_count,
                mask
            }
        );
    }

    /*
     * Extend through all remaining primer positions,
     * from right to left.
     */
    for (
        std::size_t remaining =
            last_query_position;
        remaining > 0;
        --remaining
    ) {
        const std::size_t query_position =
            remaining - 1;

        std::vector<BranchState>
            next;

        /*
         * Four possible genomic symbols can be
         * explored for every surviving state.
         */
        next.reserve(
            branches.size() * 2
        );

        for (const auto& branch : branches) {

            for (
                const char genomic_base :
                kCanonicalBases
            ) {
                const bool mismatch =
                    genomic_base !=
                    query.at(
                        query_position
                    );

                const std::size_t
                    mismatch_count =
                        branch.mismatches +
                        (
                            mismatch
                                ? 1
                                : 0
                        );

                if (
                    mismatch_count >
                    max_mismatches
                ) {
                    continue;
                }

                const auto interval =
                    index_.extend_left(
                        branch.interval,
                        genomic_base
                    );

                if (interval.empty()) {
                    continue;
                }

                std::uint64_t mask =
                    branch.mismatch_mask;

                if (mismatch) {

                    const std::size_t
                        position =
                            original_position(
                                query_position,
                                query.size(),
                                reverse_to_original
                            );

                    mask |=
                        (
                            std::uint64_t{1}
                            <<
                            position
                        );
                }

                next.push_back(
                    BranchState{
                        interval,
                        mismatch_count,
                        mask
                    }
                );
            }
        }

        branches =
            std::move(
                next
            );

        if (branches.empty()) {
            break;
        }
    }

    std::vector<
        OrientedPrimerSearchHit
    > hits;

    /*
     * Every completed branch describes a concrete
     * full-length genomic sequence within the Hamming
     * budget.
     */
    for (const auto& branch : branches) {

        const auto positions =
            index_.locate(
                branch.interval
            );

        for (
            const std::uint64_t position :
            positions
        ) {
            hits.push_back(
                OrientedPrimerSearchHit{
                    position,
                    branch.mismatches,
                    orientation,
                    branch.mismatch_mask
                }
            );
        }
    }

    std::sort(
        hits.begin(),
        hits.end(),
        hit_less
    );

    /*
     * Defensive duplicate elimination.
     */
    hits.erase(
        std::unique(
            hits.begin(),
            hits.end()
        ),
        hits.end()
    );

    return hits;
}

SensitivePrimerSearchResult
SensitivePrimerSearchEngine::search(
    const std::string_view primer,
    const std::size_t max_mismatches
) const {
    /*
     * Keep the current MVP budget bounded.
     */
    if (max_mismatches > 3) {
        throw std::invalid_argument(
            "SENSITIVE MVP supports at most "
            "3 mismatches."
        );
    }

    const std::string normalized =
        normalize_primer(
            primer
        );

    /*
     * Forward genomic orientation.
     */
    auto forward_hits =
        search_reference_oriented_query(
            normalized,
            PrimerOrientation::Forward,
            false,
            max_mismatches
        );

    /*
     * Reverse genomic orientation:
     *
     * Search the reverse-complement sequence against
     * the reference, but map mismatch positions back
     * into ORIGINAL primer coordinates.
     */
    const std::string reverse_query =
        reverse_complement(
            normalized
        );

    auto reverse_hits =
        search_reference_oriented_query(
            reverse_query,
            PrimerOrientation::Reverse,
            true,
            max_mismatches
        );

    std::vector<
        OrientedPrimerSearchHit
    > hits;

    hits.reserve(
        forward_hits.size() +
        reverse_hits.size()
    );

    hits.insert(
        hits.end(),
        forward_hits.begin(),
        forward_hits.end()
    );

    hits.insert(
        hits.end(),
        reverse_hits.begin(),
        reverse_hits.end()
    );

    std::sort(
        hits.begin(),
        hits.end(),
        hit_less
    );

    hits.erase(
        std::unique(
            hits.begin(),
            hits.end()
        ),
        hits.end()
    );

    return SensitivePrimerSearchResult{
        normalized.size(),
        max_mismatches,
        std::move(
            hits
        )
    };
}

}  // namespace primerpair
