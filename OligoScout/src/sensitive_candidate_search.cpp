#include "primerpair/sensitive_candidate_search.hpp"

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
    kBases{
        'A',
        'C',
        'G',
        'T'
    };

struct SeedBranch {
    BidirectionalInterval interval{};
    std::size_t mismatches{0};
};

std::string normalize(
    const std::string_view primer
) {
    if (
        primer.size() < 2 ||
        primer.size() > 64
    ) {
        throw std::invalid_argument(
            "SENSITIVE candidate primer length "
            "must be in the range 2..64."
        );
    }

    std::string output;

    output.reserve(
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

        if (
            base != 'A' &&
            base != 'C' &&
            base != 'G' &&
            base != 'T'
        ) {
            throw std::invalid_argument(
                "Primer must contain only A/C/G/T."
            );
        }

        output.push_back(
            base
        );
    }

    return output;
}

bool hit_less(
    const OrientedPrimerSearchHit& lhs,
    const OrientedPrimerSearchHit& rhs
) noexcept {
    if (lhs.position != rhs.position) {
        return lhs.position < rhs.position;
    }

    if (lhs.orientation != rhs.orientation) {
        return
            static_cast<int>(lhs.orientation) <
            static_cast<int>(rhs.orientation);
    }

    if (lhs.mismatches != rhs.mismatches) {
        return lhs.mismatches < rhs.mismatches;
    }

    return
        lhs.mismatch_mask <
        rhs.mismatch_mask;
}

}  // namespace

SensitiveCandidateSearchEngine::
SensitiveCandidateSearchEngine(
    const BidirectionalFMIndex& index,
    const PackedReference& reference
)
    : index_(index),
      reference_(reference) {
}

std::vector<std::uint64_t>
SensitiveCandidateSearchEngine::search_seed(
    const std::string_view seed,
    const std::size_t max_mismatches
) const {
    if (seed.empty()) {
        throw std::invalid_argument(
            "Seed cannot be empty."
        );
    }

    /*
     * Exact seed fast path.
     */
    if (max_mismatches == 0) {

        const auto interval =
            index_.search(
                seed
            );

        if (interval.empty()) {
            return {};
        }

        auto positions =
            index_.locate(
                interval
            );

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

        return positions;
    }

    /*
     * Current optimized backend only needs
     * seed budget 0 or 1 for k <= 3.
     */
    if (max_mismatches > 1) {
        throw std::invalid_argument(
            "Half-seed backend currently supports "
            "seed mismatch budget <= 1."
        );
    }

    const std::size_t last =
        seed.size() - 1;

    std::vector<SeedBranch>
        branches;

    for (const char genomic : kBases) {

        const std::size_t mismatches =
            genomic == seed.at(last)
                ? 0
                : 1;

        if (
            mismatches >
            max_mismatches
        ) {
            continue;
        }

        const std::string one_base(
            1,
            genomic
        );

        const auto interval =
            index_.search(
                one_base
            );

        if (!interval.empty()) {
            branches.push_back(
                SeedBranch{
                    interval,
                    mismatches
                }
            );
        }
    }

    for (
        std::size_t remaining = last;
        remaining > 0;
        --remaining
    ) {
        const std::size_t position =
            remaining - 1;

        std::vector<SeedBranch>
            next;

        for (const auto& branch : branches) {

            for (const char genomic : kBases) {

                const std::size_t mismatches =
                    branch.mismatches +
                    (
                        genomic ==
                        seed.at(position)
                            ? 0
                            : 1
                    );

                if (
                    mismatches >
                    max_mismatches
                ) {
                    continue;
                }

                const auto interval =
                    index_.extend_left(
                        branch.interval,
                        genomic
                    );

                if (!interval.empty()) {
                    next.push_back(
                        SeedBranch{
                            interval,
                            mismatches
                        }
                    );
                }
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

    std::vector<std::uint64_t>
        positions;

    for (const auto& branch : branches) {

        const auto located =
            index_.locate(
                branch.interval
            );

        positions.insert(
            positions.end(),
            located.begin(),
            located.end()
        );
    }

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

    return positions;
}

bool SensitiveCandidateSearchEngine::
verify_candidate(
    const std::uint64_t start,
    const std::string_view query,
    const bool reverse_to_original,
    const std::size_t max_mismatches,
    std::size_t& mismatches,
    std::uint64_t& mismatch_mask
) const {
    mismatches = 0;
    mismatch_mask = 0;

    if (
        start >
        reference_.size()
    ) {
        return false;
    }

    if (
        query.size() >
        reference_.size() -
        start
    ) {
        return false;
    }

    for (
        std::size_t i = 0;
        i < query.size();
        ++i
    ) {
        const char reference_base =
            reference_.base_at(
                start +
                static_cast<std::uint64_t>(
                    i
                )
            );

        /*
         * PackedReference deliberately rejects
         * ambiguous reference sequence.
         */
        if (
            reference_base != 'A' &&
            reference_base != 'C' &&
            reference_base != 'G' &&
            reference_base != 'T'
        ) {
            return false;
        }

        if (
            reference_base ==
            query.at(i)
        ) {
            continue;
        }

        ++mismatches;

        if (
            mismatches >
            max_mismatches
        ) {
            return false;
        }

        const std::size_t original_position =
            reverse_to_original
                ? query.size() - 1 - i
                : i;

        mismatch_mask |=
            (
                std::uint64_t{1}
                <<
                original_position
            );
    }

    return true;
}

std::vector<OrientedPrimerSearchHit>
SensitiveCandidateSearchEngine::
search_oriented(
    const std::string_view query,
    const PrimerOrientation orientation,
    const bool reverse_to_original,
    const std::size_t max_mismatches
) const {
    /*
     * k=0: full exact search is cheaper than
     * generating half-seed candidates.
     */
    if (max_mismatches == 0) {

        const auto interval =
            index_.search(
                query
            );

        std::vector<
            OrientedPrimerSearchHit
        > exact_hits;

        if (interval.empty()) {
            return exact_hits;
        }

        const auto positions =
            index_.locate(
                interval
            );

        exact_hits.reserve(
            positions.size()
        );

        for (
            const std::uint64_t position :
            positions
        ) {
            exact_hits.push_back(
                OrientedPrimerSearchHit{
                    position,
                    0,
                    orientation,
                    0
                }
            );
        }

        std::sort(
            exact_hits.begin(),
            exact_hits.end(),
            hit_less
        );

        exact_hits.erase(
            std::unique(
                exact_hits.begin(),
                exact_hits.end()
            ),
            exact_hits.end()
        );

        return exact_hits;
    }

    const std::size_t split =
        query.size() / 2;

    const std::string_view left_seed =
        query.substr(
            0,
            split
        );

    const std::string_view right_seed =
        query.substr(
            split
        );

    /*
     * Lossless bound:
     *
     * d_left + d_right <= k
     *
     * therefore at least one half has
     *
     * d_half <= floor(k / 2).
     */
    const std::size_t seed_budget =
        max_mismatches / 2;

    const auto left_positions =
        search_seed(
            left_seed,
            seed_budget
        );

    const auto right_positions =
        search_seed(
            right_seed,
            seed_budget
        );

    std::vector<std::uint64_t>
        candidate_starts;

    candidate_starts.reserve(
        left_positions.size() +
        right_positions.size()
    );

    /*
     * Left seed begins at full-primer offset 0.
     */
    for (
        const std::uint64_t position :
        left_positions
    ) {
        candidate_starts.push_back(
            position
        );
    }

    /*
     * Right seed begins at full-primer offset split.
     */
    for (
        const std::uint64_t position :
        right_positions
    ) {
        if (
            position <
            split
        ) {
            continue;
        }

        candidate_starts.push_back(
            position -
            static_cast<std::uint64_t>(
                split
            )
        );
    }

    std::sort(
        candidate_starts.begin(),
        candidate_starts.end()
    );

    candidate_starts.erase(
        std::unique(
            candidate_starts.begin(),
            candidate_starts.end()
        ),
        candidate_starts.end()
    );

    std::vector<
        OrientedPrimerSearchHit
    > hits;

    for (
        const std::uint64_t start :
        candidate_starts
    ) {
        std::size_t mismatches = 0;
        std::uint64_t mismatch_mask = 0;

        if (
            !verify_candidate(
                start,
                query,
                reverse_to_original,
                max_mismatches,
                mismatches,
                mismatch_mask
            )
        ) {
            continue;
        }

        hits.push_back(
            OrientedPrimerSearchHit{
                start,
                mismatches,
                orientation,
                mismatch_mask
            }
        );
    }

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

    return hits;
}

SensitivePrimerSearchResult
SensitiveCandidateSearchEngine::search(
    const std::string_view primer,
    const std::size_t max_mismatches
) const {
    if (max_mismatches > 3) {
        throw std::invalid_argument(
            "SENSITIVE candidate MVP supports "
            "at most 3 mismatches."
        );
    }

    const std::string normalized =
        normalize(
            primer
        );

    auto forward =
        search_oriented(
            normalized,
            PrimerOrientation::Forward,
            false,
            max_mismatches
        );

    const std::string reverse_query =
        reverse_complement(
            normalized
        );

    auto reverse =
        search_oriented(
            reverse_query,
            PrimerOrientation::Reverse,
            true,
            max_mismatches
        );

    std::vector<
        OrientedPrimerSearchHit
    > hits;

    hits.reserve(
        forward.size() +
        reverse.size()
    );

    hits.insert(
        hits.end(),
        forward.begin(),
        forward.end()
    );

    hits.insert(
        hits.end(),
        reverse.begin(),
        reverse.end()
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
