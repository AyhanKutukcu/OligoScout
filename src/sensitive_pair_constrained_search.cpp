#include "primerpair/sensitive_pair_constrained_search.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace primerpair {

/*
 * Narrow internal bridge used by the pair-constrained
 * engine. SensitiveCandidateSearchEngine keeps its
 * validated primitives private to normal callers.
 */
struct SensitiveCandidateSearchAccess {

    [[nodiscard]]
    static std::vector<std::uint64_t>
    search_seed(
        const SensitiveCandidateSearchEngine& engine,
        const std::string_view seed,
        const std::size_t max_mismatches
    ) {
        return engine.search_seed(
            seed,
            max_mismatches
        );
    }

    [[nodiscard]]
    static bool verify_candidate(
        const SensitiveCandidateSearchEngine& engine,
        const std::uint64_t start,
        const std::string_view query,
        const bool reverse_to_original,
        const std::size_t max_mismatches,
        std::size_t& mismatches,
        std::uint64_t& mismatch_mask
    ) {
        return engine.verify_candidate(
            start,
            query,
            reverse_to_original,
            max_mismatches,
            mismatches,
            mismatch_mask
        );
    }
};

namespace {

double elapsed_us(
    const std::chrono::steady_clock::time_point start,
    const std::chrono::steady_clock::time_point stop
) {
    return static_cast<double>(
        std::chrono::duration_cast<
            std::chrono::nanoseconds
        >(
            stop - start
        ).count()
    ) / 1000.0;
}


struct StartWindow {
    std::uint64_t begin{0};
    std::uint64_t end{0};  // inclusive
};


bool window_less(
    const StartWindow& lhs,
    const StartWindow& rhs
) noexcept {
    if (lhs.begin != rhs.begin) {
        return lhs.begin < rhs.begin;
    }

    return lhs.end < rhs.end;
}


std::uint64_t saturating_add(
    const std::uint64_t lhs,
    const std::uint64_t rhs
) noexcept {
    if (
        lhs >
        std::numeric_limits<
            std::uint64_t
        >::max() -
        rhs
    ) {
        return
            std::numeric_limits<
                std::uint64_t
            >::max();
    }

    return lhs + rhs;
}


std::uint64_t saturating_sub(
    const std::uint64_t lhs,
    const std::uint64_t rhs
) noexcept {
    return
        lhs >= rhs
            ? lhs - rhs
            : 0;
}


std::string normalize_primer(
    const std::string_view primer
) {
    /*
     * v1 threshold has only been validated on
     * 20-mers.
     */
    if (
        primer.size() < 18
        ||
        primer.size() > 35
    ) {
        throw std::invalid_argument(
            "Pair-constrained SENSITIVE MVP supports "
            "primer lengths from 18 through 35 nt."
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
            static_cast<int>(
                lhs.orientation
            )
            <
            static_cast<int>(
                rhs.orientation
            );
    }

    if (lhs.mismatches != rhs.mismatches) {
        return
            lhs.mismatches <
            rhs.mismatches;
    }

    return
        lhs.mismatch_mask <
        rhs.mismatch_mask;
}


void normalize_hits(
    std::vector<
        OrientedPrimerSearchHit
    >& hits
) {
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
}


std::vector<StartWindow>
merge_windows(
    std::vector<StartWindow> windows
) {
    if (windows.empty()) {
        return {};
    }

    std::sort(
        windows.begin(),
        windows.end(),
        window_less
    );

    std::vector<StartWindow>
        merged;

    merged.reserve(
        windows.size()
    );

    merged.push_back(
        windows.front()
    );

    for (
        std::size_t i = 1;
        i < windows.size();
        ++i
    ) {
        auto& current =
            merged.back();

        const auto& next =
            windows.at(i);

        const std::uint64_t
            current_plus_one =
                current.end ==
                std::numeric_limits<
                    std::uint64_t
                >::max()
                    ? current.end
                    : current.end + 1;

        if (
            next.begin <=
            current_plus_one
        ) {
            current.end =
                std::max(
                    current.end,
                    next.end
                );

        } else {

            merged.push_back(
                next
            );
        }
    }

    return merged;
}


/*
 * Anchor FORWARD:
 *
 * ANCHOR --->                <--- MATE
 *
 * Return allowed REVERSE-mate start positions.
 */
void add_downstream_window(
    const OrientedPrimerSearchHit& anchor,
    const std::size_t anchor_length,
    const std::size_t mate_length,
    const std::uint64_t reference_size,
    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon,
    std::vector<StartWindow>& output
) {
    if (
        max_amplicon <
        mate_length
    ) {
        return;
    }

    if (
        reference_size <
        mate_length
    ) {
        return;
    }

    const std::uint64_t
        max_reference_start =
            reference_size -
            static_cast<std::uint64_t>(
                mate_length
            );

    std::uint64_t lower =
        saturating_add(
            anchor.position,
            static_cast<std::uint64_t>(
                anchor_length
            )
        );

    if (
        min_amplicon >
        mate_length
    ) {
        const auto by_min =
            saturating_add(
                anchor.position,
                min_amplicon -
                static_cast<std::uint64_t>(
                    mate_length
                )
            );

        lower =
            std::max(
                lower,
                by_min
            );
    }

    std::uint64_t upper =
        saturating_add(
            anchor.position,
            max_amplicon -
            static_cast<std::uint64_t>(
                mate_length
            )
        );

    upper =
        std::min(
            upper,
            max_reference_start
        );

    if (lower > upper) {
        return;
    }

    output.push_back(
        StartWindow{
            lower,
            upper
        }
    );
}


/*
 * Anchor REVERSE:
 *
 * MATE --->                 <--- ANCHOR
 *
 * Return allowed FORWARD-mate start positions.
 */
void add_upstream_window(
    const OrientedPrimerSearchHit& anchor,
    const std::size_t anchor_length,
    const std::size_t mate_length,
    const std::uint64_t reference_size,
    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon,
    std::vector<StartWindow>& output
) {
    if (
        reference_size <
        mate_length
    ) {
        return;
    }

    /*
     * Non-overlap requires:
     *
     * mate_start + mate_length
     * <= anchor.position
     */
    if (
        anchor.position <
        mate_length
    ) {
        return;
    }

    const std::uint64_t
        anchor_end =
            saturating_add(
                anchor.position,
                static_cast<std::uint64_t>(
                    anchor_length
                )
            );

    if (
        anchor_end <
        min_amplicon
    ) {
        return;
    }

    std::uint64_t lower =
        saturating_sub(
            anchor_end,
            max_amplicon
        );

    std::uint64_t upper =
        anchor_end -
        min_amplicon;

    const std::uint64_t
        by_nonoverlap =
            anchor.position -
            static_cast<std::uint64_t>(
                mate_length
            );

    upper =
        std::min(
            upper,
            by_nonoverlap
        );

    const std::uint64_t
        max_reference_start =
            reference_size -
            static_cast<std::uint64_t>(
                mate_length
            );

    upper =
        std::min(
            upper,
            max_reference_start
        );

    if (lower > upper) {
        return;
    }

    output.push_back(
        StartWindow{
            lower,
            upper
        }
    );
}


std::vector<std::uint64_t>
build_half_seed_candidate_starts(
    const SensitiveCandidateSearchEngine& searcher,
    const std::string_view query,
    const std::size_t max_mismatches
) {
    if (query.size() < 2) {
        throw std::invalid_argument(
            "Candidate query must contain "
            "at least two bases."
        );
    }

    /*
     * Pigeonhole theorem:
     *
     * d(left) + d(right) <= k
     *
     * therefore at least one half satisfies
     *
     * d(half) <= floor(k / 2).
     *
     * For the validated 20-mer / k=3 backend:
     *
     * 10-mer left  <= 1 mismatch
     * OR
     * 10-mer right <= 1 mismatch.
     */
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

    const std::size_t seed_budget =
        max_mismatches / 2;

    const auto left_positions =
        SensitiveCandidateSearchAccess::search_seed(
            searcher,
            left_seed,
            seed_budget
        );

    const auto right_positions =
        SensitiveCandidateSearchAccess::search_seed(
            searcher,
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
     * Left half starts at full-primer offset 0.
     */
    candidate_starts.insert(
        candidate_starts.end(),
        left_positions.begin(),
        left_positions.end()
    );

    /*
     * Right half begins at full-primer offset split.
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

    return candidate_starts;
}


std::vector<OrientedPrimerSearchHit>
search_candidates_in_windows(
    const SensitiveCandidateSearchEngine& searcher,
    const std::vector<StartWindow>& windows,
    const std::string_view query,
    const PrimerOrientation orientation,
    const bool reverse_to_original,
    const std::size_t max_mismatches,
    std::uint64_t& verified_candidate_positions,
    SensitivePairConstrainedTiming& timing
) {
    if (windows.empty()) {
        return {};
    }

    /*
     * Generate lossless full-primer candidate starts
     * with the already validated half-seed backend.
     */
    const auto seed_start =
        std::chrono::steady_clock::now();

    const auto candidate_starts =
        build_half_seed_candidate_starts(
            searcher,
            query,
            max_mismatches
        );

    timing.mate_seed_generation_us +=
        elapsed_us(
            seed_start,
            std::chrono::steady_clock::now()
        );

    std::vector<
        OrientedPrimerSearchHit
    > hits;

    std::size_t window_index = 0;

    /*
     * Both vectors are sorted.
     *
     * Therefore PCR-window filtering is O(C + W),
     * rather than one binary search per candidate.
     */
    const auto filter_verify_start =
        std::chrono::steady_clock::now();

    for (
        const std::uint64_t start :
        candidate_starts
    ) {
        while (
            window_index <
                windows.size()
            &&
            windows.at(
                window_index
            ).end <
                start
        ) {
            ++window_index;
        }

        if (
            window_index ==
            windows.size()
        ) {
            break;
        }

        const auto& window =
            windows.at(
                window_index
            );

        if (
            start <
            window.begin
        ) {
            continue;
        }

        /*
         * Only starts already satisfying the pair
         * constraint reach expensive full-primer
         * verification.
         */
        ++verified_candidate_positions;

        std::size_t mismatches = 0;
        std::uint64_t mismatch_mask = 0;

        if (
            !SensitiveCandidateSearchAccess::verify_candidate(
                searcher,
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

    normalize_hits(
        hits
    );

    timing.mate_filter_verify_us +=
        elapsed_us(
            filter_verify_start,
            std::chrono::steady_clock::now()
        );

    return hits;
}


PrimerIdentity choose_anchor(
    const SensitiveCandidateCostEstimate& p1,
    const SensitiveCandidateCostEstimate& p2,
    const SensitivePairAnchorPolicy policy
) noexcept {
    if (
        policy ==
        SensitivePairAnchorPolicy::ForcePrimer1
    ) {
        return PrimerIdentity::Primer1;
    }

    if (
        policy ==
        SensitivePairAnchorPolicy::ForcePrimer2
    ) {
        return PrimerIdentity::Primer2;
    }

    /*
     * Validated default policy:
     * select the lower estimated-cost primer
     * as the global anchor.
     */
    if (
        p1.max_seed_occurrences <
        p2.max_seed_occurrences
    ) {
        return PrimerIdentity::Primer1;
    }

    if (
        p2.max_seed_occurrences <
        p1.max_seed_occurrences
    ) {
        return PrimerIdentity::Primer2;
    }

    if (
        p1.total_seed_occurrences <=
        p2.total_seed_occurrences
    ) {
        return PrimerIdentity::Primer1;
    }

    return PrimerIdentity::Primer2;
}

}  // namespace


SensitivePairConstrainedSearchEngine::
SensitivePairConstrainedSearchEngine(
    const BidirectionalFMIndex& index,
    const PackedReference& reference
)
    : anchor_searcher_(
          index,
          reference
      ),
      candidate_searcher_(
          index,
          reference
      ),
      cost_estimator_(
          index
      ),
      reference_(
          reference
      ) {
}


SensitivePairConstrainedSearchResult
SensitivePairConstrainedSearchEngine::search(
    const std::string_view primer1_raw,
    const std::string_view primer2_raw,
    const std::size_t max_mismatches,
    const std::uint64_t min_amplicon_length,
    const std::uint64_t max_amplicon_length,
    const SensitivePairAnchorPolicy anchor_policy
) const {
    if (max_mismatches != 3) {
        throw std::invalid_argument(
            "Pair-constrained SENSITIVE v1 "
            "currently supports exactly k=3."
        );
    }

    if (min_amplicon_length == 0) {
        throw std::invalid_argument(
            "Minimum amplicon length must be > 0."
        );
    }

    if (
        min_amplicon_length >
        max_amplicon_length
    ) {
        throw std::invalid_argument(
            "Minimum amplicon length cannot "
            "exceed maximum amplicon length."
        );
    }

    SensitivePairConstrainedTiming timing{};

    const std::string primer1 =
        normalize_primer(
            primer1_raw
        );

    const std::string primer2 =
        normalize_primer(
            primer2_raw
        );

    /*
     * --------------------------------------------------
     * Cost estimation / anchor selection
     * --------------------------------------------------
     */

    const auto cost_start =
        std::chrono::steady_clock::now();

    const auto p1_cost =
        cost_estimator_.estimate_k3(
            primer1
        );

    const auto p2_cost =
        cost_estimator_.estimate_k3(
            primer2
        );

    const PrimerIdentity anchor_identity =
        choose_anchor(
            p1_cost,
            p2_cost,
            anchor_policy
        );

    timing.cost_estimator_us =
        elapsed_us(
            cost_start,
            std::chrono::steady_clock::now()
        );

    const std::string_view anchor_primer =
        anchor_identity ==
            PrimerIdentity::Primer1
            ? std::string_view(
                  primer1
              )
            : std::string_view(
                  primer2
              );

    const std::string_view mate_primer =
        anchor_identity ==
            PrimerIdentity::Primer1
            ? std::string_view(
                  primer2
              )
            : std::string_view(
                  primer1
              );

    /*
     * --------------------------------------------------
     * Global anchor search
     * --------------------------------------------------
     */

    const auto anchor_start =
        std::chrono::steady_clock::now();

    const auto anchor_result =
        anchor_searcher_.search(
            anchor_primer,
            max_mismatches
        );

    timing.anchor_search_us =
        elapsed_us(
            anchor_start,
            std::chrono::steady_clock::now()
        );

    /*
     * --------------------------------------------------
     * PCR-compatible mate windows
     * --------------------------------------------------
     */

    const auto windows_start =
        std::chrono::steady_clock::now();

    std::vector<StartWindow>
        forward_mate_windows;

    std::vector<StartWindow>
        reverse_mate_windows;

    for (
        const auto& hit :
        anchor_result
            .search_result
            .hits
    ) {
        switch (hit.orientation) {

            case PrimerOrientation::Forward:

                add_downstream_window(
                    hit,
                    anchor_primer.size(),
                    mate_primer.size(),
                    reference_.size(),
                    min_amplicon_length,
                    max_amplicon_length,
                    reverse_mate_windows
                );

                break;


            case PrimerOrientation::Reverse:

                add_upstream_window(
                    hit,
                    anchor_primer.size(),
                    mate_primer.size(),
                    reference_.size(),
                    min_amplicon_length,
                    max_amplicon_length,
                    forward_mate_windows
                );

                break;
        }
    }

    forward_mate_windows =
        merge_windows(
            std::move(
                forward_mate_windows
            )
        );

    reverse_mate_windows =
        merge_windows(
            std::move(
                reverse_mate_windows
            )
        );

    timing.window_build_merge_us =
        elapsed_us(
            windows_start,
            std::chrono::steady_clock::now()
        );

    /*
     * --------------------------------------------------
     * Window-filtered half-seed mate search
     * --------------------------------------------------
     */

    const std::string mate_reverse =
        reverse_complement(
            mate_primer
        );

    std::uint64_t verified_positions =
        0;

    auto mate_forward_hits =
        search_candidates_in_windows(
            candidate_searcher_,
            forward_mate_windows,
            mate_primer,
            PrimerOrientation::Forward,
            false,
            max_mismatches,
            verified_positions,
            timing
        );

    auto mate_reverse_hits =
        search_candidates_in_windows(
            candidate_searcher_,
            reverse_mate_windows,
            mate_reverse,
            PrimerOrientation::Reverse,
            true,
            max_mismatches,
            verified_positions,
            timing
        );

    std::vector<
        OrientedPrimerSearchHit
    > mate_hits;

    mate_hits.reserve(
        mate_forward_hits.size() +
        mate_reverse_hits.size()
    );

    mate_hits.insert(
        mate_hits.end(),
        mate_forward_hits.begin(),
        mate_forward_hits.end()
    );

    mate_hits.insert(
        mate_hits.end(),
        mate_reverse_hits.begin(),
        mate_reverse_hits.end()
    );

    normalize_hits(
        mate_hits
    );

    /*
     * --------------------------------------------------
     * Validated pair assembly kernel
     * --------------------------------------------------
     */

    const auto pair_start =
        std::chrono::steady_clock::now();

    const auto& anchor_hits =
        anchor_result
            .search_result
            .hits;

    PrimerPairSearchResult
        pair_result;

    if (
        anchor_identity ==
        PrimerIdentity::Primer1
    ) {
        pair_result =
            PrimerPairSearchEngine::
                assemble_pairs(
                    primer1,
                    anchor_hits,

                    primer2,
                    mate_hits,

                    min_amplicon_length,
                    max_amplicon_length
                );

    } else {

        pair_result =
            PrimerPairSearchEngine::
                assemble_pairs(
                    primer1,
                    mate_hits,

                    primer2,
                    anchor_hits,

                    min_amplicon_length,
                    max_amplicon_length
                );
    }

    timing.pair_assembly_us =
        elapsed_us(
            pair_start,
            std::chrono::steady_clock::now()
        );

    SensitivePairConstrainedSearchResult
        output{};

    output.anchor_primer =
        anchor_identity;

    output.anchor_global_hit_count =
        anchor_hits.size();

    output.mate_local_hit_count =
        mate_hits.size();

    output.merged_forward_mate_windows =
        forward_mate_windows.size();

    output.merged_reverse_mate_windows =
        reverse_mate_windows.size();

    output.scanned_mate_start_positions =
        verified_positions;

    output.primer1_cost =
        p1_cost;

    output.primer2_cost =
        p2_cost;

    output.timing =
        timing;

    output.pair_result =
        std::move(
            pair_result
        );

    return output;
}

}  // namespace primerpair
