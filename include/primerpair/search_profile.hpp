#pragma once

#include <cstddef>
#include <string_view>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/sensitive_adaptive_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {

enum class SearchProfile {
    Strict,
    Sensitive
};

[[nodiscard]]
const char* to_string(
    SearchProfile profile
) noexcept;


/*
 * Profile façade result.
 *
 * Backward compatibility:
 *
 * STRICT continues to populate search_result exactly
 * as before.
 *
 * SENSITIVE populates sensitive_result because its
 * backend metadata is intentionally different from
 * StrandAwarePrimerSearchResult.
 *
 * New generic callers should use:
 *
 *     hits()
 *     hit_count()
 *     primer_length()
 *
 * rather than assuming one concrete backend result.
 */
struct ProfiledPrimerSearchResult {
    SearchProfile profile{
        SearchProfile::Strict
    };

    /*
     * Existing STRICT result.
     *
     * Kept unchanged for backward compatibility.
     */
    StrandAwarePrimerSearchResult
        search_result{};

    /*
     * Populated only for SENSITIVE.
     */
    SensitiveAdaptiveSearchResult
        sensitive_result{};

    bool sensitive_result_available{
        false
    };

    [[nodiscard]]
    const auto& hits() const noexcept {
        if (sensitive_result_available) {
            return
                sensitive_result
                    .search_result
                    .hits;
        }

        return
            search_result.hits;
    }

    [[nodiscard]]
    std::size_t hit_count() const noexcept {
        return hits().size();
    }

    [[nodiscard]]
    std::size_t primer_length() const noexcept {
        if (sensitive_result_available) {
            return
                sensitive_result
                    .search_result
                    .primer_length;
        }

        return
            search_result
                .primer_length;
    }
};


class ProfiledPrimerSearchEngine {
public:
    ProfiledPrimerSearchEngine(
        const BidirectionalFMIndex& index,
        const PackedReference& reference
    );

    /*
     * STRICT
     * ------
     *
     * Uses the already validated biological model:
     *
     *     exact 3-prime anchor
     *     +
     *     5-prime mismatch allowance
     *
     *
     * SENSITIVE
     * ---------
     *
     * Uses the independently validated adaptive
     * full-primer Hamming search:
     *
     *     k=0,1,2
     *         exhaustive reference backend
     *
     *     k=3
     *         adaptive MAX-seed router
     *
     *         max_seed <= 12275
     *             candidate backend
     *
     *         max_seed > 12275
     *             exhaustive backend
     *
     *
     * anchor_length is meaningful only for STRICT.
     * SENSITIVE deliberately searches the entire
     * primer and therefore does not use this value.
     */
    [[nodiscard]]
    ProfiledPrimerSearchResult search(
        std::string_view primer,
        SearchProfile profile,
        std::size_t anchor_length = 12,
        std::size_t max_mismatches = 3
    ) const;

private:
    StrandAwarePrimerSearchEngine
        strict_searcher_;

    SensitiveAdaptiveSearchEngine
        sensitive_searcher_;
};

}  // namespace primerpair
