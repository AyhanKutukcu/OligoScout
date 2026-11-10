#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <span>
#include <vector>

#include "primerpair/fm_index.hpp"

namespace primerpair {

/*
 * Interleaved Prefix BWT (IP-BWT).
 *
 * Construction:
 *
 *   shared suffix array
 *          ↓
 *   K-symbol prefix
 *          +
 *   K-shifted suffix-array row
 *
 * Prefixes are encoded numerically using 3 bits
 * per symbol:
 *
 *   $ = 0
 *   A = 1
 *   C = 2
 *   G = 3
 *   N = 4
 *   T = 5
 *
 * Therefore K <= 21 fits inside uint64_t:
 *
 *   21 * 3 = 63 bits.
 *
 * Storage is structure-of-arrays:
 *
 *   prefix_codes_  uint64_t
 *   paired_rows_   uint32_t
 *   suffix_array_  uint32_t
 *
 * No std::string is stored per IP-BWT row.
 */
struct IPBWTCertifiedWindowResult {
    Interval interval{};
    bool used_global_fallback{false};
    std::uint64_t window_begin{0};
    std::uint64_t window_end{0};
};


class IPBWTIndex {
public:
    static constexpr std::size_t
        max_chunk_length = 21;

    explicit IPBWTIndex(
        std::string reference,
        std::size_t chunk_length
    );

    [[nodiscard]]
    Interval exact_search(
        std::string_view query
    ) const;

    [[nodiscard]]
    std::vector<std::uint64_t>
    locate(
        const Interval& interval
    ) const;

    [[nodiscard]]
    std::uint64_t lower_bound(
        std::string_view chunk,
        std::uint64_t row
    ) const;


    /*
     * Four independent IP-BWT lower-bound
     * searches executed in an interleaved loop.
     *
     * Phase 4A is deliberately scalar.
     * Phase 4B will replace the comparison stage
     * with AVX2 while preserving this API.
     */
    [[nodiscard]]
    std::array<std::uint64_t, 4>
    lower_bound_batch4(
        const std::array<std::string_view, 4>& chunks,
        const std::array<std::uint64_t, 4>& rows
    ) const;


    /*
     * Exact one-chunk prefix search accelerated by
     * a predicted suffix-array lower-bound.
     *
     * The supplied window is accepted only when its
     * left/right boundaries mathematically bracket the
     * true lower-bound. Otherwise the implementation
     * falls back to the global exact binary search.
     *
     * Correctness therefore does not depend on the
     * predictor.
     */
    [[nodiscard]]
    IPBWTCertifiedWindowResult
    exact_prefix_search_certified_window(
        std::string_view query,
        std::uint64_t predicted_lower,
        std::uint64_t radius
    ) const;


    /*
     * Eight independent exact 21-mer prefix searches
     * with interleaved binary-search progress.
     *
     * This is the fair throughput baseline for the
     * neural Batch8 path.
     */
    [[nodiscard]]
    std::array<Interval, 8>
    exact_prefix_search_batch8(
        const std::array<std::string_view, 8>& queries
    ) const;


    /*
     * Production throughput API.
     *
     * Queries are processed in interleaved groups
     * of eight. A final partial group is handled
     * exactly with the scalar path.
     *
     * Every query must contain exactly one complete
     * IP-BWT chunk.
     */
    [[nodiscard]]
    std::vector<Interval>
    exact_prefix_search_many(
        const std::vector<std::string_view>& queries
    ) const;


    /*
     * Allocation-free production hot path.
     *
     * results.size() must equal queries.size().
     * The caller owns and may reuse the output buffer.
     */
    void
    exact_prefix_search_many(
        std::span<const std::string_view> queries,
        std::span<Interval> results
    ) const;


    /*
     * Batch8 version of certified neural-window search.
     *
     * Every lane is exact. A failed neural window
     * certification turns that lane into a global
     * binary search while the other lanes continue
     * within their local windows.
     */
    [[nodiscard]]
    std::array<IPBWTCertifiedWindowResult, 8>
    exact_prefix_search_certified_window_batch8(
        const std::array<std::string_view, 8>& queries,
        const std::array<std::uint64_t, 8>& predicted_lowers,
        std::uint64_t radius
    ) const;

    [[nodiscard]]
    std::size_t chunk_length() const noexcept {
        return chunk_length_;
    }

    [[nodiscard]]
    std::uint64_t row_count() const noexcept {
        return static_cast<std::uint64_t>(
            prefix_codes_.size()
        );
    }

    [[nodiscard]]
    std::size_t reference_length() const noexcept {
        return reference_length_;
    }

    [[nodiscard]]
    std::size_t compact_storage_bytes() const noexcept {
        return
            prefix_codes_.capacity() *
                sizeof(std::uint64_t)
            +
            paired_rows_.capacity() *
                sizeof(std::uint32_t)
            +
            suffix_array_.capacity() *
                sizeof(std::uint32_t);
    }

private:
    friend class IPBWTRMI;
    std::string reference_;
    std::string text_;

    std::size_t reference_length_{0};
    std::size_t chunk_length_{0};

    /*
     * Structure-of-arrays avoids padding that a
     * {uint64_t,uint32_t} struct may introduce.
     */
    std::vector<std::uint64_t>
        prefix_codes_;

    std::vector<std::uint32_t>
        paired_rows_;

    /*
     * Used by locate().
     *
     * Later this can be replaced by/reconnected to
     * sampled-SA resolution when persistent learned
     * indexes are implemented.
     */
    std::vector<std::uint32_t>
        suffix_array_;

    [[nodiscard]]
    char cyclic_base(
        std::uint64_t rotation_start,
        std::size_t offset
    ) const;

    [[nodiscard]]
    std::uint64_t prefix_code(
        std::uint64_t rotation_start
    ) const;

    [[nodiscard]]
    std::uint64_t encode_chunk(
        std::string_view chunk
    ) const;

    [[nodiscard]]
    static std::uint64_t symbol_code(
        char base
    );

    static void validate_reference(
        std::string_view reference
    );

    static void validate_query(
        std::string_view query
    );
};

}  // namespace primerpair
