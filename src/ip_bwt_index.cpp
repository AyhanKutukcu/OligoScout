#include "primerpair/ip_bwt_index.hpp"

#include "primerpair/suffix_array_builder.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace primerpair {


std::uint64_t IPBWTIndex::symbol_code(
    const char base
) {
    switch (base) {

        case '$':
            return 0;

        case 'A':
            return 1;

        case 'C':
            return 2;

        case 'G':
            return 3;

        case 'N':
            return 4;

        case 'T':
            return 5;

        default:
            throw std::invalid_argument(
                "Unsupported IP-BWT nucleotide."
            );
    }
}


IPBWTIndex::IPBWTIndex(
    std::string reference,
    const std::size_t chunk_length
)
    : reference_(
          std::move(
              reference
          )
      ),
      text_(
          reference_ + "$"
      ),
      reference_length_(
          reference_.size()
      ),
      chunk_length_(
          chunk_length
      ) {

    validate_reference(
        reference_
    );


    if (chunk_length_ == 0) {
        throw std::invalid_argument(
            "IP-BWT chunk length cannot be zero."
        );
    }


    if (
        chunk_length_ >
        max_chunk_length
    ) {
        throw std::invalid_argument(
            "IP-BWT chunk length exceeds "
            "63-bit packed-prefix capacity."
        );
    }


    if (
        chunk_length_ >
        text_.size()
    ) {
        throw std::invalid_argument(
            "IP-BWT chunk length exceeds "
            "indexed text length."
        );
    }


    /*
     * Shared scalable suffix-array construction.
     */
    suffix_array_ =
        build_suffix_array_prefix_doubling(
            text_
        );


    const std::size_t n =
        suffix_array_.size();


    /*
     * inverse_row[start] = SA/BW row.
     */
    std::vector<std::uint32_t>
        inverse_row(
            n,
            std::uint32_t{0}
        );


    for (
        std::size_t row = 0;
        row < n;
        ++row
    ) {
        const std::uint32_t start =
            suffix_array_.at(
                row
            );

        inverse_row.at(
            static_cast<std::size_t>(
                start
            )
        ) =
            static_cast<std::uint32_t>(
                row
            );
    }


    prefix_codes_.resize(
        n
    );

    paired_rows_.resize(
        n
    );


    /*
     * Build numeric IP-BWT compound keys:
     *
     *   (packed K-prefix, paired SA row)
     */
    for (
        std::size_t row = 0;
        row < n;
        ++row
    ) {
        const std::uint32_t start =
            suffix_array_.at(
                row
            );


        prefix_codes_.at(
            row
        ) =
            prefix_code(
                start
            );


        const std::size_t shifted_start =
            (
                static_cast<std::size_t>(
                    start
                )
                +
                chunk_length_
            )
            %
            n;


        paired_rows_.at(
            row
        ) =
            inverse_row.at(
                shifted_start
            );
    }


    /*
     * IP-BWT invariant:
     *
     * Compound keys must already be sorted by:
     *
     *   prefix_code
     *   paired_row
     *
     * SA ordering guarantees this if construction
     * is correct.
     *
     * Do not silently sort here.
     */
    for (
        std::size_t row = 1;
        row < n;
        ++row
    ) {
        const std::uint64_t previous_prefix =
            prefix_codes_.at(
                row - 1
            );

        const std::uint64_t current_prefix =
            prefix_codes_.at(
                row
            );


        const std::uint32_t previous_paired =
            paired_rows_.at(
                row - 1
            );

        const std::uint32_t current_paired =
            paired_rows_.at(
                row
            );


        const bool sorted =
            (
                previous_prefix <
                current_prefix
            )
            ||
            (
                previous_prefix ==
                    current_prefix
                &&
                previous_paired <=
                    current_paired
            );


        if (!sorted) {
            throw std::logic_error(
                "Compact IP-BWT compound keys "
                "are not sorted."
            );
        }
    }


    /*
     * Build-only text/reference strings are no longer
     * required after packed prefixes are materialized.
     *
     * Keeping reference_length_ is sufficient for
     * filtering the sentinel row during locate().
     */
    std::string{}.swap(
        reference_
    );

    std::string{}.swap(
        text_
    );
}


char IPBWTIndex::cyclic_base(
    const std::uint64_t rotation_start,
    const std::size_t offset
) const {
    if (text_.empty()) {
        throw std::logic_error(
            "IP-BWT build text has already "
            "been released."
        );
    }


    const std::uint64_t n =
        static_cast<std::uint64_t>(
            text_.size()
        );


    const std::uint64_t position =
        (
            rotation_start
            +
            static_cast<std::uint64_t>(
                offset
            )
        )
        %
        n;


    return
        text_.at(
            static_cast<std::size_t>(
                position
            )
        );
}


std::uint64_t IPBWTIndex::prefix_code(
    const std::uint64_t rotation_start
) const {
    std::uint64_t code = 0;


    for (
        std::size_t offset = 0;
        offset < chunk_length_;
        ++offset
    ) {
        code <<= 3;

        code |=
            symbol_code(
                cyclic_base(
                    rotation_start,
                    offset
                )
            );
    }


    return code;
}


std::uint64_t IPBWTIndex::encode_chunk(
    const std::string_view chunk
) const {
    if (
        chunk.size() !=
        chunk_length_
    ) {
        throw std::invalid_argument(
            "IP-BWT lookup chunk has "
            "wrong length."
        );
    }


    std::uint64_t code = 0;


    for (
        const char base :
        chunk
    ) {
        code <<= 3;

        code |=
            symbol_code(
                base
            );
    }


    return code;
}



std::array<std::uint64_t, 4>
IPBWTIndex::lower_bound_batch4(
    const std::array<std::string_view, 4>& chunks,
    const std::array<std::uint64_t, 4>& rows
) const {

    std::array<std::uint64_t, 4>
        query_prefixes{};

    std::array<std::size_t, 4>
        low{};

    std::array<std::size_t, 4>
        high{};

    std::array<bool, 4>
        active{};


    const std::size_t n =
        prefix_codes_.size();


    for (
        std::size_t lane = 0;
        lane < 4;
        ++lane
    ) {
        if (
            chunks.at(lane).size() !=
            chunk_length_
        ) {
            throw std::invalid_argument(
                "IP-BWT batch lookup chunk "
                "has wrong length."
            );
        }


        if (
            rows.at(lane) >
            prefix_codes_.size()
        ) {
            throw std::out_of_range(
                "IP-BWT batch lookup row "
                "is out of range."
            );
        }


        query_prefixes.at(lane) =
            encode_chunk(
                chunks.at(lane)
            );

        low.at(lane) = 0;

        high.at(lane) = n;

        active.at(lane) =
            n != 0;
    }


    bool any_active = true;


    while (
        any_active
    ) {
        any_active = false;


        for (
            std::size_t lane = 0;
            lane < 4;
            ++lane
        ) {
            if (
                !active.at(lane)
            ) {
                continue;
            }


            const std::size_t middle =
                low.at(lane)
                +
                (
                    high.at(lane)
                    -
                    low.at(lane)
                )
                /
                2;


            const std::uint64_t prefix =
                prefix_codes_.at(
                    middle
                );


            const std::uint32_t paired =
                paired_rows_.at(
                    middle
                );


            const bool less =
                (
                    prefix <
                    query_prefixes.at(lane)
                )
                ||
                (
                    prefix ==
                        query_prefixes.at(lane)
                    &&
                    static_cast<std::uint64_t>(
                        paired
                    )
                    <
                    rows.at(lane)
                );


            if (less) {
                low.at(lane) =
                    middle + 1;
            } else {
                high.at(lane) =
                    middle;
            }


            if (
                low.at(lane) <
                high.at(lane)
            ) {
                any_active = true;
            } else {
                active.at(lane) = false;
            }
        }
    }


    return {
        static_cast<std::uint64_t>(
            low.at(0)
        ),
        static_cast<std::uint64_t>(
            low.at(1)
        ),
        static_cast<std::uint64_t>(
            low.at(2)
        ),
        static_cast<std::uint64_t>(
            low.at(3)
        )
    };
}


std::uint64_t IPBWTIndex::lower_bound(
    const std::string_view chunk,
    const std::uint64_t row
) const {
    if (
        chunk.size() !=
        chunk_length_
    ) {
        throw std::invalid_argument(
            "IP-BWT lookup chunk has "
            "wrong length."
        );
    }


    if (
        row >
        prefix_codes_.size()
    ) {
        throw std::out_of_range(
            "IP-BWT lookup row is out of range."
        );
    }


    const std::uint64_t query_prefix =
        encode_chunk(
            chunk
        );


    /*
     * Manual binary search over the numeric
     * compound key:
     *
     *   (prefix_code, paired_row)
     */
    std::size_t low = 0;

    std::size_t high =
        prefix_codes_.size();


    while (low < high) {

        const std::size_t middle =
            low +
            (
                high - low
            ) /
            2;


        const std::uint64_t prefix =
            prefix_codes_.at(
                middle
            );

        const std::uint32_t paired =
            paired_rows_.at(
                middle
            );


        const bool less =
            (
                prefix <
                query_prefix
            )
            ||
            (
                prefix ==
                    query_prefix
                &&
                static_cast<std::uint64_t>(
                    paired
                )
                <
                row
            );


        if (less) {
            low =
                middle + 1;
        } else {
            high =
                middle;
        }
    }


    return
        static_cast<std::uint64_t>(
            low
        );
}



std::array<Interval, 8>
IPBWTIndex::exact_prefix_search_batch8(
    const std::array<std::string_view, 8>& queries
) const {
    std::array<Interval, 8>
        results{};


    const std::uint64_t n =
        static_cast<std::uint64_t>(
            prefix_codes_.size()
        );


    if (
        n == 0
    ) {
        return results;
    }


    std::array<std::uint64_t, 8>
        prefixes{};

    std::array<std::uint64_t, 8>
        low{};

    std::array<std::uint64_t, 8>
        high{};


    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        if (
            queries[lane].size() !=
            chunk_length_
        ) {
            throw std::invalid_argument(
                "Batch8 exact prefix search requires "
                "one complete IP-BWT chunk per lane."
            );
        }


        prefixes[lane] =
            encode_chunk(
                queries[lane]
            );


        low[lane] = 0;
        high[lane] = n;
    }


    /*
     * Interleave the eight independent binary searches.
     *
     * This creates instruction-level and memory-level
     * parallelism without changing exact semantics.
     */
    bool active = true;


    while (
        active
    ) {
        active = false;


        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            if (
                low[lane] >=
                high[lane]
            ) {
                continue;
            }


            active = true;


            const std::uint64_t middle =
                low[lane]
                +
                (
                    high[lane]
                    -
                    low[lane]
                )
                /
                2;


            if (
                prefix_codes_[
                    static_cast<std::size_t>(
                        middle
                    )
                ]
                <
                prefixes[lane]
            ) {
                low[lane] =
                    middle + 1;

            } else {

                high[lane] =
                    middle;
            }
        }
    }


    /*
     * Lower bounds are now exact.
     *
     * Most 21-mers are unique. Repetitive blocks use
     * exponential expansion followed by binary search.
     */
    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        const std::uint64_t lower =
            low[lane];

        const std::uint64_t prefix =
            prefixes[lane];


        if (
            lower >= n
            ||
            prefix_codes_[
                static_cast<std::size_t>(
                    lower
                )
            ]
            !=
            prefix
        ) {
            results[lane] =
                Interval{
                    lower,
                    lower
                };

            continue;
        }


        if (
            lower + 1 >= n
            ||
            prefix_codes_[
                static_cast<std::size_t>(
                    lower + 1
                )
            ]
            !=
            prefix
        ) {
            results[lane] =
                Interval{
                    lower,
                    lower + 1
                };

            continue;
        }


        std::uint64_t upper_low =
            lower + 2;

        std::uint64_t step = 2;

        std::uint64_t upper_high =
            n;


        while (
            true
        ) {
            if (
                step >
                n -
                lower -
                1
            ) {
                upper_high =
                    n;

                break;
            }


            const std::uint64_t probe =
                lower +
                step;


            if (
                prefix_codes_[
                    static_cast<std::size_t>(
                        probe
                    )
                ]
                !=
                prefix
            ) {
                upper_high =
                    probe;

                break;
            }


            upper_low =
                probe + 1;


            if (
                upper_low >= n
            ) {
                upper_high =
                    n;

                break;
            }


            if (
                step >
                (
                    n -
                    lower
                )
                /
                2
            ) {
                upper_high =
                    n;

                break;
            }


            step *= 2;
        }


        while (
            upper_low <
            upper_high
        ) {
            const std::uint64_t middle =
                upper_low
                +
                (
                    upper_high
                    -
                    upper_low
                )
                /
                2;


            if (
                prefix_codes_[
                    static_cast<std::size_t>(
                        middle
                    )
                ]
                <=
                prefix
            ) {
                upper_low =
                    middle + 1;

            } else {

                upper_high =
                    middle;
            }
        }


        results[lane] =
            Interval{
                lower,
                upper_low
            };
    }


    return results;
}



void
IPBWTIndex::exact_prefix_search_many(
    const std::span<const std::string_view> queries,
    const std::span<Interval> results
) const {
    if (
        results.size() !=
        queries.size()
    ) {
        throw std::invalid_argument(
            "exact_prefix_search_many requires "
            "results.size() == queries.size()."
        );
    }


    const std::size_t full =
        (
            queries.size()
            /
            8
        )
        *
        8;


    std::size_t i = 0;


    for (
        ;
        i < full;
        i += 8
    ) {
        const std::array<std::string_view, 8>
        batch{
            queries[i + 0],
            queries[i + 1],
            queries[i + 2],
            queries[i + 3],
            queries[i + 4],
            queries[i + 5],
            queries[i + 6],
            queries[i + 7]
        };


        const auto batch_results =
            exact_prefix_search_batch8(
                batch
            );


        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            results[
                i + lane
            ] =
                batch_results[
                    lane
                ];
        }
    }


    /*
     * Final partial group.
     */
    for (
        ;
        i < queries.size();
        ++i
    ) {
        if (
            queries[i].size() !=
            chunk_length_
        ) {
            throw std::invalid_argument(
                "exact_prefix_search_many requires "
                "one complete IP-BWT chunk per query."
            );
        }


        results[i] =
            exact_search(
                queries[i]
            );
    }
}


std::vector<Interval>
IPBWTIndex::exact_prefix_search_many(
    const std::vector<std::string_view>& queries
) const {
    std::vector<Interval>
        results(
            queries.size()
        );


    exact_prefix_search_many(
        std::span<const std::string_view>(
            queries.data(),
            queries.size()
        ),
        std::span<Interval>(
            results.data(),
            results.size()
        )
    );


    return results;
}


std::array<IPBWTCertifiedWindowResult, 8>
IPBWTIndex::exact_prefix_search_certified_window_batch8(
    const std::array<std::string_view, 8>& queries,
    const std::array<std::uint64_t, 8>& predicted_lowers,
    const std::uint64_t radius
) const {
    std::array<IPBWTCertifiedWindowResult, 8>
        results{};


    const std::uint64_t n =
        static_cast<std::uint64_t>(
            prefix_codes_.size()
        );


    if (
        n == 0
    ) {
        return results;
    }


    std::array<std::uint64_t, 8>
        prefixes{};

    std::array<std::uint64_t, 8>
        low{};

    std::array<std::uint64_t, 8>
        high{};

    std::array<bool, 8>
        fallback{};


    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        if (
            queries[lane].size() !=
            chunk_length_
        ) {
            throw std::invalid_argument(
                "Certified Batch8 search requires "
                "one complete IP-BWT chunk per lane."
            );
        }


        const std::uint64_t prefix =
            encode_chunk(
                queries[lane]
            );


        prefixes[lane] =
            prefix;


        const std::uint64_t predicted =
            std::min(
                predicted_lowers[lane],
                n
            );


        const std::uint64_t begin =
            (
                predicted >
                radius
            )
            ?
            predicted -
                radius
            :
            0;


        const std::uint64_t end =
            std::min(
                n,
                predicted +
                    radius +
                    1
            );


        const bool left_certified =
            (
                begin == 0
            )
            ||
            (
                prefix_codes_[
                    static_cast<std::size_t>(
                        begin - 1
                    )
                ]
                <
                prefix
            );


        const bool right_certified =
            (
                end == n
            )
            ||
            (
                prefix_codes_[
                    static_cast<std::size_t>(
                        end
                    )
                ]
                >=
                prefix
            );


        fallback[lane] =
            !(
                left_certified
                &&
                right_certified
            );


        if (
            fallback[lane]
        ) {
            /*
             * Exact global fallback, but still
             * interleaved with all other lanes.
             */
            low[lane] = 0;
            high[lane] = n;

        } else {

            low[lane] = begin;
            high[lane] = end;
        }


        results[lane].used_global_fallback =
            fallback[lane];

        results[lane].window_begin =
            begin;

        results[lane].window_end =
            end;
    }


    /*
     * Local and fallback lanes execute together.
     *
     * A fallback lane simply has [0,n) as its search
     * region; certified lanes retain their narrow
     * neural windows.
     */
    bool active = true;


    while (
        active
    ) {
        active = false;


        for (
            std::size_t lane = 0;
            lane < 8;
            ++lane
        ) {
            if (
                low[lane] >=
                high[lane]
            ) {
                continue;
            }


            active = true;


            const std::uint64_t middle =
                low[lane]
                +
                (
                    high[lane]
                    -
                    low[lane]
                )
                /
                2;


            if (
                prefix_codes_[
                    static_cast<std::size_t>(
                        middle
                    )
                ]
                <
                prefixes[lane]
            ) {
                low[lane] =
                    middle + 1;

            } else {

                high[lane] =
                    middle;
            }
        }
    }


    for (
        std::size_t lane = 0;
        lane < 8;
        ++lane
    ) {
        const std::uint64_t lower =
            low[lane];

        const std::uint64_t prefix =
            prefixes[lane];


        if (
            lower >= n
            ||
            prefix_codes_[
                static_cast<std::size_t>(
                    lower
                )
            ]
            !=
            prefix
        ) {
            results[lane].interval =
                Interval{
                    lower,
                    lower
                };

            continue;
        }


        if (
            lower + 1 >= n
            ||
            prefix_codes_[
                static_cast<std::size_t>(
                    lower + 1
                )
            ]
            !=
            prefix
        ) {
            results[lane].interval =
                Interval{
                    lower,
                    lower + 1
                };

            continue;
        }


        std::uint64_t upper_low =
            lower + 2;

        std::uint64_t step = 2;

        std::uint64_t upper_high =
            n;


        while (
            true
        ) {
            if (
                step >
                n -
                lower -
                1
            ) {
                upper_high =
                    n;

                break;
            }


            const std::uint64_t probe =
                lower +
                step;


            if (
                prefix_codes_[
                    static_cast<std::size_t>(
                        probe
                    )
                ]
                !=
                prefix
            ) {
                upper_high =
                    probe;

                break;
            }


            upper_low =
                probe + 1;


            if (
                upper_low >= n
            ) {
                upper_high =
                    n;

                break;
            }


            if (
                step >
                (
                    n -
                    lower
                )
                /
                2
            ) {
                upper_high =
                    n;

                break;
            }


            step *= 2;
        }


        while (
            upper_low <
            upper_high
        ) {
            const std::uint64_t middle =
                upper_low
                +
                (
                    upper_high
                    -
                    upper_low
                )
                /
                2;


            if (
                prefix_codes_[
                    static_cast<std::size_t>(
                        middle
                    )
                ]
                <=
                prefix
            ) {
                upper_low =
                    middle + 1;

            } else {

                upper_high =
                    middle;
            }
        }


        results[lane].interval =
            Interval{
                lower,
                upper_low
            };
    }


    return results;
}


IPBWTCertifiedWindowResult
IPBWTIndex::exact_prefix_search_certified_window(
    const std::string_view query,
    const std::uint64_t predicted_lower,
    const std::uint64_t radius
) const {

    if (
        query.size() !=
        chunk_length_
    ) {
        throw std::invalid_argument(
            "Certified neural-window search requires "
            "exactly one full IP-BWT chunk."
        );
    }


    const std::uint64_t prefix =
        encode_chunk(
            query
        );


    const std::uint64_t n =
        static_cast<std::uint64_t>(
            prefix_codes_.size()
        );


    if (
        n == 0
    ) {
        return {
            Interval{0, 0},
            false,
            0,
            0
        };
    }


    const std::uint64_t predicted =
        std::min(
            predicted_lower,
            n
        );


    const std::uint64_t begin =
        (
            predicted >
            radius
        )
        ?
        predicted -
            radius
        :
        0;


    const std::uint64_t end =
        std::min(
            n,
            predicted +
                radius +
                1
        );


    /*
     * For row=0, the compound lower-bound is the
     * first entry whose prefix is >= query prefix.
     *
     * We can certify that the true lower-bound lies
     * in [begin,end] using only the two boundaries.
     */
    const bool left_certified =
        (
            begin == 0
        )
        ||
        (
            prefix_codes_.at(
                static_cast<std::size_t>(
                    begin - 1
                )
            )
            <
            prefix
        );


    const bool right_certified =
        (
            end == n
        )
        ||
        (
            prefix_codes_.at(
                static_cast<std::size_t>(
                    end
                )
            )
            >=
            prefix
        );


    bool used_global_fallback =
        !(
            left_certified
            &&
            right_certified
        );


    std::uint64_t lower = 0;


    if (
        !used_global_fallback
    ) {
        std::uint64_t low =
            begin;

        std::uint64_t high =
            end;


        while (
            low <
            high
        ) {
            const std::uint64_t middle =
                low +
                (
                    high -
                    low
                )
                /
                2;


            if (
                prefix_codes_.at(
                    static_cast<std::size_t>(
                        middle
                    )
                )
                <
                prefix
            ) {
                low =
                    middle + 1;

            } else {

                high =
                    middle;
            }
        }


        lower =
            low;

    } else {

        /*
         * Correctness fallback.
         *
         * lower_bound(query,0) is the exact global
         * beginning of this prefix block.
         */
        lower =
            lower_bound(
                query,
                0
            );
    }


    /*
     * Query absent:
     *
     * lower is still the exact suffix-array insertion
     * point, but there is no equal-prefix interval.
     */
    if (
        lower >= n
        ||
        prefix_codes_.at(
            static_cast<std::size_t>(
                lower
            )
        )
        !=
        prefix
    ) {
        return {
            Interval{
                lower,
                lower
            },
            used_global_fallback,
            begin,
            end
        };
    }


    /*
     * Find the END of the equal-prefix block.
     *
     * Most 21-mers are unique, so the common case
     * terminates after one neighboring-row check.
     *
     * Repetitive k-mers use exponential expansion,
     * followed by a small binary search.
     */
    if (
        lower + 1 >= n
        ||
        prefix_codes_.at(
            static_cast<std::size_t>(
                lower + 1
            )
        )
        !=
        prefix
    ) {
        return {
            Interval{
                lower,
                lower + 1
            },
            used_global_fallback,
            begin,
            end
        };
    }


    std::uint64_t upper_low =
        lower + 2;

    std::uint64_t step =
        2;

    std::uint64_t upper_high =
        n;


    while (true) {

        if (
            step >
            n -
            lower -
            1
        ) {
            upper_high =
                n;

            break;
        }


        const std::uint64_t probe =
            lower +
            step;


        if (
            prefix_codes_.at(
                static_cast<std::size_t>(
                    probe
                )
            )
            !=
            prefix
        ) {
            upper_high =
                probe;

            break;
        }


        upper_low =
            probe + 1;


        if (
            upper_low >=
            n
        ) {
            return {
                Interval{
                    lower,
                    n
                },
                used_global_fallback,
                begin,
                end
            };
        }


        if (
            step >
            (
                n -
                lower
            )
            /
            2
        ) {
            upper_high =
                n;

            break;
        }


        step *= 2;
    }


    while (
        upper_low <
        upper_high
    ) {
        const std::uint64_t middle =
            upper_low +
            (
                upper_high -
                upper_low
            )
            /
            2;


        if (
            prefix_codes_.at(
                static_cast<std::size_t>(
                    middle
                )
            )
            <=
            prefix
        ) {
            upper_low =
                middle + 1;

        } else {

            upper_high =
                middle;
        }
    }


    return {
        Interval{
            lower,
            upper_low
        },
        used_global_fallback,
        begin,
        end
    };
}


Interval IPBWTIndex::exact_search(
    const std::string_view query
) const {
    validate_query(
        query
    );


    const std::uint64_t n =
        static_cast<std::uint64_t>(
            prefix_codes_.size()
        );


    std::uint64_t low = 0;
    std::uint64_t high = n;


    /*
     * Split query into K-sized chunks.
     */
    std::vector<std::string_view>
        chunks;


    for (
        std::size_t position = 0;
        position < query.size();
        position += chunk_length_
    ) {
        const std::size_t remaining =
            query.size() -
            position;


        const std::size_t width =
            std::min(
                chunk_length_,
                remaining
            );


        chunks.push_back(
            query.substr(
                position,
                width
            )
        );
    }


    /*
     * Process from query suffix toward query prefix,
     * analogous to FM backward search.
     */
    for (
        auto it = chunks.rbegin();
        it != chunks.rend();
        ++it
    ) {
        const std::string_view chunk =
            *it;


        if (
            chunk.size() <
            chunk_length_
        ) {

            /*
             * Lowest possible completion.
             *
             * Sentinel is lexicographically lowest.
             * Following positions use A, the lowest
             * non-sentinel reference symbol.
             */
            std::string low_chunk(
                chunk
            );

            low_chunk.push_back(
                '$'
            );

            low_chunk.append(
                chunk_length_
                    -
                    chunk.size()
                    -
                    1,
                'A'
            );


            /*
             * Highest possible completion.
             */
            std::string high_chunk(
                chunk
            );

            high_chunk.append(
                chunk_length_
                    -
                    chunk.size(),
                'T'
            );


            low =
                lower_bound(
                    low_chunk,
                    low
                );


            high =
                lower_bound(
                    high_chunk,
                    high
                );

        } else {

            low =
                lower_bound(
                    chunk,
                    low
                );


            high =
                lower_bound(
                    chunk,
                    high
                );
        }


        if (low >= high) {
            return
                Interval{
                    low,
                    low
                };
        }
    }


    return
        Interval{
            low,
            high
        };
}


std::vector<std::uint64_t>
IPBWTIndex::locate(
    const Interval& interval
) const {
    if (
        interval.begin >
            interval.end
        ||
        interval.end >
            suffix_array_.size()
    ) {
        throw std::out_of_range(
            "Invalid IP-BWT interval."
        );
    }


    std::vector<std::uint64_t>
        positions;


    positions.reserve(
        static_cast<std::size_t>(
            interval.size()
        )
    );


    for (
        std::uint64_t row =
            interval.begin;
        row <
            interval.end;
        ++row
    ) {
        const std::uint32_t start =
            suffix_array_.at(
                static_cast<std::size_t>(
                    row
                )
            );


        /*
         * Sentinel-only suffix is not a genomic
         * reference position.
         */
        if (
            static_cast<std::size_t>(
                start
            )
            <
            reference_length_
        ) {
            positions.push_back(
                static_cast<std::uint64_t>(
                    start
                )
            );
        }
    }


    std::sort(
        positions.begin(),
        positions.end()
    );


    return positions;
}


void IPBWTIndex::validate_reference(
    const std::string_view reference
) {
    if (reference.empty()) {
        throw std::invalid_argument(
            "IP-BWT reference cannot be empty."
        );
    }


    for (
        const char base :
        reference
    ) {
        switch (base) {

            case 'A':
            case 'C':
            case 'G':
            case 'N':
            case 'T':
                break;

            default:
                throw std::invalid_argument(
                    "IP-BWT reference must contain "
                    "only A/C/G/N/T."
                );
        }
    }
}


void IPBWTIndex::validate_query(
    const std::string_view query
) {
    if (query.empty()) {
        throw std::invalid_argument(
            "IP-BWT query cannot be empty."
        );
    }


    for (
        const char base :
        query
    ) {
        switch (base) {

            case 'A':
            case 'C':
            case 'G':
            case 'T':
                break;

            default:
                throw std::invalid_argument(
                    "IP-BWT query must contain "
                    "only A/C/G/T."
                );
        }
    }
}

}  // namespace primerpair
