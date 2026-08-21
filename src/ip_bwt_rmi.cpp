#include "primerpair/ip_bwt_rmi.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace primerpair {

namespace {

constexpr double
    kRowScale =
        4294967296.0;  // 2^32

}  // namespace


IPBWTRMI::IPBWTRMI(
    const IPBWTIndex& index,
    const std::size_t requested_leaf_count,
    const std::size_t correction_margin
)
    : index_(
          index
      ),
      correction_margin_(
          correction_margin
      ) {

    const std::size_t n =
        static_cast<std::size_t>(
            index_.row_count()
        );

    if (n == 0) {
        throw std::invalid_argument(
            "Cannot train RMI on empty IP-BWT."
        );
    }

    if (requested_leaf_count == 0) {
        throw std::invalid_argument(
            "RMI leaf count cannot be zero."
        );
    }


    /*
     * Root model sees the complete IP-BWT.
     */
    root_ =
        train_model(
            0,
            n
        );


    const std::size_t actual_leaf_count =
        std::min(
            requested_leaf_count,
            n
        );


    leaves_.reserve(
        actual_leaf_count
    );


    /*
     * Equal-rank partitions for RMI v1.
     *
     * Later we can replace this with LISA-style
     * error-bounded bottom-up partitioning without
     * changing the public lookup semantics.
     */
    for (
        std::size_t leaf = 0;
        leaf < actual_leaf_count;
        ++leaf
    ) {
        const std::size_t begin =
            (
                leaf *
                n
            )
            /
            actual_leaf_count;

        const std::size_t end =
            (
                (
                    leaf + 1
                )
                *
                n
            )
            /
            actual_leaf_count;


        leaves_.push_back(
            train_model(
                begin,
                end
            )
        );
    }


    double error_sum = 0.0;

    maximum_leaf_training_error_ = 0;


    for (
        const auto& leaf :
        leaves_
    ) {
        error_sum +=
            leaf.mean_abs_error;

        maximum_leaf_training_error_ =
            std::max(
                maximum_leaf_training_error_,
                leaf.max_abs_error
            );
    }


    mean_leaf_training_error_ =
        error_sum /
        static_cast<double>(
            leaves_.size()
        );
}


bool IPBWTRMI::key_less(
    const std::uint64_t lhs_prefix,
    const std::uint32_t lhs_row,
    const std::uint64_t rhs_prefix,
    const std::uint64_t rhs_row
) noexcept {

    if (
        lhs_prefix <
        rhs_prefix
    ) {
        return true;
    }

    if (
        lhs_prefix >
        rhs_prefix
    ) {
        return false;
    }

    return
        static_cast<std::uint64_t>(
            lhs_row
        )
        <
        rhs_row;
}


double IPBWTRMI::model_input(
    const std::uint64_t prefix,
    const std::uint64_t row,
    const std::uint64_t base_prefix
) const noexcept {

    double prefix_delta = 0.0;

    if (
        prefix >=
        base_prefix
    ) {
        prefix_delta =
            static_cast<double>(
                prefix -
                base_prefix
            );

    } else {

        prefix_delta =
            -
            static_cast<double>(
                base_prefix -
                prefix
            );
    }


    const double row_fraction =
        static_cast<double>(
            row
        )
        /
        kRowScale;


    return
        prefix_delta +
        row_fraction;
}


IPBWTRMI::LinearModel
IPBWTRMI::train_model(
    const std::size_t begin,
    const std::size_t end
) const {

    if (
        begin >= end ||
        end >
            index_.prefix_codes_.size()
    ) {
        throw std::invalid_argument(
            "Invalid RMI training range."
        );
    }


    LinearModel model;

    model.begin =
        begin;

    model.end =
        end;

    model.base_prefix =
        index_.prefix_codes_.at(
            begin
        );


    const std::size_t count =
        end -
        begin;


    double mean_x = 0.0;
    double mean_y = 0.0;


    for (
        std::size_t position = begin;
        position < end;
        ++position
    ) {
        const double x =
            model_input(
                index_.prefix_codes_.at(
                    position
                ),
                index_.paired_rows_.at(
                    position
                ),
                model.base_prefix
            );

        const double y =
            static_cast<double>(
                position
            );

        mean_x += x;
        mean_y += y;
    }


    mean_x /=
        static_cast<double>(
            count
        );

    mean_y /=
        static_cast<double>(
            count
        );


    double numerator = 0.0;
    double denominator = 0.0;


    for (
        std::size_t position = begin;
        position < end;
        ++position
    ) {
        const double x =
            model_input(
                index_.prefix_codes_.at(
                    position
                ),
                index_.paired_rows_.at(
                    position
                ),
                model.base_prefix
            );

        const double y =
            static_cast<double>(
                position
            );


        const double dx =
            x -
            mean_x;

        const double dy =
            y -
            mean_y;


        numerator +=
            dx *
            dy;

        denominator +=
            dx *
            dx;
    }


    if (
        denominator >
        std::numeric_limits<
            double
        >::epsilon()
    ) {
        model.slope =
            numerator /
            denominator;

    } else {

        model.slope =
            0.0;
    }


    model.intercept =
        mean_y -
        model.slope *
        mean_x;


    double absolute_error_sum =
        0.0;

    std::size_t max_error = 0;


    for (
        std::size_t position = begin;
        position < end;
        ++position
    ) {
        const std::size_t predicted =
            predict(
                model,
                index_.prefix_codes_.at(
                    position
                ),
                index_.paired_rows_.at(
                    position
                )
            );


        const std::size_t error =
            (
                predicted >
                position
            )
            ?
            predicted -
                position
            :
            position -
                predicted;


        max_error =
            std::max(
                max_error,
                error
            );


        absolute_error_sum +=
            static_cast<double>(
                error
            );
    }


    model.max_abs_error =
        max_error;


    model.mean_abs_error =
        absolute_error_sum /
        static_cast<double>(
            count
        );


    return model;
}


std::size_t IPBWTRMI::predict(
    const LinearModel& model,
    const std::uint64_t prefix,
    const std::uint64_t row
) const noexcept {

    const std::size_t n =
        index_.prefix_codes_.size();


    const double x =
        model_input(
            prefix,
            row,
            model.base_prefix
        );


    const double estimate =
        model.slope *
        x
        +
        model.intercept;


    if (
        !std::isfinite(
            estimate
        )
    ) {
        return model.begin;
    }


    if (
        estimate <=
        0.0
    ) {
        return 0;
    }


    const double upper =
        static_cast<double>(
            n - 1
        );


    if (
        estimate >=
        upper
    ) {
        return
            n - 1;
    }


    return
        static_cast<std::size_t>(
            estimate
        );
}


std::size_t IPBWTRMI::select_leaf(
    const std::uint64_t prefix,
    const std::uint64_t row
) const noexcept {

    const std::size_t n =
        index_.prefix_codes_.size();


    const std::size_t root_prediction =
        predict(
            root_,
            prefix,
            row
        );


    std::size_t leaf =
        (
            root_prediction *
            leaves_.size()
        )
        /
        n;


    if (
        leaf >=
        leaves_.size()
    ) {
        leaf =
            leaves_.size() - 1;
    }


    return leaf;
}


std::size_t
IPBWTRMI::select_leaf_by_boundaries(
    const std::uint64_t prefix,
    const std::uint64_t row
) const noexcept {

    /*
     * Find the first leaf whose LAST compound key
     * is >= the query key.
     *
     * Leaves partition the sorted IP-BWT into
     * contiguous rank ranges.
     */
    std::size_t low = 0;

    std::size_t high =
        leaves_.size();


    while (
        low <
        high
    ) {
        const std::size_t middle =
            low +
            (
                high -
                low
            )
            /
            2;


        const LinearModel& leaf =
            leaves_[
                middle
            ];


        const std::size_t last_position =
            leaf.end -
            1;


        const bool leaf_entirely_before_query =
            key_less(
                index_.prefix_codes_[
                    last_position
                ],
                index_.paired_rows_[
                    last_position
                ],
                prefix,
                row
            );


        if (
            leaf_entirely_before_query
        ) {
            low =
                middle + 1;

        } else {
            high =
                middle;
        }
    }


    if (
        low >=
        leaves_.size()
    ) {
        return
            leaves_.size() - 1;
    }


    return low;
}


std::uint64_t
IPBWTRMI::exact_lower_bound_range(
    const std::uint64_t prefix,
    const std::uint64_t row,
    std::size_t begin,
    std::size_t end
) const noexcept {

    while (
        begin <
        end
    ) {
        const std::size_t middle =
            begin +
            (
                end -
                begin
            )
            /
            2;


        const bool less =
            key_less(
                index_.prefix_codes_[middle],
                index_.paired_rows_[middle],
                prefix,
                row
            );


        if (less) {
            begin =
                middle + 1;
        } else {
            end =
                middle;
        }
    }


    return
        static_cast<std::uint64_t>(
            begin
        );
}


std::uint64_t
IPBWTRMI::exact_lower_bound_galloping(
    const std::uint64_t prefix,
    const std::uint64_t row,
    const std::size_t predicted
) const noexcept {

    const std::size_t n =
        index_.prefix_codes_.size();


    if (n == 0) {
        return 0;
    }


    const std::size_t position =
        std::min(
            predicted,
            n - 1
        );


    const auto entry_less =
        [
            this,
            prefix,
            row
        ](
            const std::size_t index
        ) noexcept {

            return
                key_less(
                    index_.prefix_codes_[index],
                    index_.paired_rows_[index],
                    prefix,
                    row
                );
        };


    /*
     * ====================================================
     * CASE 1
     *
     * entry[position] < query key
     *
     * Therefore:
     *
     *     lower_bound > position
     *
     * Gallop right until an entry >= query key is found.
     * ====================================================
     */
    if (
        entry_less(
            position
        )
    ) {
        std::size_t lower =
            position + 1;


        if (
            lower >= n
        ) {
            return
                static_cast<std::uint64_t>(
                    n
                );
        }


        std::size_t step = 1;


        while (true) {

            /*
             * Avoid overflow and clamp to the final row.
             */
            if (
                step >
                n - 1 - position
            ) {
                return
                    exact_lower_bound_range(
                        prefix,
                        row,
                        lower,
                        n
                    );
            }


            const std::size_t probe =
                position +
                step;


            if (
                !entry_less(
                    probe
                )
            ) {
                return
                    exact_lower_bound_range(
                        prefix,
                        row,
                        lower,
                        probe + 1
                    );
            }


            lower =
                probe + 1;


            if (
                lower >= n
            ) {
                return
                    static_cast<std::uint64_t>(
                        n
                    );
            }


            if (
                step >
                n / 2
            ) {
                return
                    exact_lower_bound_range(
                        prefix,
                        row,
                        lower,
                        n
                    );
            }


            step *= 2;
        }
    }


    /*
     * ====================================================
     * CASE 2
     *
     * entry[position] >= query key
     *
     * Therefore:
     *
     *     lower_bound <= position
     *
     * Gallop left until an entry < query key is found.
     * ====================================================
     */

    std::size_t upper =
        position + 1;


    if (
        position == 0
    ) {
        return 0;
    }


    std::size_t step = 1;


    while (true) {

        const std::size_t probe =
            (
                step >
                position
            )
            ?
            0
            :
            position -
                step;


        if (
            entry_less(
                probe
            )
        ) {
            return
                exact_lower_bound_range(
                    prefix,
                    row,
                    probe + 1,
                    upper
                );
        }


        upper =
            probe + 1;


        if (
            probe == 0
        ) {
            return 0;
        }


        if (
            step >
            position / 2
        ) {
            step =
                position + 1;
        } else {
            step *= 2;
        }
    }
}


bool
IPBWTRMI::correction_window_contains_answer(
    const std::uint64_t prefix,
    const std::uint64_t row,
    const std::size_t begin,
    const std::size_t end
) const noexcept {

    const std::size_t n =
        index_.prefix_codes_.size();


    /*
     * lower_bound >= begin iff
     * the entry immediately before begin is < key.
     */
    const bool left_ok =
        (
            begin == 0
        )
        ||
        key_less(
            index_.prefix_codes_[
                begin - 1
            ],
            index_.paired_rows_[
                begin - 1
            ],
            prefix,
            row
        );


    /*
     * lower_bound <= end iff
     * end == n or entry[end] >= key.
     *
     * Note that 'end' is an exclusive search bound.
     */
    const bool right_ok =
        (
            end >= n
        )
        ||
        !key_less(
            index_.prefix_codes_[
                end
            ],
            index_.paired_rows_[
                end
            ],
            prefix,
            row
        );


    return
        left_ok &&
        right_ok;
}


std::uint64_t
IPBWTRMI::lower_bound_galloping(
    const std::string_view chunk,
    const std::uint64_t row
) const {

    if (
        chunk.size() !=
        index_.chunk_length_
    ) {
        throw std::invalid_argument(
            "Galloping RMI lookup chunk "
            "has wrong length."
        );
    }


    if (
        row >
        index_.prefix_codes_.size()
    ) {
        throw std::out_of_range(
            "Galloping RMI lookup row "
            "is out of range."
        );
    }


    const std::uint64_t prefix =
        index_.encode_chunk(
            chunk
        );


    const std::size_t leaf_index =
        select_leaf(
            prefix,
            row
        );


    const LinearModel& leaf =
        leaves_.at(
            leaf_index
        );


    const std::size_t predicted =
        predict(
            leaf,
            prefix,
            row
        );


    return
        exact_lower_bound_galloping(
            prefix,
            row,
            predicted
        );
}


std::uint64_t IPBWTRMI::lower_bound(
    const std::string_view chunk,
    const std::uint64_t row
) const {

    return
        lower_bound_diagnostics(
            chunk,
            row
        ).position;
}


IPBWTRMILookupDiagnostics
IPBWTRMI::lower_bound_diagnostics(
    const std::string_view chunk,
    const std::uint64_t row
) const {

    if (
        chunk.size() !=
        index_.chunk_length_
    ) {
        throw std::invalid_argument(
            "RMI lookup chunk has wrong length."
        );
    }


    if (
        row >
        index_.prefix_codes_.size()
    ) {
        throw std::out_of_range(
            "RMI lookup row is out of range."
        );
    }


    const std::uint64_t prefix =
        index_.encode_chunk(
            chunk
        );


    const std::size_t n =
        index_.prefix_codes_.size();


    const std::size_t leaf_index =
        select_leaf(
            prefix,
            row
        );


    const LinearModel& leaf =
        leaves_.at(
            leaf_index
        );


    const std::size_t predicted =
        predict(
            leaf,
            prefix,
            row
        );


    const std::size_t radius =
        leaf.max_abs_error
        +
        correction_margin_;


    const std::size_t local_begin =
        (
            predicted >
            radius
        )
        ?
        predicted -
            radius
        :
        0;


    const std::size_t local_end =
        std::min(
            n,
            predicted +
                radius +
                1
        );


    IPBWTRMILookupDiagnostics
        diagnostics;

    diagnostics.selected_leaf =
        leaf_index;

    diagnostics.predicted_position =
        predicted;

    diagnostics.correction_radius =
        radius;

    diagnostics.correction_begin =
        local_begin;

    diagnostics.correction_end =
        local_end;


    if (
        correction_window_contains_answer(
            prefix,
            row,
            local_begin,
            local_end
        )
    ) {
        diagnostics.position =
            exact_lower_bound_range(
                prefix,
                row,
                local_begin,
                local_end
            );

        diagnostics.used_global_fallback =
            false;

        return diagnostics;
    }


    diagnostics.position =
        exact_lower_bound_range(
            prefix,
            row,
            0,
            n
        );

    diagnostics.used_global_fallback =
        true;

    return diagnostics;
}


std::uint64_t
IPBWTRMI::lower_bound_boundary_routed(
    const std::string_view chunk,
    const std::uint64_t row
) const {

    return
        lower_bound_boundary_routed_diagnostics(
            chunk,
            row
        ).position;
}


IPBWTRMILookupDiagnostics
IPBWTRMI::lower_bound_boundary_routed_diagnostics(
    const std::string_view chunk,
    const std::uint64_t row
) const {

    if (
        chunk.size() !=
        index_.chunk_length_
    ) {
        throw std::invalid_argument(
            "Boundary-routed RMI lookup chunk "
            "has wrong length."
        );
    }


    if (
        row >
        index_.prefix_codes_.size()
    ) {
        throw std::out_of_range(
            "Boundary-routed RMI row "
            "is out of range."
        );
    }


    const std::uint64_t prefix =
        index_.encode_chunk(
            chunk
        );


    const std::size_t n =
        index_.prefix_codes_.size();


    const std::size_t leaf_index =
        select_leaf_by_boundaries(
            prefix,
            row
        );


    const LinearModel& leaf =
        leaves_.at(
            leaf_index
        );


    const std::size_t predicted =
        predict(
            leaf,
            prefix,
            row
        );


    const std::size_t radius =
        leaf.max_abs_error
        +
        correction_margin_;


    const std::size_t local_begin =
        (
            predicted >
            radius
        )
        ?
        predicted -
            radius
        :
        0;


    const std::size_t local_end =
        std::min(
            n,
            predicted +
                radius +
                1
        );


    IPBWTRMILookupDiagnostics
        diagnostics;


    diagnostics.selected_leaf =
        leaf_index;

    diagnostics.predicted_position =
        predicted;

    diagnostics.correction_radius =
        radius;

    diagnostics.correction_begin =
        local_begin;

    diagnostics.correction_end =
        local_end;


    if (
        correction_window_contains_answer(
            prefix,
            row,
            local_begin,
            local_end
        )
    ) {
        diagnostics.position =
            exact_lower_bound_range(
                prefix,
                row,
                local_begin,
                local_end
            );

        diagnostics.used_global_fallback =
            false;

        return diagnostics;
    }


    /*
     * Correctness remains absolute.
     */
    diagnostics.position =
        exact_lower_bound_range(
            prefix,
            row,
            0,
            n
        );

    diagnostics.used_global_fallback =
        true;


    return diagnostics;
}


Interval
IPBWTRMI::exact_search_boundary_routed(
    const std::string_view query
) const {

    validate_query(
        query
    );


    const std::uint64_t n =
        index_.row_count();


    std::uint64_t low = 0;
    std::uint64_t high = n;


    std::vector<std::string_view>
        chunks;


    for (
        std::size_t position = 0;
        position < query.size();
        position += index_.chunk_length_
    ) {
        const std::size_t remaining =
            query.size() -
            position;


        const std::size_t width =
            std::min(
                index_.chunk_length_,
                remaining
            );


        chunks.push_back(
            query.substr(
                position,
                width
            )
        );
    }


    for (
        auto it = chunks.rbegin();
        it != chunks.rend();
        ++it
    ) {
        const std::string_view chunk =
            *it;


        if (
            chunk.size() <
            index_.chunk_length_
        ) {
            std::string low_chunk(
                chunk
            );

            low_chunk.push_back(
                '$'
            );

            low_chunk.append(
                index_.chunk_length_
                    -
                    chunk.size()
                    -
                    1,
                'A'
            );


            std::string high_chunk(
                chunk
            );

            high_chunk.append(
                index_.chunk_length_
                    -
                    chunk.size(),
                'T'
            );


            low =
                lower_bound_boundary_routed(
                    low_chunk,
                    low
                );


            high =
                lower_bound_boundary_routed(
                    high_chunk,
                    high
                );

        } else {

            low =
                lower_bound_boundary_routed(
                    chunk,
                    low
                );


            high =
                lower_bound_boundary_routed(
                    chunk,
                    high
                );
        }


        if (
            low >=
            high
        ) {
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


Interval
IPBWTRMI::exact_search_galloping(
    const std::string_view query
) const {

    validate_query(
        query
    );


    const std::uint64_t n =
        index_.row_count();


    std::uint64_t low = 0;
    std::uint64_t high = n;


    std::vector<std::string_view>
        chunks;


    for (
        std::size_t position = 0;
        position < query.size();
        position += index_.chunk_length_
    ) {
        const std::size_t remaining =
            query.size() -
            position;


        const std::size_t width =
            std::min(
                index_.chunk_length_,
                remaining
            );


        chunks.push_back(
            query.substr(
                position,
                width
            )
        );
    }


    for (
        auto it = chunks.rbegin();
        it != chunks.rend();
        ++it
    ) {
        const std::string_view chunk =
            *it;


        if (
            chunk.size() <
            index_.chunk_length_
        ) {
            std::string low_chunk(
                chunk
            );


            low_chunk.push_back(
                '$'
            );


            low_chunk.append(
                index_.chunk_length_
                    -
                    chunk.size()
                    -
                    1,
                'A'
            );


            std::string high_chunk(
                chunk
            );


            high_chunk.append(
                index_.chunk_length_
                    -
                    chunk.size(),
                'T'
            );


            low =
                lower_bound_galloping(
                    low_chunk,
                    low
                );


            high =
                lower_bound_galloping(
                    high_chunk,
                    high
                );

        } else {

            low =
                lower_bound_galloping(
                    chunk,
                    low
                );


            high =
                lower_bound_galloping(
                    chunk,
                    high
                );
        }


        if (
            low >=
            high
        ) {
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


Interval IPBWTRMI::exact_search(
    const std::string_view query
) const {

    validate_query(
        query
    );


    const std::uint64_t n =
        index_.row_count();


    std::uint64_t low = 0;
    std::uint64_t high = n;


    std::vector<std::string_view>
        chunks;


    for (
        std::size_t position = 0;
        position < query.size();
        position += index_.chunk_length_
    ) {
        const std::size_t remaining =
            query.size() -
            position;


        const std::size_t width =
            std::min(
                index_.chunk_length_,
                remaining
            );


        chunks.push_back(
            query.substr(
                position,
                width
            )
        );
    }


    for (
        auto it = chunks.rbegin();
        it != chunks.rend();
        ++it
    ) {
        const std::string_view chunk =
            *it;


        if (
            chunk.size() <
            index_.chunk_length_
        ) {
            std::string low_chunk(
                chunk
            );


            low_chunk.push_back(
                '$'
            );


            low_chunk.append(
                index_.chunk_length_
                    -
                    chunk.size()
                    -
                    1,
                'A'
            );


            std::string high_chunk(
                chunk
            );


            high_chunk.append(
                index_.chunk_length_
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


        if (
            low >=
            high
        ) {
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


void IPBWTRMI::validate_query(
    const std::string_view query
) {
    if (
        query.empty()
    ) {
        throw std::invalid_argument(
            "RMI query cannot be empty."
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
                    "RMI query must contain "
                    "only A/C/G/T."
                );
        }
    }
}

}  // namespace primerpair
