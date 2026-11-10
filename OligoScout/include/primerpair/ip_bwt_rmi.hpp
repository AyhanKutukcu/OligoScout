#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"

namespace primerpair {


struct IPBWTRMILookupDiagnostics {
    std::uint64_t position{0};

    std::size_t selected_leaf{0};
    std::size_t predicted_position{0};

    std::size_t correction_radius{0};
    std::size_t correction_begin{0};
    std::size_t correction_end{0};

    bool used_global_fallback{false};
};


/*
 * Correctness-first Recursive Model Index for
 * IP-BWT lower-bound prediction.
 *
 * Architecture:
 *
 *     root linear regression
 *              ↓
 *       leaf selection
 *              ↓
 *     leaf linear regression
 *              ↓
 *      predicted position
 *              ↓
 *    exact local correction
 *              ↓
 * global binary fallback if necessary
 *
 * Prediction can be inaccurate.
 * Returned lower_bound cannot be approximate.
 */
class IPBWTRMI {
public:
    explicit IPBWTRMI(
        const IPBWTIndex& index,
        std::size_t requested_leaf_count = 16,
        std::size_t correction_margin = 8
    );

    [[nodiscard]]
    std::uint64_t lower_bound(
        std::string_view chunk,
        std::uint64_t row
    ) const;


    /*
     * Prediction-guided exact lower_bound.
     *
     * Uses the learned prediction as a starting
     * point and exponentially expands left/right
     * until the true lower-bound is bracketed.
     *
     * No global binary fallback is required.
     */
    [[nodiscard]]
    std::uint64_t lower_bound_galloping(
        std::string_view chunk,
        std::uint64_t row
    ) const;

    [[nodiscard]]
    Interval exact_search_galloping(
        std::string_view query
    ) const;


    [[nodiscard]]
    IPBWTRMILookupDiagnostics
    lower_bound_diagnostics(
        std::string_view chunk,
        std::uint64_t row
    ) const;


    /*
     * Diagnostic/ablation backend.
     *
     * Selects the exact leaf from sorted IP-BWT
     * leaf boundaries, then uses the SAME learned
     * leaf model and correction logic as RMI-v1.
     *
     * This isolates leaf-routing quality from
     * leaf-model quality.
     */
    [[nodiscard]]
    std::uint64_t
    lower_bound_boundary_routed(
        std::string_view chunk,
        std::uint64_t row
    ) const;

    [[nodiscard]]
    IPBWTRMILookupDiagnostics
    lower_bound_boundary_routed_diagnostics(
        std::string_view chunk,
        std::uint64_t row
    ) const;

    [[nodiscard]]
    Interval exact_search_boundary_routed(
        std::string_view query
    ) const;

    [[nodiscard]]
    Interval exact_search(
        std::string_view query
    ) const;

    [[nodiscard]]
    std::size_t leaf_count() const noexcept {
        return leaves_.size();
    }

    [[nodiscard]]
    std::size_t correction_margin() const noexcept {
        return correction_margin_;
    }

    [[nodiscard]]
    double
    mean_leaf_training_error() const noexcept {
        return mean_leaf_training_error_;
    }

    [[nodiscard]]
    std::size_t
    maximum_leaf_training_error() const noexcept {
        return maximum_leaf_training_error_;
    }

private:
    struct LinearModel {
        std::size_t begin{0};
        std::size_t end{0};

        std::uint64_t base_prefix{0};

        double slope{0.0};
        double intercept{0.0};

        std::size_t max_abs_error{0};
        double mean_abs_error{0.0};
    };

    const IPBWTIndex& index_;

    std::size_t correction_margin_{0};

    LinearModel root_;
    std::vector<LinearModel> leaves_;

    double mean_leaf_training_error_{0.0};
    std::size_t maximum_leaf_training_error_{0};

    [[nodiscard]]
    static bool key_less(
        std::uint64_t lhs_prefix,
        std::uint32_t lhs_row,
        std::uint64_t rhs_prefix,
        std::uint64_t rhs_row
    ) noexcept;

    [[nodiscard]]
    double model_input(
        std::uint64_t prefix,
        std::uint64_t row,
        std::uint64_t base_prefix
    ) const noexcept;

    [[nodiscard]]
    LinearModel train_model(
        std::size_t begin,
        std::size_t end
    ) const;

    [[nodiscard]]
    std::size_t predict(
        const LinearModel& model,
        std::uint64_t prefix,
        std::uint64_t row
    ) const noexcept;

    [[nodiscard]]
    std::size_t select_leaf(
        std::uint64_t prefix,
        std::uint64_t row
    ) const noexcept;


    [[nodiscard]]
    std::size_t select_leaf_by_boundaries(
        std::uint64_t prefix,
        std::uint64_t row
    ) const noexcept;

    [[nodiscard]]
    std::uint64_t exact_lower_bound_range(
        std::uint64_t prefix,
        std::uint64_t row,
        std::size_t begin,
        std::size_t end
    ) const noexcept;


    [[nodiscard]]
    std::uint64_t exact_lower_bound_galloping(
        std::uint64_t prefix,
        std::uint64_t row,
        std::size_t predicted
    ) const noexcept;

    [[nodiscard]]
    bool correction_window_contains_answer(
        std::uint64_t prefix,
        std::uint64_t row,
        std::size_t begin,
        std::size_t end
    ) const noexcept;

    static void validate_query(
        std::string_view query
    );
};

}  // namespace primerpair
