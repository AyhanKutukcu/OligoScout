#pragma once

/*
 * Sensitive ("terminal rule") primer site search.
 *
 * The default OligoScout contract enumerates sites whose 3'-terminal
 * `anchor_length` bases match exactly. Primer-BLAST's default specificity
 * rule is weaker: a target is still considered amplifiable unless it carries
 * at least `terminal_min_to_ignore` (2) mismatches within the 3'-terminal
 * `terminal_window` (5) bases. This searcher enumerates EVERY site with
 *
 *     total mismatches            <= max_mismatches          (0..3)
 *     mismatches in last window   <= max_terminal_mismatches (default 1)
 *
 * on both strands, with no exact-anchor requirement. It walks the
 * bidirectional FM-index from the primer's 3' end (left-extension for the
 * plus strand, right-extension for the minus strand) so the terminal budget
 * prunes first, then verifies every located site against the packed
 * reference (N rejects) and recomputes the mismatch mask.
 *
 * The anchored contract is a subset: every anchored site has zero
 * mismatches in its last 5 nt. With terminal_window = anchor length and
 * max_terminal_mismatches = 0 the search reproduces the anchored contract.
 *
 * Primers may carry IUPAC codes (R Y S W K M B D H V N): a template base
 * matches when it is in the code's set. The reference itself is A/C/G/T/N
 * (N rejects), so no template ambiguity arises here.
 */

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {

struct TerminalRuleOptions {
    std::size_t max_mismatches{3};
    std::size_t terminal_window{5};
    std::size_t max_terminal_mismatches{1};
};

struct TerminalRuleStats {
    std::uint64_t states_visited{0};
    std::uint64_t located{0};
    std::uint64_t rejected_on_verify{0};
};

[[nodiscard]]
std::vector<OrientedPrimerSearchHit>
search_terminal_rule(
    const BidirectionalFMIndex& index,
    const PackedReference& reference,
    std::string_view primer,
    const TerminalRuleOptions& options,
    TerminalRuleStats* stats = nullptr
);

}  // namespace primerpair
