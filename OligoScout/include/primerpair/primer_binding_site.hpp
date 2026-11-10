#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"

namespace primerpair {

struct PrimerBindingMismatch {
    /*
     * Zero-based position from primer 5-prime end.
     *
     * 0 = first 5-prime base.
     */
    std::size_t primer_position{0};

    char primer_base{'N'};
    char target_base{'N'};

    /*
     * 0 = terminal 3-prime base,
     * 1 = penultimate base, ...
     */
    std::size_t distance_from_3prime{0};

    bool operator==(
        const PrimerBindingMismatch&
    ) const = default;
};


struct PrimerBindingSite {
    PrimerIdentity primer{
        PrimerIdentity::Primer1
    };

    bool reverse_strand{false};

    std::uint64_t genomic_start{0};
    std::uint64_t genomic_end_exclusive{0};

    /*
     * Original primer, always 5' -> 3'.
     */
    std::string primer_sequence;

    /*
     * Reference genome sequence in forward
     * reference orientation.
     */
    std::string genomic_sequence;

    /*
     * Reference binding site oriented in the same
     * 5' -> 3' direction as the primer.
     *
     * left primer  -> genomic_sequence
     * right primer -> reverse complement
     */
    std::string binding_sequence;

    std::uint64_t mismatch_mask{0};

    std::vector<PrimerBindingMismatch>
        mismatches;

    [[nodiscard]]
    std::size_t mismatch_count() const noexcept {
        return mismatches.size();
    }
};


struct PrimerPairBindingSites {
    PrimerBindingSite left;
    PrimerBindingSite right;
};


[[nodiscard]]
std::string packed_reference_subsequence(
    const PackedReference& reference,
    std::uint64_t start,
    std::size_t length
);


[[nodiscard]]
std::string reverse_complement_binding_sequence(
    std::string_view sequence
);


[[nodiscard]]
PrimerPairBindingSites extract_pair_binding_sites(
    const PackedReference& reference,
    const PrimerPairHit& hit,
    std::string_view primer1,
    std::string_view primer2
);

}  // namespace primerpair
