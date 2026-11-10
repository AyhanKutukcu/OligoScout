#pragma once

/*
 * BINDING_SITE_THERMODYNAMICS_V1
 *
 * Thermodynamic interpretation of an already
 * discovered PrimerBindingSite.
 *
 * Important:
 *
 * PrimerBindingSite::binding_sequence is stored
 * in the same biological 5' -> 3' orientation as
 * the original primer.
 *
 * ntthal requires the actual antiparallel partner.
 * Therefore:
 *
 * target_for_ntthal =
 *     reverse_complement(binding_sequence)
 *
 * Delta-Tm convention:
 *
 * delta_tm_celsius =
 *     perfect_match_tm_celsius
 *     -
 *     observed_binding_tm_celsius
 *
 * Positive delta-Tm therefore represents a loss
 * of duplex stability relative to the perfect
 * primer-target duplex.
 *
 * This feature is NOT:
 * - PCR probability
 * - calibrated off-target risk
 * - delta-G
 */

#include <cstddef>

#include "primerpair/primer_binding_site.hpp"
#include "primerpair/thermodynamic_backend.hpp"

namespace primerpair {


struct BindingSiteThermodynamicProfile {

    PrimerIdentity primer{
        PrimerIdentity::Primer1
    };

    bool reverse_strand{false};

    std::size_t primer_length{0};

    std::size_t mismatch_count{0};

    /*
     * Primer hybridized against its ideal
     * perfectly complementary target.
     */
    double perfect_match_tm_celsius{0.0};

    /*
     * Primer hybridized against the actual
     * genomic binding site.
     */
    double observed_binding_tm_celsius{0.0};

    /*
     * perfect - observed
     *
     * Positive:
     * observed site is thermodynamically weaker.
     *
     * Zero:
     * no Tm penalty relative to perfect match.
     *
     * Negative:
     * observed interaction is predicted stronger.
     * This is preserved rather than artificially
     * clipped.
     */
    double delta_tm_celsius{0.0};

    [[nodiscard]]
    bool destabilized() const noexcept
    {
        return
            delta_tm_celsius >
            0.0;
    }

    bool operator==(
        const BindingSiteThermodynamicProfile&
    ) const = default;
};


struct PrimerPairBindingThermodynamicProfile {

    BindingSiteThermodynamicProfile left{};
    BindingSiteThermodynamicProfile right{};

    double mean_delta_tm_celsius{0.0};

    double max_delta_tm_celsius{0.0};
};


[[nodiscard]]
BindingSiteThermodynamicProfile
profile_binding_site_thermodynamics(
    const PrimerBindingSite& site,
    const Primer3ThermodynamicBackend& backend,
    const ThermodynamicConditions& conditions = {}
);


[[nodiscard]]
PrimerPairBindingThermodynamicProfile
profile_pair_binding_thermodynamics(
    const PrimerPairBindingSites& sites,
    const Primer3ThermodynamicBackend& backend,
    const ThermodynamicConditions& conditions = {}
);


}  // namespace primerpair
