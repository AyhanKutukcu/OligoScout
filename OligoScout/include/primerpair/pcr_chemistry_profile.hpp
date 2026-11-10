#pragma once

/*
 * PCR_CHEMISTRY_PROFILE_V1
 *
 * OligoScout
 *
 * Source-backed PCR chemistry / annealing
 * guidance metadata.
 *
 * IMPORTANT ARCHITECTURE:
 *
 * PCR chemistry is separate from:
 *
 * - search strategy
 * - Sensitive search mode
 * - biological risk scoring weights
 * - PCR probability
 *
 * No search hits are modified here.
 *
 * No biological score coefficients are
 * introduced here.
 */

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace primerpair {


enum class PcrChemistryKind : std::uint8_t {
    Custom = 0,
    StandardTaq = 1,
    Q5HighFidelity = 2,
    PhusionHighFidelity = 3
};


enum class PrimerTmBasis : std::uint8_t {

    /*
     * Generic calculated primer Tm.
     *
     * Appropriate only where the profile does not
     * explicitly require a chemistry-specific Tm.
     */
    GenericCalculated = 0,


    /*
     * Tm obtained using the chemistry-specific
     * manufacturer method / calculator / buffer
     * assumptions.
     */
    ChemistrySpecific = 1
};


enum class AnnealingTemperatureRule
    : std::uint8_t {

    UserSpecified = 0,

    /*
     * Ta = lower primer Tm - offset.
     */
    LowerTmMinusOffset = 1,

    /*
     * Ta = lower primer Tm + offset.
     */
    LowerTmPlusOffset = 2,

    /*
     * Phusion source-backed length-aware rule:
     *
     * >20 nt:
     *     Ta = lower primer Tm + offset
     *
     * <20 nt:
     *     Ta = lower primer Tm
     *
     * exactly 20 nt:
     *     not automatically inferred in V1.
     */
    PhusionLengthAware = 3
};


enum class ChemistryEvidenceSource
    : std::uint8_t {

    UserDefined = 0,

    ThermoFisherStandardPcrProtocol = 1,

    NebQ5Protocol = 2,

    NebPhusionProtocol = 3
};


struct PcrChemistryProfile {

    PcrChemistryKind kind{
        PcrChemistryKind::Custom
    };

    std::string_view name{
        "custom"
    };

    AnnealingTemperatureRule
        annealing_rule{
            AnnealingTemperatureRule::
                UserSpecified
        };

    ChemistryEvidenceSource
        evidence_source{
            ChemistryEvidenceSource::
                UserDefined
        };


    /*
     * Magnitude applied by the selected rule.
     */
    double annealing_offset_celsius{
        0.0
    };


    /*
     * Q5 / Phusion recommendations depend on
     * chemistry-specific Tm calculations.
     *
     * V1 therefore refuses to silently substitute
     * a generic Primer3 Tm for such profiles.
     */
    bool requires_chemistry_specific_tm{
        false
    };


    /*
     * True only for built-in profiles whose
     * annealing rule is based on an external
     * manufacturer protocol.
     */
    bool source_backed{
        false
    };


    /*
     * Must remain false in this milestone.
     */
    bool biological_risk_weights_calibrated{
        false
    };

    bool pcr_probability_calibrated{
        false
    };


    bool operator==(
        const PcrChemistryProfile&
    ) const = default;
};


[[nodiscard]]
PcrChemistryProfile
custom_pcr_chemistry_profile() noexcept;


[[nodiscard]]
PcrChemistryProfile
standard_taq_chemistry_profile() noexcept;


[[nodiscard]]
PcrChemistryProfile
q5_high_fidelity_chemistry_profile() noexcept;


[[nodiscard]]
PcrChemistryProfile
phusion_high_fidelity_chemistry_profile() noexcept;


/*
 * Calculate source-guided annealing temperature.
 *
 * lower_primer_tm_celsius:
 *     Tm of the lower-Tm primer.
 *
 * lower_primer_length:
 *     length of that primer.
 *
 * supplied_tm_basis:
 *     explicitly identifies whether the provided
 *     Tm is generic or chemistry-specific.
 *
 * Throws when automatic calculation would make
 * an unsupported assumption.
 */
[[nodiscard]]
double recommended_annealing_temperature(
    const PcrChemistryProfile& profile,
    double lower_primer_tm_celsius,
    std::size_t lower_primer_length,
    PrimerTmBasis supplied_tm_basis
);


}  // namespace primerpair
