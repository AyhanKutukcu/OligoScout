/*
 * PCR_CHEMISTRY_PROFILE_V1
 */

#include "primerpair/pcr_chemistry_profile.hpp"

#include <cmath>
#include <stdexcept>

namespace primerpair {


PcrChemistryProfile
custom_pcr_chemistry_profile() noexcept
{
    PcrChemistryProfile output;

    output.kind =
        PcrChemistryKind::Custom;

    output.name =
        "custom";

    output.annealing_rule =
        AnnealingTemperatureRule::
            UserSpecified;

    output.evidence_source =
        ChemistryEvidenceSource::
            UserDefined;

    output.annealing_offset_celsius =
        0.0;

    output.requires_chemistry_specific_tm =
        false;

    output.source_backed =
        false;

    output.biological_risk_weights_calibrated =
        false;

    output.pcr_probability_calibrated =
        false;

    return output;
}


PcrChemistryProfile
standard_taq_chemistry_profile() noexcept
{
    PcrChemistryProfile output;

    output.kind =
        PcrChemistryKind::StandardTaq;

    output.name =
        "standard-taq";

    output.annealing_rule =
        AnnealingTemperatureRule::
            LowerTmMinusOffset;

    output.evidence_source =
        ChemistryEvidenceSource::
            ThermoFisherStandardPcrProtocol;

    /*
     * Source-backed Standard PCR guidance:
     *
     * Ta = lower primer Tm - 5 C.
     */
    output.annealing_offset_celsius =
        5.0;

    output.requires_chemistry_specific_tm =
        false;

    output.source_backed =
        true;

    output.biological_risk_weights_calibrated =
        false;

    output.pcr_probability_calibrated =
        false;

    return output;
}


PcrChemistryProfile
q5_high_fidelity_chemistry_profile() noexcept
{
    PcrChemistryProfile output;

    output.kind =
        PcrChemistryKind::Q5HighFidelity;

    output.name =
        "q5-high-fidelity";

    output.annealing_rule =
        AnnealingTemperatureRule::
            LowerTmPlusOffset;

    output.evidence_source =
        ChemistryEvidenceSource::
            NebQ5Protocol;

    /*
     * NEB Q5 guidance:
     *
     * typically lower-primer Tm + 3 C.
     *
     * The Tm itself must be determined using
     * chemistry-appropriate assumptions.
     */
    output.annealing_offset_celsius =
        3.0;

    output.requires_chemistry_specific_tm =
        true;

    output.source_backed =
        true;

    output.biological_risk_weights_calibrated =
        false;

    output.pcr_probability_calibrated =
        false;

    return output;
}


PcrChemistryProfile
phusion_high_fidelity_chemistry_profile() noexcept
{
    PcrChemistryProfile output;

    output.kind =
        PcrChemistryKind::
            PhusionHighFidelity;

    output.name =
        "phusion-high-fidelity";

    output.annealing_rule =
        AnnealingTemperatureRule::
            PhusionLengthAware;

    output.evidence_source =
        ChemistryEvidenceSource::
            NebPhusionProtocol;

    /*
     * Source-backed rule for primers >20 nt:
     *
     * lower-primer Tm + 3 C.
     *
     * Primers <20 nt:
     *
     * anneal at approximately lower-primer Tm.
     */
    output.annealing_offset_celsius =
        3.0;

    output.requires_chemistry_specific_tm =
        true;

    output.source_backed =
        true;

    output.biological_risk_weights_calibrated =
        false;

    output.pcr_probability_calibrated =
        false;

    return output;
}


double recommended_annealing_temperature(
    const PcrChemistryProfile& profile,
    const double lower_primer_tm_celsius,
    const std::size_t lower_primer_length,
    const PrimerTmBasis supplied_tm_basis
)
{
    if (
        !std::isfinite(
            lower_primer_tm_celsius
        )
    ) {
        throw std::invalid_argument(
            "Lower-primer Tm must be finite."
        );
    }


    if (lower_primer_length == 0) {
        throw std::invalid_argument(
            "Lower-primer length must be > 0."
        );
    }


    if (
        !std::isfinite(
            profile.annealing_offset_celsius
        ) ||
        profile.annealing_offset_celsius <
            0.0
    ) {
        throw std::invalid_argument(
            "Annealing offset must be finite "
            "and non-negative."
        );
    }


    if (
        profile.requires_chemistry_specific_tm &&
        supplied_tm_basis !=
            PrimerTmBasis::ChemistrySpecific
    ) {
        throw std::invalid_argument(
            "Selected PCR chemistry requires "
            "a chemistry-specific primer Tm. "
            "Generic Primer3 Tm is not silently "
            "substituted."
        );
    }


    switch (
        profile.annealing_rule
    ) {

        case AnnealingTemperatureRule::
            UserSpecified:

            throw std::invalid_argument(
                "Custom PCR chemistry requires "
                "an explicitly supplied annealing "
                "temperature."
            );


        case AnnealingTemperatureRule::
            LowerTmMinusOffset:

            return
                lower_primer_tm_celsius -
                profile
                    .annealing_offset_celsius;


        case AnnealingTemperatureRule::
            LowerTmPlusOffset:

            return
                lower_primer_tm_celsius +
                profile
                    .annealing_offset_celsius;


        case AnnealingTemperatureRule::
            PhusionLengthAware:

            if (lower_primer_length > 20) {

                return
                    lower_primer_tm_celsius +
                    profile
                        .annealing_offset_celsius;
            }


            if (lower_primer_length < 20) {

                return
                    lower_primer_tm_celsius;
            }


            /*
             * Manufacturer wording explicitly
             * distinguishes >20 and <20 nt.
             *
             * V1 therefore refuses to invent an
             * automatic rule for exactly 20 nt.
             */
            throw std::invalid_argument(
                "Phusion V1 automatic annealing "
                "guidance does not infer a rule "
                "for an exactly 20-nt lower-Tm "
                "primer. Supply an explicit "
                "annealing temperature."
            );
    }


    throw std::logic_error(
        "Unknown annealing-temperature rule."
    );
}


}  // namespace primerpair
