/*
 * PCR_CHEMISTRY_PROFILE_V1
 * Test #74
 */

#include "primerpair/pcr_chemistry_profile.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <iostream>

namespace {


void expect(
    const bool condition,
    const std::string_view message
)
{
    if (!condition) {
        throw std::runtime_error(
            std::string(message)
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}


bool almost_equal(
    const double lhs,
    const double rhs,
    const double tolerance = 1.0e-12
)
{
    return
        std::abs(
            lhs -
            rhs
        ) <=
        tolerance;
}


}  // namespace


int main()
{
    try {

        using namespace primerpair;


        const auto custom =
            custom_pcr_chemistry_profile();

        const auto taq =
            standard_taq_chemistry_profile();

        const auto q5 =
            q5_high_fidelity_chemistry_profile();

        const auto phusion =
            phusion_high_fidelity_chemistry_profile();


        expect(
            custom.kind ==
                PcrChemistryKind::Custom,
            "Custom chemistry identity valid"
        );


        expect(
            taq.kind ==
                PcrChemistryKind::StandardTaq,
            "Standard Taq chemistry identity valid"
        );


        expect(
            q5.kind ==
                PcrChemistryKind::Q5HighFidelity,
            "Q5 chemistry identity valid"
        );


        expect(
            phusion.kind ==
                PcrChemistryKind::
                    PhusionHighFidelity,
            "Phusion chemistry identity valid"
        );


        expect(
            !custom.source_backed,
            "Custom profile is user-defined"
        );


        expect(
            taq.source_backed &&
            q5.source_backed &&
            phusion.source_backed,
            "Built-in chemistry guidance source-backed"
        );


        expect(
            !taq
                .requires_chemistry_specific_tm,
            "Standard Taq accepts generic calculated Tm"
        );


        expect(
            q5
                .requires_chemistry_specific_tm &&
            phusion
                .requires_chemistry_specific_tm,
            "High-fidelity profiles require chemistry-specific Tm"
        );


        expect(
            !taq
                .biological_risk_weights_calibrated &&
            !q5
                .biological_risk_weights_calibrated &&
            !phusion
                .biological_risk_weights_calibrated,
            "Chemistry profiles do not claim calibrated risk weights"
        );


        expect(
            !taq.pcr_probability_calibrated &&
            !q5.pcr_probability_calibrated &&
            !phusion.pcr_probability_calibrated,
            "Chemistry profiles do not claim PCR probability"
        );


        /*
         * Standard Taq:
         *
         * lower primer Tm 60 -> annealing 55.
         */
        const double taq_ta =
            recommended_annealing_temperature(
                taq,
                60.0,
                20,
                PrimerTmBasis::
                    GenericCalculated
            );


        expect(
            almost_equal(
                taq_ta,
                55.0
            ),
            "Standard Taq Tm-minus-5 rule exact"
        );


        /*
         * Q5:
         *
         * chemistry-specific lower primer
         * Tm 60 -> annealing 63.
         */
        const double q5_ta =
            recommended_annealing_temperature(
                q5,
                60.0,
                24,
                PrimerTmBasis::
                    ChemistrySpecific
            );


        expect(
            almost_equal(
                q5_ta,
                63.0
            ),
            "Q5 chemistry-specific Tm-plus-3 rule exact"
        );


        bool q5_generic_rejected =
            false;


        try {

            static_cast<void>(
                recommended_annealing_temperature(
                    q5,
                    60.0,
                    24,
                    PrimerTmBasis::
                        GenericCalculated
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            q5_generic_rejected =
                true;
        }


        expect(
            q5_generic_rejected,
            "Q5 generic Primer3 Tm substitution rejected"
        );


        /*
         * Phusion >20 nt:
         * 60 -> 63.
         */
        const double phusion_long_ta =
            recommended_annealing_temperature(
                phusion,
                60.0,
                24,
                PrimerTmBasis::
                    ChemistrySpecific
            );


        expect(
            almost_equal(
                phusion_long_ta,
                63.0
            ),
            "Phusion >20-nt Tm-plus-3 rule exact"
        );


        /*
         * Phusion <20 nt:
         * use lower primer Tm.
         */
        const double phusion_short_ta =
            recommended_annealing_temperature(
                phusion,
                60.0,
                19,
                PrimerTmBasis::
                    ChemistrySpecific
            );


        expect(
            almost_equal(
                phusion_short_ta,
                60.0
            ),
            "Phusion <20-nt Tm rule exact"
        );


        bool phusion_20_rejected =
            false;


        try {

            static_cast<void>(
                recommended_annealing_temperature(
                    phusion,
                    60.0,
                    20,
                    PrimerTmBasis::
                        ChemistrySpecific
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            phusion_20_rejected =
                true;
        }


        expect(
            phusion_20_rejected,
            "Phusion exactly-20-nt unsupported inference rejected"
        );


        bool phusion_generic_rejected =
            false;


        try {

            static_cast<void>(
                recommended_annealing_temperature(
                    phusion,
                    60.0,
                    24,
                    PrimerTmBasis::
                        GenericCalculated
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            phusion_generic_rejected =
                true;
        }


        expect(
            phusion_generic_rejected,
            "Phusion generic Primer3 Tm substitution rejected"
        );


        bool custom_auto_rejected =
            false;


        try {

            static_cast<void>(
                recommended_annealing_temperature(
                    custom,
                    60.0,
                    20,
                    PrimerTmBasis::
                        GenericCalculated
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            custom_auto_rejected =
                true;
        }


        expect(
            custom_auto_rejected,
            "Custom chemistry requires explicit annealing temperature"
        );


        bool zero_length_rejected =
            false;


        try {

            static_cast<void>(
                recommended_annealing_temperature(
                    taq,
                    60.0,
                    0,
                    PrimerTmBasis::
                        GenericCalculated
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            zero_length_rejected =
                true;
        }


        expect(
            zero_length_rejected,
            "Zero primer length rejected"
        );


        bool nonfinite_tm_rejected =
            false;


        try {

            static_cast<void>(
                recommended_annealing_temperature(
                    taq,
                    std::numeric_limits<
                        double
                    >::quiet_NaN(),
                    20,
                    PrimerTmBasis::
                        GenericCalculated
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            nonfinite_tm_rejected =
                true;
        }


        expect(
            nonfinite_tm_rejected,
            "Non-finite primer Tm rejected"
        );


        /*
         * Critical architecture:
         *
         * Sensitive is intentionally absent from
         * PcrChemistryKind because it belongs to
         * search behavior rather than polymerase
         * chemistry.
         */
        expect(
            static_cast<int>(
                PcrChemistryKind::
                    PhusionHighFidelity
            ) == 3,
            "PCR chemistry enum remains chemistry-only"
        );


        std::cout
            << "taq_lower_tm_c\t"
            << 60.0
            << '\n';

        std::cout
            << "taq_annealing_c\t"
            << taq_ta
            << '\n';

        std::cout
            << "q5_lower_tm_c\t"
            << 60.0
            << '\n';

        std::cout
            << "q5_annealing_c\t"
            << q5_ta
            << '\n';

        std::cout
            << "phusion_long_annealing_c\t"
            << phusion_long_ta
            << '\n';

        std::cout
            << "phusion_short_annealing_c\t"
            << phusion_short_ta
            << '\n';


        std::cout
            << "PCR_CHEMISTRY_STANDARD_TAQ_VALID\tYES\n";

        std::cout
            << "PCR_CHEMISTRY_Q5_VALID\tYES\n";

        std::cout
            << "PCR_CHEMISTRY_PHUSION_VALID\tYES\n";

        std::cout
            << "PCR_CHEMISTRY_TM_BASIS_GUARD\tYES\n";

        std::cout
            << "PCR_CHEMISTRY_SOURCE_BACKED\tYES\n";

        std::cout
            << "PCR_CHEMISTRY_SENSITIVE_SEARCH_SEPARATE\tYES\n";

        std::cout
            << "PCR_CHEMISTRY_RISK_WEIGHTS\tNO\n";

        std::cout
            << "PCR_CHEMISTRY_PCR_PROBABILITY\tNO\n";

        std::cout
            << "PCR_CHEMISTRY_SEARCH_CHANGE\tNO\n";

        std::cout
            << "PCR_CHEMISTRY_PROFILE_V1_COMPLETE\tYES\n";

        std::cout
            << "ALL_CHECKS\tYES\n";


        return 0;

    } catch (
        const std::exception& e
    ) {

        std::cerr
            << "TEST_FAILURE\t"
            << e.what()
            << '\n';

        return 1;
    }
}
