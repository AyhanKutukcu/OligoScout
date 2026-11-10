/*
 * THERMODYNAMIC_BACKEND_V1
 * Test #70
 */

#include "primerpair/thermodynamic_backend.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

void expect(
    const bool condition,
    const std::string_view message
) {
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


/*
 * ntthal may legitimately report a very low,
 * but finite, melting temperature for a weak
 * partial interaction.
 *
 * The adapter preserves the raw Primer3 value.
 * Biological interpretation / thresholding belongs
 * to the later risk-scoring layer.
 */
bool finite_tm(
    const double value
) {
    return
        std::isfinite(
            value
        );
}


std::string reverse_complement(
    const std::string_view sequence
) {
    std::string result;

    result.reserve(
        sequence.size()
    );

    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {
        switch (*it) {
            case 'A':
                result.push_back('T');
                break;

            case 'C':
                result.push_back('G');
                break;

            case 'G':
                result.push_back('C');
                break;

            case 'T':
                result.push_back('A');
                break;

            default:
                throw std::runtime_error(
                    "Invalid test DNA."
                );
        }
    }

    return result;
}


std::string required_environment(
    const char* name
) {
    const char* value =
        std::getenv(name);

    if (
        value == nullptr ||
        *value == '\0'
    ) {
        throw std::runtime_error(
            std::string(
                "Missing environment variable: "
            ) +
            name
        );
    }

    return value;
}

}  // namespace


int main()
{
    try {

        using namespace primerpair;

        const std::filesystem::path ntthal =
            required_environment(
                "PRIMERPAIR_NTTHAL"
            );

        const std::filesystem::path oligotm =
            required_environment(
                "PRIMERPAIR_OLIGOTM"
            );

        const std::filesystem::path config =
            required_environment(
                "PRIMERPAIR_THERMO_CONFIG"
            );

        Primer3ThermodynamicBackend backend(
            ntthal,
            oligotm,
            config
        );

        expect(
            backend.available(),
            "Primer3 backend available"
        );

        const std::string primer =
            "GCGTACGATCGTACGCATGC";

        const std::string complement =
            reverse_complement(
                primer
            );

        const std::string weak_partner =
            "AAAAAAAAAAAAAAAAAAAA";


        ThermodynamicConditions conditions;

        expect(
            conditions.monovalent_mM == 50.0,
            "Default monovalent concentration"
        );

        expect(
            conditions.divalent_mM == 1.5,
            "Default divalent concentration"
        );

        expect(
            conditions.dntp_mM == 0.6,
            "Default dNTP concentration"
        );

        expect(
            conditions.dna_nM == 50.0,
            "Default DNA concentration"
        );


        const double oligo_tm =
            backend.oligo_tm(
                primer,
                conditions
            );

        expect(
            finite_tm(
                oligo_tm
            ),
            "oligotm returns finite Tm"
        );


        const double oligo_tm_repeat =
            backend.oligo_tm(
                primer,
                conditions
            );

        expect(
            std::abs(
                oligo_tm -
                oligo_tm_repeat
            ) < 1.0e-9,
            "oligotm deterministic"
        );


        const double strong_any =
            backend.duplex_tm(
                primer,
                complement,
                ThermodynamicAlignment::Any,
                conditions
            );

        const double weak_any =
            backend.duplex_tm(
                primer,
                weak_partner,
                ThermodynamicAlignment::Any,
                conditions
            );

        expect(
            finite_tm(
                strong_any
            ),
            "Strong heterodimer ANY Tm finite"
        );

        expect(
            strong_any > 0.0,
            "Complementary duplex has positive Tm"
        );

        expect(
            finite_tm(
                weak_any
            ),
            "Weak heterodimer ANY raw Tm finite"
        );

        expect(
            strong_any >
                weak_any,
            "Complementary duplex stronger than weak partner"
        );


        const double strong_end1 =
            backend.duplex_tm(
                primer,
                complement,
                ThermodynamicAlignment::End1,
                conditions
            );

        const double strong_end2 =
            backend.duplex_tm(
                primer,
                complement,
                ThermodynamicAlignment::End2,
                conditions
            );

        expect(
            finite_tm(
                strong_end1
            ),
            "END1 interaction Tm finite"
        );

        expect(
            finite_tm(
                strong_end2
            ),
            "END2 interaction Tm finite"
        );


        const double hairpin_tm =
            backend.hairpin_tm(
                primer,
                conditions
            );

        expect(
            finite_tm(
                hairpin_tm
            ),
            "Hairpin Tm finite"
        );


        const auto primer_profile =
            backend.profile_primer(
                primer,
                conditions
            );

        expect(
            std::abs(
                primer_profile
                    .oligo_tm_celsius -
                oligo_tm
            ) < 1.0e-9,
            "Primer profile preserves oligo Tm"
        );

        expect(
            finite_tm(
                primer_profile
                    .homodimer_any_tm_celsius
            ),
            "Primer profile homodimer ANY finite"
        );

        expect(
            finite_tm(
                primer_profile
                    .hairpin_tm_celsius
            ),
            "Primer profile hairpin finite"
        );


        const auto pair_profile =
            backend.profile_pair(
                primer,
                complement,
                conditions
            );

        expect(
            finite_tm(
                pair_profile
                    .left
                    .oligo_tm_celsius
            ),
            "Pair profile left primer Tm finite"
        );

        expect(
            finite_tm(
                pair_profile
                    .right
                    .oligo_tm_celsius
            ),
            "Pair profile right primer Tm finite"
        );

        expect(
            std::abs(
                pair_profile
                    .heterodimer_any_tm_celsius -
                strong_any
            ) < 1.0e-9,
            "Pair profile heterodimer ANY exact"
        );


        bool invalid_base_rejected = false;

        try {
            static_cast<void>(
                backend.oligo_tm(
                    "ACGTXCGT"
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            invalid_base_rejected = true;
        }

        expect(
            invalid_base_rejected,
            "Ambiguous/invalid thermodynamic DNA rejected"
        );


        bool invalid_conditions_rejected = false;

        try {
            ThermodynamicConditions invalid;
            invalid.dna_nM = 0.0;

            static_cast<void>(
                backend.oligo_tm(
                    primer,
                    invalid
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            invalid_conditions_rejected = true;
        }

        expect(
            invalid_conditions_rejected,
            "Invalid DNA concentration rejected"
        );


        bool missing_backend_rejected = false;

        try {
            Primer3ThermodynamicBackend bad(
                "/definitely/missing/ntthal",
                oligotm,
                config
            );

            static_cast<void>(
                bad.available()
            );

        } catch (
            const std::invalid_argument&
        ) {
            missing_backend_rejected = true;
        }

        expect(
            missing_backend_rejected,
            "Missing backend executable rejected"
        );


        std::cout
            << "oligo_tm\t"
            << oligo_tm
            << '\n';

        std::cout
            << "strong_duplex_tm\t"
            << strong_any
            << '\n';

        std::cout
            << "weak_duplex_tm\t"
            << weak_any
            << '\n';

        std::cout
            << "hairpin_tm\t"
            << hairpin_tm
            << '\n';


        std::cout
            << "THERMODYNAMIC_BACKEND_AVAILABLE\tYES\n";

        std::cout
            << "THERMODYNAMIC_OLIGOTM_VALID\tYES\n";

        std::cout
            << "THERMODYNAMIC_NTTHAL_DUPLEX_VALID\tYES\n";

        std::cout
            << "THERMODYNAMIC_NTTHAL_HAIRPIN_VALID\tYES\n";

        std::cout
            << "THERMODYNAMIC_PROFILE_PAIR_VALID\tYES\n";

        std::cout
            << "THERMODYNAMIC_BACKEND_DETERMINISTIC\tYES\n";

        std::cout
            << "THERMODYNAMIC_BACKEND_DELTA_G\tNO\n";

        std::cout
            << "THERMODYNAMIC_BACKEND_PCR_PROBABILITY\tNO\n";

        std::cout
            << "THERMODYNAMIC_BACKEND_V1_COMPLETE\tYES\n";

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
