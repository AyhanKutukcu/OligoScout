#pragma once

/*
 * THERMODYNAMIC_BACKEND_V1
 *
 * External Primer3 thermodynamic adapter.
 *
 * Scientific scope:
 * - oligo Tm via oligotm
 * - duplex interaction Tm via ntthal
 * - hairpin Tm via ntthal
 * - homodimer / heterodimer interaction Tm
 *
 * This layer does NOT yet represent:
 * - PCR probability
 * - calibrated off-target risk
 * - delta-G scoring
 * - mismatch-specific thermodynamic calibration
 */

#include <cstddef>
#include <filesystem>
#include <string_view>

namespace primerpair {

struct ThermodynamicConditions {
    double monovalent_mM{50.0};
    double divalent_mM{1.5};
    double dntp_mM{0.6};
    double dna_nM{50.0};

    double temperature_celsius{37.0};

    std::size_t max_loop{30};

    /*
     * oligotm:
     * 1 = SantaLucia 1998 nearest-neighbor parameters.
     */
    int tm_method{1};

    /*
     * oligotm:
     * 1 = SantaLucia 1998 salt correction.
     */
    int salt_correction{1};

    bool operator==(
        const ThermodynamicConditions&
    ) const = default;
};


enum class ThermodynamicAlignment {
    Any,
    End1,
    End2
};


struct PrimerThermodynamicProfile {
    double oligo_tm_celsius{0.0};

    double hairpin_tm_celsius{0.0};

    double homodimer_any_tm_celsius{0.0};
    double homodimer_end1_tm_celsius{0.0};
    double homodimer_end2_tm_celsius{0.0};

    bool operator==(
        const PrimerThermodynamicProfile&
    ) const = default;
};


struct PrimerPairThermodynamicProfile {
    PrimerThermodynamicProfile left{};
    PrimerThermodynamicProfile right{};

    double heterodimer_any_tm_celsius{0.0};
    double heterodimer_end1_tm_celsius{0.0};
    double heterodimer_end2_tm_celsius{0.0};

    bool operator==(
        const PrimerPairThermodynamicProfile&
    ) const = default;
};


class Primer3ThermodynamicBackend {
public:
    Primer3ThermodynamicBackend(
        std::filesystem::path ntthal_path,
        std::filesystem::path oligotm_path,
        std::filesystem::path
            thermodynamic_parameters_path
    );

    [[nodiscard]]
    bool available() const noexcept;


    [[nodiscard]]
    double oligo_tm(
        std::string_view sequence,
        const ThermodynamicConditions& conditions = {}
    ) const;


    [[nodiscard]]
    double duplex_tm(
        std::string_view sequence1,
        std::string_view sequence2,
        ThermodynamicAlignment alignment =
            ThermodynamicAlignment::Any,
        const ThermodynamicConditions& conditions = {}
    ) const;


    [[nodiscard]]
    double hairpin_tm(
        std::string_view sequence,
        const ThermodynamicConditions& conditions = {}
    ) const;


    [[nodiscard]]
    PrimerThermodynamicProfile
    profile_primer(
        std::string_view sequence,
        const ThermodynamicConditions& conditions = {}
    ) const;


    [[nodiscard]]
    PrimerPairThermodynamicProfile
    profile_pair(
        std::string_view left_primer,
        std::string_view right_primer,
        const ThermodynamicConditions& conditions = {}
    ) const;


    [[nodiscard]]
    const std::filesystem::path&
    ntthal_path() const noexcept {
        return ntthal_path_;
    }


    [[nodiscard]]
    const std::filesystem::path&
    oligotm_path() const noexcept {
        return oligotm_path_;
    }


    [[nodiscard]]
    const std::filesystem::path&
    thermodynamic_parameters_path() const noexcept {
        return thermodynamic_parameters_path_;
    }


private:
    std::filesystem::path ntthal_path_;
    std::filesystem::path oligotm_path_;

    std::filesystem::path
        thermodynamic_parameters_path_;
};

}  // namespace primerpair
