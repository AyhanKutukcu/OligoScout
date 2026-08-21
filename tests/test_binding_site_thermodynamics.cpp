/*
 * BINDING_SITE_THERMODYNAMICS_V1
 * Test #71
 */

#include "primerpair/binding_site_thermodynamics.hpp"

#include <bit>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

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


std::string required_environment(
    const char* name
)
{
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


char alternative_base(
    const char original,
    const std::size_t variant_index
)
{
    constexpr char bases[] = {
        'A',
        'C',
        'G',
        'T'
    };

    std::size_t seen = 0;

    for (const char base : bases) {

        if (base == original) {
            continue;
        }

        if (seen == variant_index) {
            return base;
        }

        ++seen;
    }

    throw std::runtime_error(
        "Could not select alternative base."
    );
}


primerpair::PrimerBindingSite
make_site(
    const std::string& primer,
    const std::string& binding,
    const bool reverse_strand,
    const primerpair::PrimerIdentity identity
)
{
    using namespace primerpair;

    if (
        primer.size() !=
        binding.size()
    ) {
        throw std::runtime_error(
            "Test site sequence lengths disagree."
        );
    }

    PrimerBindingSite site;

    site.primer =
        identity;

    site.reverse_strand =
        reverse_strand;

    site.genomic_start =
        1000;

    site.genomic_end_exclusive =
        site.genomic_start +
        primer.size();

    site.primer_sequence =
        primer;

    site.binding_sequence =
        binding;

    site.genomic_sequence =
        reverse_strand
            ?
            reverse_complement_binding_sequence(
                binding
            )
            :
            binding;

    site.mismatch_mask =
        0;

    for (
        std::size_t i = 0;
        i < primer.size();
        ++i
    ) {
        if (
            primer[i] ==
            binding[i]
        ) {
            continue;
        }

        site.mismatch_mask |=
            (
                std::uint64_t{1}
                <<
                i
            );

        PrimerBindingMismatch mismatch;

        mismatch.primer_position =
            i;

        mismatch.primer_base =
            primer[i];

        mismatch.target_base =
            binding[i];

        mismatch.distance_from_3prime =
            primer.size() -
            1 -
            i;

        site.mismatches.push_back(
            mismatch
        );
    }

    return site;
}


bool almost_equal(
    const double lhs,
    const double rhs,
    const double tolerance = 1.0e-6
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


        ThermodynamicConditions conditions;


        const std::string primer =
            "GCGTACGATCGTACGCATGC";


        const auto perfect_site =
            make_site(
                primer,
                primer,
                false,
                PrimerIdentity::Primer1
            );


        const auto perfect_profile =
            profile_binding_site_thermodynamics(
                perfect_site,
                backend,
                conditions
            );


        expect(
            perfect_profile.primer_length ==
                primer.size(),
            "Perfect site primer length preserved"
        );

        expect(
            perfect_profile.mismatch_count == 0,
            "Perfect site mismatch count zero"
        );

        expect(
            std::isfinite(
                perfect_profile
                    .perfect_match_tm_celsius
            ),
            "Perfect reference Tm finite"
        );

        expect(
            std::isfinite(
                perfect_profile
                    .observed_binding_tm_celsius
            ),
            "Perfect observed Tm finite"
        );

        expect(
            almost_equal(
                perfect_profile
                    .perfect_match_tm_celsius,
                perfect_profile
                    .observed_binding_tm_celsius
            ),
            "Perfect site observed Tm equals perfect Tm"
        );

        expect(
            almost_equal(
                perfect_profile
                    .delta_tm_celsius,
                0.0
            ),
            "Perfect site delta-Tm equals zero"
        );


        /*
         * Find a deterministic single mismatch outside
         * the strict terminal 12-base anchor that lowers
         * duplex Tm.
         *
         * This avoids hard-coding a thermodynamic claim
         * about one particular mismatch substitution.
         */
        const std::size_t strict_anchor_length =
            12;

        const std::size_t mutable_prefix =
            primer.size() -
            strict_anchor_length;

        bool destabilizing_found =
            false;

        PrimerBindingSite mismatch_site;

        BindingSiteThermodynamicProfile
            mismatch_profile;

        std::size_t selected_position =
            std::numeric_limits<
                std::size_t
            >::max();

        char selected_base =
            'N';


        for (
            std::size_t position = 0;
            position < mutable_prefix &&
            !destabilizing_found;
            ++position
        ) {

            for (
                std::size_t variant = 0;
                variant < 3;
                ++variant
            ) {

                std::string binding =
                    primer;

                binding[position] =
                    alternative_base(
                        primer[position],
                        variant
                    );

                const auto candidate_site =
                    make_site(
                        primer,
                        binding,
                        false,
                        PrimerIdentity::Primer1
                    );

                const auto candidate_profile =
                    profile_binding_site_thermodynamics(
                        candidate_site,
                        backend,
                        conditions
                    );

                if (
                    candidate_profile
                        .delta_tm_celsius >
                    0.01
                ) {
                    mismatch_site =
                        candidate_site;

                    mismatch_profile =
                        candidate_profile;

                    selected_position =
                        position;

                    selected_base =
                        binding[position];

                    destabilizing_found =
                        true;

                    break;
                }
            }
        }


        expect(
            destabilizing_found,
            "Destabilizing non-anchor mismatch found"
        );

        expect(
            selected_position <
                mutable_prefix,
            "Selected mismatch remains outside strict 3-prime anchor"
        );

        expect(
            mismatch_profile.mismatch_count == 1,
            "Mismatch profile contains one mismatch"
        );

        expect(
            std::popcount(
                mismatch_site.mismatch_mask
            ) == 1,
            "Mismatch mask contains one mismatch"
        );

        expect(
            mismatch_profile
                .perfect_match_tm_celsius >
            mismatch_profile
                .observed_binding_tm_celsius,
            "Selected mismatch lowers duplex Tm"
        );

        expect(
            mismatch_profile
                .delta_tm_celsius >
            0.0,
            "Selected mismatch has positive delta-Tm"
        );

        expect(
            almost_equal(
                mismatch_profile
                    .delta_tm_celsius,
                mismatch_profile
                    .perfect_match_tm_celsius -
                mismatch_profile
                    .observed_binding_tm_celsius
            ),
            "Delta-Tm convention is perfect minus observed"
        );


        /*
         * Direct backend equivalence:
         *
         * PrimerBindingSite stores the binding sequence
         * in primer orientation.
         *
         * ntthal must receive the antiparallel target.
         */
        const std::string observed_target =
            reverse_complement_binding_sequence(
                mismatch_site
                    .binding_sequence
            );

        const double direct_observed =
            backend.duplex_tm(
                primer,
                observed_target,
                ThermodynamicAlignment::Any,
                conditions
            );

        expect(
            almost_equal(
                direct_observed,
                mismatch_profile
                    .observed_binding_tm_celsius
            ),
            "Binding-site adapter equals direct ntthal orientation"
        );


        /*
         * Reverse-strand metadata must not change
         * thermodynamic sequence interpretation because
         * binding_sequence is already normalized to the
         * primer's biological orientation.
         */
        auto reverse_site =
            mismatch_site;

        reverse_site.reverse_strand =
            true;

        reverse_site.primer =
            PrimerIdentity::Primer2;

        reverse_site.genomic_sequence =
            reverse_complement_binding_sequence(
                reverse_site
                    .binding_sequence
            );


        const auto reverse_profile =
            profile_binding_site_thermodynamics(
                reverse_site,
                backend,
                conditions
            );


        expect(
            reverse_profile.reverse_strand,
            "Reverse-strand metadata preserved"
        );

        expect(
            reverse_profile.primer ==
                PrimerIdentity::Primer2,
            "Primer identity preserved"
        );

        expect(
            almost_equal(
                reverse_profile
                    .observed_binding_tm_celsius,
                mismatch_profile
                    .observed_binding_tm_celsius
            ),
            "Normalized binding sequence gives strand-invariant Tm"
        );


        PrimerPairBindingSites pair_sites;

        pair_sites.left =
            perfect_site;

        pair_sites.right =
            reverse_site;


        const auto pair_profile =
            profile_pair_binding_thermodynamics(
                pair_sites,
                backend,
                conditions
            );


        expect(
            almost_equal(
                pair_profile
                    .left
                    .delta_tm_celsius,
                0.0
            ),
            "Pair left perfect delta-Tm zero"
        );

        expect(
            pair_profile
                .right
                .delta_tm_celsius >
            0.0,
            "Pair right mismatch delta-Tm positive"
        );

        expect(
            almost_equal(
                pair_profile
                    .mean_delta_tm_celsius,
                (
                    pair_profile
                        .left
                        .delta_tm_celsius +
                    pair_profile
                        .right
                        .delta_tm_celsius
                ) /
                2.0
            ),
            "Pair mean delta-Tm aggregation exact"
        );

        expect(
            almost_equal(
                pair_profile
                    .max_delta_tm_celsius,
                pair_profile
                    .right
                    .delta_tm_celsius
            ),
            "Pair max delta-Tm aggregation exact"
        );


        bool length_mismatch_rejected =
            false;

        try {

            auto invalid =
                perfect_site;

            invalid.binding_sequence.pop_back();

            static_cast<void>(
                profile_binding_site_thermodynamics(
                    invalid,
                    backend,
                    conditions
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            length_mismatch_rejected =
                true;
        }

        expect(
            length_mismatch_rejected,
            "Primer/binding length disagreement rejected"
        );


        bool mask_disagreement_rejected =
            false;

        try {

            auto invalid =
                mismatch_site;

            invalid.mismatch_mask =
                0;

            static_cast<void>(
                profile_binding_site_thermodynamics(
                    invalid,
                    backend,
                    conditions
                )
            );

        } catch (
            const std::logic_error&
        ) {
            mask_disagreement_rejected =
                true;
        }

        expect(
            mask_disagreement_rejected,
            "Mismatch-mask disagreement rejected"
        );


        bool vector_disagreement_rejected =
            false;

        try {

            auto invalid =
                mismatch_site;

            invalid.mismatches.clear();

            static_cast<void>(
                profile_binding_site_thermodynamics(
                    invalid,
                    backend,
                    conditions
                )
            );

        } catch (
            const std::logic_error&
        ) {
            vector_disagreement_rejected =
                true;
        }

        expect(
            vector_disagreement_rejected,
            "Mismatch-vector disagreement rejected"
        );


        std::cout
            << "perfect_tm_c\t"
            << perfect_profile
                   .perfect_match_tm_celsius
            << '\n';

        std::cout
            << "offtarget_tm_c\t"
            << mismatch_profile
                   .observed_binding_tm_celsius
            << '\n';

        std::cout
            << "delta_tm_c\t"
            << mismatch_profile
                   .delta_tm_celsius
            << '\n';

        std::cout
            << "selected_mismatch_position_5prime\t"
            << selected_position
            << '\n';

        std::cout
            << "selected_mismatch_distance_3prime\t"
            << (
                primer.size() -
                1 -
                selected_position
            )
            << '\n';

        std::cout
            << "selected_mismatch_base\t"
            << selected_base
            << '\n';


        std::cout
            << "BINDING_SITE_PERFECT_TM_VALID\tYES\n";

        std::cout
            << "BINDING_SITE_OBSERVED_TM_VALID\tYES\n";

        std::cout
            << "BINDING_SITE_DELTA_TM_VALID\tYES\n";

        std::cout
            << "BINDING_SITE_TARGET_ORIENTATION_VALID\tYES\n";

        std::cout
            << "BINDING_SITE_STRAND_NORMALIZATION_VALID\tYES\n";

        std::cout
            << "BINDING_SITE_PAIR_AGGREGATION_VALID\tYES\n";

        std::cout
            << "BINDING_SITE_THERMODYNAMICS_PCR_PROBABILITY\tNO\n";

        std::cout
            << "BINDING_SITE_THERMODYNAMICS_RISK_SCORE\tNO\n";

        std::cout
            << "BINDING_SITE_THERMODYNAMICS_V1_COMPLETE\tYES\n";

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
