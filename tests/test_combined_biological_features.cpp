/*
 * COMBINED_BIOLOGICAL_FEATURE_VECTOR_V1
 * Test #72
 */

#include "primerpair/combined_biological_features.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

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
    const double tolerance = 1.0e-5
)
{
    return
        std::abs(
            lhs -
            rhs
        ) <=
        tolerance;
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

        if (
            seen ==
            variant_index
        ) {
            return base;
        }

        ++seen;
    }

    throw std::runtime_error(
        "Could not select alternative DNA base."
    );
}


primerpair::PrimerBindingSite
make_site(
    const std::string& primer,
    const std::string& binding,
    const bool reverse_strand,
    const primerpair::PrimerIdentity identity,
    const std::uint64_t genomic_start
)
{
    using namespace primerpair;

    if (
        primer.size() !=
        binding.size()
    ) {
        throw std::runtime_error(
            "Test primer/binding lengths disagree."
        );
    }


    PrimerBindingSite site;

    site.primer =
        identity;

    site.reverse_strand =
        reverse_strand;

    site.genomic_start =
        genomic_start;

    site.genomic_end_exclusive =
        genomic_start +
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


    site.mismatch_mask = 0;


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


        const std::string primer1 =
            "GCGTACGATCGTACGCATGC";


        const std::string primer2 =
            "CTGACGTAGCTAGTCGATGC";


        constexpr std::uint64_t
            left_start = 1000;

        constexpr std::uint64_t
            right_start = 1200;


        const PrimerBindingSite left_site =
            make_site(
                primer1,
                primer1,
                false,
                PrimerIdentity::Primer1,
                left_start
            );


        /*
         * Search for a deterministic mismatch in
         * the 5-prime prefix, outside the strict
         * terminal-12 biological anchor.
         */
        constexpr std::size_t
            strict_anchor_length = 12;


        const std::size_t mutable_prefix =
            primer2.size() -
            strict_anchor_length;


        bool destabilizing_found =
            false;


        PrimerBindingSite right_site;


        BindingSiteThermodynamicProfile
            selected_binding_profile;


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
                    primer2;


                binding[position] =
                    alternative_base(
                        primer2[position],
                        variant
                    );


                const PrimerBindingSite candidate =
                    make_site(
                        primer2,
                        binding,
                        true,
                        PrimerIdentity::Primer2,
                        right_start
                    );


                const auto candidate_profile =
                    profile_binding_site_thermodynamics(
                        candidate,
                        backend,
                        conditions
                    );


                if (
                    candidate_profile
                        .delta_tm_celsius >
                    0.01
                ) {

                    right_site =
                        candidate;

                    selected_binding_profile =
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
            "Selected mismatch outside strict 3-prime anchor"
        );


        PrimerPairBindingSites sites;

        sites.left =
            left_site;

        sites.right =
            right_site;


        PrimerPairHit hit;

        hit.left_primer =
            PrimerIdentity::Primer1;

        hit.right_primer =
            PrimerIdentity::Primer2;

        hit.left_position =
            left_site.genomic_start;

        hit.right_position =
            right_site.genomic_start;

        hit.left_mismatches =
            left_site.mismatch_count();

        hit.right_mismatches =
            right_site.mismatch_count();

        hit.amplicon_start =
            left_site.genomic_start;

        hit.amplicon_end_exclusive =
            right_site.genomic_end_exclusive;

        hit.amplicon_length =
            hit.amplicon_end_exclusive -
            hit.amplicon_start;

        hit.left_mismatch_mask =
            left_site.mismatch_mask;

        hit.right_mismatch_mask =
            right_site.mismatch_mask;


        const auto combined =
            build_combined_biological_features(
                hit,
                sites,
                backend,
                conditions
            );


        expect(
            combined.left.primer ==
                PrimerIdentity::Primer1,
            "Left primer identity preserved"
        );


        expect(
            combined.right.primer ==
                PrimerIdentity::Primer2,
            "Right primer identity preserved"
        );


        expect(
            !combined.left.reverse_strand,
            "Left strand remains forward"
        );


        expect(
            combined.right.reverse_strand,
            "Right strand remains reverse"
        );


        expect(
            combined.left.mismatch_count == 0,
            "Perfect left primer has zero mismatches"
        );


        expect(
            combined.right.mismatch_count == 1,
            "Right off-target contains one mismatch"
        );


        expect(
            combined.total_mismatches == 1,
            "Pair total mismatch count exact"
        );


        expect(
            combined.amplicon_start ==
                left_start,
            "Amplicon start preserved"
        );


        expect(
            combined.amplicon_end_exclusive ==
                right_site.genomic_end_exclusive,
            "Amplicon end preserved"
        );


        expect(
            combined.amplicon_length ==
                hit.amplicon_length,
            "Amplicon length preserved"
        );


        /*
         * Independent biological mismatch profile.
         */
        const auto right_biological =
            build_biological_mismatch_profile(
                right_site.mismatch_mask,
                primer2.size()
            );


        expect(
            almost_equal(
                combined.right
                    .mismatch_fraction,
                right_biological
                    .mismatch_fraction
            ),
            "Mismatch fraction equals biological profile"
        );


        expect(
            combined.right.last_1_count ==
                right_biological.last_1_count &&
            combined.right.last_2_count ==
                right_biological.last_2_count &&
            combined.right.last_3_count ==
                right_biological.last_3_count &&
            combined.right.last_5_count ==
                right_biological.last_5_count &&
            combined.right.last_8_count ==
                right_biological.last_8_count &&
            combined.right.last_12_count ==
                right_biological.last_12_count,
            "3-prime mismatch windows equal biological profile"
        );


        expect(
            combined.right
                .nearest_mismatch_to_3prime ==
            right_biological
                .nearest_mismatch_to_3prime,
            "Nearest 3-prime mismatch distance preserved"
        );


        expect(
            combined.right
                .exact_3prime_run_length ==
            right_biological
                .exact_3prime_run_length,
            "Exact 3-prime run length preserved"
        );


        expect(
            almost_equal(
                combined.right
                    .normalized_3prime_positional_burden,
                right_biological
                    .normalized_3prime_positional_burden
            ),
            "3-prime positional burden preserved"
        );


        expect(
            combined.right.last_12_count == 0,
            "Strict terminal-12 anchor remains mismatch-free"
        );


        /*
         * Independent binding-site thermodynamics.
         */
        expect(
            almost_equal(
                combined.right
                    .perfect_match_tm_celsius,
                selected_binding_profile
                    .perfect_match_tm_celsius
            ),
            "Perfect binding Tm equals Test #71 layer"
        );


        expect(
            almost_equal(
                combined.right
                    .observed_binding_tm_celsius,
                selected_binding_profile
                    .observed_binding_tm_celsius
            ),
            "Observed binding Tm equals Test #71 layer"
        );


        expect(
            almost_equal(
                combined.right
                    .delta_tm_celsius,
                selected_binding_profile
                    .delta_tm_celsius
            ),
            "Delta-Tm equals Test #71 layer"
        );


        expect(
            combined.right
                .delta_tm_celsius >
            0.0,
            "Selected off-target has positive delta-Tm"
        );


        /*
         * Independent intrinsic primer/pair
         * thermodynamic calculation.
         */
        const auto direct_pair_thermo =
            backend.profile_pair(
                left_site.primer_sequence,
                right_site.primer_sequence,
                conditions
            );


        expect(
            almost_equal(
                combined.left
                    .oligo_tm_celsius,
                direct_pair_thermo
                    .left
                    .oligo_tm_celsius
            ),
            "Left oligo Tm exact"
        );


        expect(
            almost_equal(
                combined.right
                    .oligo_tm_celsius,
                direct_pair_thermo
                    .right
                    .oligo_tm_celsius
            ),
            "Right oligo Tm exact"
        );


        expect(
            almost_equal(
                combined.right
                    .hairpin_tm_celsius,
                direct_pair_thermo
                    .right
                    .hairpin_tm_celsius
            ),
            "Right hairpin Tm exact"
        );


        expect(
            almost_equal(
                combined.right
                    .homodimer_any_tm_celsius,
                direct_pair_thermo
                    .right
                    .homodimer_any_tm_celsius
            ),
            "Right homodimer ANY Tm exact"
        );


        expect(
            almost_equal(
                combined.right
                    .homodimer_end1_tm_celsius,
                direct_pair_thermo
                    .right
                    .homodimer_end1_tm_celsius
            ),
            "Right homodimer END1 Tm exact"
        );


        expect(
            almost_equal(
                combined.right
                    .homodimer_end2_tm_celsius,
                direct_pair_thermo
                    .right
                    .homodimer_end2_tm_celsius
            ),
            "Right homodimer END2 Tm exact"
        );


        expect(
            almost_equal(
                combined
                    .heterodimer_any_tm_celsius,
                direct_pair_thermo
                    .heterodimer_any_tm_celsius
            ),
            "Pair heterodimer ANY Tm exact"
        );


        expect(
            almost_equal(
                combined
                    .heterodimer_end1_tm_celsius,
                direct_pair_thermo
                    .heterodimer_end1_tm_celsius
            ),
            "Pair heterodimer END1 Tm exact"
        );


        expect(
            almost_equal(
                combined
                    .heterodimer_end2_tm_celsius,
                direct_pair_thermo
                    .heterodimer_end2_tm_celsius
            ),
            "Pair heterodimer END2 Tm exact"
        );


        expect(
            almost_equal(
                combined.mean_mismatch_fraction,
                (
                    combined.left.mismatch_fraction +
                    combined.right.mismatch_fraction
                ) /
                2.0
            ),
            "Mean mismatch fraction aggregation exact"
        );


        expect(
            almost_equal(
                combined.max_mismatch_fraction,
                std::max(
                    combined.left.mismatch_fraction,
                    combined.right.mismatch_fraction
                )
            ),
            "Max mismatch fraction aggregation exact"
        );


        expect(
            almost_equal(
                combined
                    .mean_normalized_3prime_positional_burden,
                (
                    combined.left
                        .normalized_3prime_positional_burden +
                    combined.right
                        .normalized_3prime_positional_burden
                ) /
                2.0
            ),
            "Mean 3-prime burden aggregation exact"
        );


        expect(
            almost_equal(
                combined
                    .max_normalized_3prime_positional_burden,
                std::max(
                    combined.left
                        .normalized_3prime_positional_burden,
                    combined.right
                        .normalized_3prime_positional_burden
                )
            ),
            "Max 3-prime burden aggregation exact"
        );


        expect(
            almost_equal(
                combined.mean_delta_tm_celsius,
                (
                    combined.left.delta_tm_celsius +
                    combined.right.delta_tm_celsius
                ) /
                2.0
            ),
            "Mean delta-Tm aggregation exact"
        );


        expect(
            almost_equal(
                combined.max_delta_tm_celsius,
                std::max(
                    combined.left.delta_tm_celsius,
                    combined.right.delta_tm_celsius
                )
            ),
            "Max delta-Tm aggregation exact"
        );


        /*
         * Existing search-hit/binding-site invariants
         * must be enforced rather than silently
         * reconciled.
         */
        bool mismatch_mask_disagreement_rejected =
            false;

        try {

            auto invalid_hit =
                hit;

            invalid_hit.right_mismatch_mask =
                0;

            static_cast<void>(
                build_combined_biological_features(
                    invalid_hit,
                    sites,
                    backend,
                    conditions
                )
            );

        } catch (
            const std::logic_error&
        ) {
            mismatch_mask_disagreement_rejected =
                true;
        }


        expect(
            mismatch_mask_disagreement_rejected,
            "Hit/site mismatch-mask disagreement rejected"
        );


        bool geometry_disagreement_rejected =
            false;

        try {

            auto invalid_hit =
                hit;

            ++invalid_hit.amplicon_length;

            static_cast<void>(
                build_combined_biological_features(
                    invalid_hit,
                    sites,
                    backend,
                    conditions
                )
            );

        } catch (
            const std::logic_error&
        ) {
            geometry_disagreement_rejected =
                true;
        }


        expect(
            geometry_disagreement_rejected,
            "Amplicon geometry disagreement rejected"
        );


        std::cout
            << "selected_mismatch_position_5prime\t"
            << selected_position
            << '\n';


        std::cout
            << "selected_mismatch_distance_3prime\t"
            << (
                primer2.size() -
                1 -
                selected_position
            )
            << '\n';


        std::cout
            << "selected_mismatch_base\t"
            << selected_base
            << '\n';


        std::cout
            << "right_mismatch_fraction\t"
            << combined.right
                   .mismatch_fraction
            << '\n';


        std::cout
            << "right_3prime_burden\t"
            << combined.right
                   .normalized_3prime_positional_burden
            << '\n';


        std::cout
            << "right_exact_3prime_run\t"
            << combined.right
                   .exact_3prime_run_length
            << '\n';


        std::cout
            << "right_delta_tm_c\t"
            << combined.right
                   .delta_tm_celsius
            << '\n';


        std::cout
            << "amplicon_length\t"
            << combined.amplicon_length
            << '\n';


        std::cout
            << "pair_heterodimer_any_tm_c\t"
            << combined
                   .heterodimer_any_tm_celsius
            << '\n';


        std::cout
            << "COMBINED_FEATURE_VECTOR_COMPONENT_EQUIVALENCE\tYES\n";

        std::cout
            << "COMBINED_FEATURE_VECTOR_3PRIME_VALID\tYES\n";

        std::cout
            << "COMBINED_FEATURE_VECTOR_THERMODYNAMICS_VALID\tYES\n";

        std::cout
            << "COMBINED_FEATURE_VECTOR_PAIR_GEOMETRY_VALID\tYES\n";

        std::cout
            << "COMBINED_FEATURE_VECTOR_NO_SEARCH_CHANGE\tYES\n";

        std::cout
            << "COMBINED_FEATURE_VECTOR_FINAL_RISK_SCORE\tNO\n";

        std::cout
            << "COMBINED_FEATURE_VECTOR_PCR_PROBABILITY\tNO\n";

        std::cout
            << "COMBINED_BIOLOGICAL_FEATURE_VECTOR_V1_COMPLETE\tYES\n";

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
