/*
 * COMBINED_BIOLOGICAL_FEATURE_VECTOR_V1
 */

#include "primerpair/combined_biological_features.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace primerpair {

namespace {


void validate_finite(
    const double value,
    const char* name
)
{
    if (!std::isfinite(value)) {
        throw std::logic_error(
            std::string(
                "Combined biological feature "
                "contains non-finite value: "
            ) +
            name
        );
    }
}


void validate_site_span(
    const PrimerBindingSite& site
)
{
    if (
        site.primer_sequence.empty()
    ) {
        throw std::invalid_argument(
            "PrimerBindingSite primer sequence "
            "must not be empty."
        );
    }

    if (
        site.genomic_end_exclusive <
        site.genomic_start
    ) {
        throw std::logic_error(
            "PrimerBindingSite genomic span "
            "is inverted."
        );
    }

    const std::uint64_t span =
        site.genomic_end_exclusive -
        site.genomic_start;

    if (
        span !=
        site.primer_sequence.size()
    ) {
        throw std::logic_error(
            "PrimerBindingSite genomic span "
            "does not equal primer length."
        );
    }
}


void validate_pair_consistency(
    const PrimerPairHit& hit,
    const PrimerPairBindingSites& sites
)
{
    validate_site_span(
        sites.left
    );

    validate_site_span(
        sites.right
    );


    /*
     * PrimerPairHit invariant:
     *
     * genomic LEFT = forward strand
     * genomic RIGHT = reverse strand
     */
    if (sites.left.reverse_strand) {
        throw std::logic_error(
            "Left PrimerBindingSite must be "
            "forward strand."
        );
    }

    if (!sites.right.reverse_strand) {
        throw std::logic_error(
            "Right PrimerBindingSite must be "
            "reverse strand."
        );
    }


    if (
        sites.left.primer !=
        hit.left_primer
    ) {
        throw std::logic_error(
            "Left primer identity disagrees "
            "with PrimerPairHit."
        );
    }

    if (
        sites.right.primer !=
        hit.right_primer
    ) {
        throw std::logic_error(
            "Right primer identity disagrees "
            "with PrimerPairHit."
        );
    }


    if (
        sites.left.genomic_start !=
        hit.left_position
    ) {
        throw std::logic_error(
            "Left binding-site position "
            "disagrees with PrimerPairHit."
        );
    }

    if (
        sites.right.genomic_start !=
        hit.right_position
    ) {
        throw std::logic_error(
            "Right binding-site position "
            "disagrees with PrimerPairHit."
        );
    }


    if (
        sites.left.mismatch_mask !=
        hit.left_mismatch_mask
    ) {
        throw std::logic_error(
            "Left mismatch mask disagrees "
            "with PrimerPairHit."
        );
    }

    if (
        sites.right.mismatch_mask !=
        hit.right_mismatch_mask
    ) {
        throw std::logic_error(
            "Right mismatch mask disagrees "
            "with PrimerPairHit."
        );
    }


    if (
        sites.left.mismatch_count() !=
        hit.left_mismatches
    ) {
        throw std::logic_error(
            "Left mismatch count disagrees "
            "with PrimerPairHit."
        );
    }

    if (
        sites.right.mismatch_count() !=
        hit.right_mismatches
    ) {
        throw std::logic_error(
            "Right mismatch count disagrees "
            "with PrimerPairHit."
        );
    }


    if (
        hit.amplicon_end_exclusive <
        hit.amplicon_start
    ) {
        throw std::logic_error(
            "PrimerPairHit amplicon span "
            "is inverted."
        );
    }


    const std::uint64_t expected_length =
        hit.amplicon_end_exclusive -
        hit.amplicon_start;


    if (
        expected_length !=
        hit.amplicon_length
    ) {
        throw std::logic_error(
            "PrimerPairHit amplicon length "
            "disagrees with coordinates."
        );
    }


    if (
        hit.amplicon_start !=
        sites.left.genomic_start
    ) {
        throw std::logic_error(
            "Amplicon start does not equal "
            "left primer start."
        );
    }


    if (
        hit.amplicon_end_exclusive !=
        sites.right.genomic_end_exclusive
    ) {
        throw std::logic_error(
            "Amplicon end does not equal "
            "right primer end."
        );
    }
}


PrimerCombinedBiologicalFeatures
combine_primer_features(
    const PrimerBindingSite& site,

    const BiologicalMismatchProfile&
        biological,

    const BindingSiteThermodynamicProfile&
        binding,

    const PrimerThermodynamicProfile&
        thermodynamic
)
{
    if (
        biological.primer_length !=
        site.primer_sequence.size()
    ) {
        throw std::logic_error(
            "Biological mismatch profile "
            "primer length disagrees with "
            "binding site."
        );
    }


    if (
        biological.mismatch_count !=
        site.mismatch_count()
    ) {
        throw std::logic_error(
            "Biological mismatch count "
            "disagrees with binding site."
        );
    }


    if (
        binding.primer_length !=
        site.primer_sequence.size()
    ) {
        throw std::logic_error(
            "Binding thermodynamic profile "
            "primer length disagrees with "
            "binding site."
        );
    }


    if (
        binding.mismatch_count !=
        site.mismatch_count()
    ) {
        throw std::logic_error(
            "Binding thermodynamic mismatch "
            "count disagrees with binding site."
        );
    }


    PrimerCombinedBiologicalFeatures output;

    output.primer =
        site.primer;

    output.reverse_strand =
        site.reverse_strand;

    output.genomic_start =
        site.genomic_start;

    output.genomic_end_exclusive =
        site.genomic_end_exclusive;

    output.primer_length =
        biological.primer_length;


    output.mismatch_count =
        biological.mismatch_count;

    output.mismatch_fraction =
        biological.mismatch_fraction;

    output.last_1_count =
        biological.last_1_count;

    output.last_2_count =
        biological.last_2_count;

    output.last_3_count =
        biological.last_3_count;

    output.last_5_count =
        biological.last_5_count;

    output.last_8_count =
        biological.last_8_count;

    output.last_12_count =
        biological.last_12_count;

    output.nearest_mismatch_to_3prime =
        biological
            .nearest_mismatch_to_3prime;

    output.exact_3prime_run_length =
        biological
            .exact_3prime_run_length;

    output
        .normalized_3prime_positional_burden =
            biological
                .normalized_3prime_positional_burden;


    output.perfect_match_tm_celsius =
        binding
            .perfect_match_tm_celsius;

    output.observed_binding_tm_celsius =
        binding
            .observed_binding_tm_celsius;

    output.delta_tm_celsius =
        binding
            .delta_tm_celsius;


    output.oligo_tm_celsius =
        thermodynamic
            .oligo_tm_celsius;

    output.hairpin_tm_celsius =
        thermodynamic
            .hairpin_tm_celsius;

    output.homodimer_any_tm_celsius =
        thermodynamic
            .homodimer_any_tm_celsius;

    output.homodimer_end1_tm_celsius =
        thermodynamic
            .homodimer_end1_tm_celsius;

    output.homodimer_end2_tm_celsius =
        thermodynamic
            .homodimer_end2_tm_celsius;


    validate_finite(
        output.mismatch_fraction,
        "mismatch_fraction"
    );

    validate_finite(
        output
            .normalized_3prime_positional_burden,
        "normalized_3prime_positional_burden"
    );

    validate_finite(
        output.perfect_match_tm_celsius,
        "perfect_match_tm_celsius"
    );

    validate_finite(
        output.observed_binding_tm_celsius,
        "observed_binding_tm_celsius"
    );

    validate_finite(
        output.delta_tm_celsius,
        "delta_tm_celsius"
    );

    validate_finite(
        output.oligo_tm_celsius,
        "oligo_tm_celsius"
    );

    validate_finite(
        output.hairpin_tm_celsius,
        "hairpin_tm_celsius"
    );

    validate_finite(
        output.homodimer_any_tm_celsius,
        "homodimer_any_tm_celsius"
    );

    validate_finite(
        output.homodimer_end1_tm_celsius,
        "homodimer_end1_tm_celsius"
    );

    validate_finite(
        output.homodimer_end2_tm_celsius,
        "homodimer_end2_tm_celsius"
    );


    return output;
}


}  // namespace


PrimerPairCombinedBiologicalFeatures
build_combined_biological_features(
    const PrimerPairHit& hit,
    const PrimerPairBindingSites& sites,
    const Primer3ThermodynamicBackend& backend,
    const ThermodynamicConditions& conditions
)
{
    validate_pair_consistency(
        hit,
        sites
    );


    const BiologicalMismatchProfile
        left_biological =
            build_biological_mismatch_profile(
                sites.left.mismatch_mask,
                sites.left
                    .primer_sequence
                    .size()
            );


    const BiologicalMismatchProfile
        right_biological =
            build_biological_mismatch_profile(
                sites.right.mismatch_mask,
                sites.right
                    .primer_sequence
                    .size()
            );


    const BindingSiteThermodynamicProfile
        left_binding =
            profile_binding_site_thermodynamics(
                sites.left,
                backend,
                conditions
            );


    const BindingSiteThermodynamicProfile
        right_binding =
            profile_binding_site_thermodynamics(
                sites.right,
                backend,
                conditions
            );


    /*
     * Compute primer-intrinsic and primer-primer
     * thermodynamics once for the complete pair.
     */
    const PrimerPairThermodynamicProfile
        pair_thermodynamic =
            backend.profile_pair(
                sites.left.primer_sequence,
                sites.right.primer_sequence,
                conditions
            );


    PrimerPairCombinedBiologicalFeatures output;


    output.left =
        combine_primer_features(
            sites.left,
            left_biological,
            left_binding,
            pair_thermodynamic.left
        );


    output.right =
        combine_primer_features(
            sites.right,
            right_biological,
            right_binding,
            pair_thermodynamic.right
        );


    output.amplicon_start =
        hit.amplicon_start;

    output.amplicon_end_exclusive =
        hit.amplicon_end_exclusive;

    output.amplicon_length =
        hit.amplicon_length;


    output.total_mismatches =
        output.left.mismatch_count +
        output.right.mismatch_count;


    if (
        output.total_mismatches !=
        hit.total_mismatches()
    ) {
        throw std::logic_error(
            "Combined feature total mismatch "
            "count disagrees with PrimerPairHit."
        );
    }


    output.mean_mismatch_fraction =
        (
            output.left.mismatch_fraction +
            output.right.mismatch_fraction
        ) /
        2.0;


    output.max_mismatch_fraction =
        std::max(
            output.left.mismatch_fraction,
            output.right.mismatch_fraction
        );


    output
        .mean_normalized_3prime_positional_burden =
            (
                output.left
                    .normalized_3prime_positional_burden +
                output.right
                    .normalized_3prime_positional_burden
            ) /
            2.0;


    output
        .max_normalized_3prime_positional_burden =
            std::max(
                output.left
                    .normalized_3prime_positional_burden,
                output.right
                    .normalized_3prime_positional_burden
            );


    output.mean_delta_tm_celsius =
        (
            output.left.delta_tm_celsius +
            output.right.delta_tm_celsius
        ) /
        2.0;


    output.max_delta_tm_celsius =
        std::max(
            output.left.delta_tm_celsius,
            output.right.delta_tm_celsius
        );


    output.heterodimer_any_tm_celsius =
        pair_thermodynamic
            .heterodimer_any_tm_celsius;

    output.heterodimer_end1_tm_celsius =
        pair_thermodynamic
            .heterodimer_end1_tm_celsius;

    output.heterodimer_end2_tm_celsius =
        pair_thermodynamic
            .heterodimer_end2_tm_celsius;


    validate_finite(
        output.mean_mismatch_fraction,
        "mean_mismatch_fraction"
    );

    validate_finite(
        output.max_mismatch_fraction,
        "max_mismatch_fraction"
    );

    validate_finite(
        output
            .mean_normalized_3prime_positional_burden,
        "mean_normalized_3prime_positional_burden"
    );

    validate_finite(
        output
            .max_normalized_3prime_positional_burden,
        "max_normalized_3prime_positional_burden"
    );

    validate_finite(
        output.mean_delta_tm_celsius,
        "mean_delta_tm_celsius"
    );

    validate_finite(
        output.max_delta_tm_celsius,
        "max_delta_tm_celsius"
    );

    validate_finite(
        output.heterodimer_any_tm_celsius,
        "heterodimer_any_tm_celsius"
    );

    validate_finite(
        output.heterodimer_end1_tm_celsius,
        "heterodimer_end1_tm_celsius"
    );

    validate_finite(
        output.heterodimer_end2_tm_celsius,
        "heterodimer_end2_tm_celsius"
    );


    return output;
}


}  // namespace primerpair
