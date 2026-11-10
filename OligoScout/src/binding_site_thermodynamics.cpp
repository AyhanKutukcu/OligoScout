/*
 * BINDING_SITE_THERMODYNAMICS_V1
 */

#include "primerpair/binding_site_thermodynamics.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace primerpair {

namespace {


bool canonical_base(
    const char base
) noexcept
{
    return
        base == 'A' ||
        base == 'C' ||
        base == 'G' ||
        base == 'T';
}


void validate_canonical_sequence(
    const std::string_view sequence,
    const char* name
)
{
    if (sequence.empty()) {
        throw std::invalid_argument(
            std::string(name) +
            " cannot be empty."
        );
    }

    for (const char base : sequence) {
        if (!canonical_base(base)) {
            throw std::invalid_argument(
                std::string(name) +
                " contains a non-canonical DNA base."
            );
        }
    }
}


std::uint64_t expected_mismatch_mask(
    const std::string_view primer,
    const std::string_view binding
)
{
    if (
        primer.size() !=
        binding.size()
    ) {
        throw std::invalid_argument(
            "Primer and binding sequence lengths "
            "must match."
        );
    }

    if (primer.size() > 64) {
        throw std::invalid_argument(
            "Binding-site mismatch mask supports "
            "primer lengths up to 64 bases."
        );
    }

    std::uint64_t mask = 0;

    for (
        std::size_t i = 0;
        i < primer.size();
        ++i
    ) {
        if (
            primer[i] !=
            binding[i]
        ) {
            mask |=
                (
                    std::uint64_t{1}
                    <<
                    i
                );
        }
    }

    return mask;
}


void validate_site(
    const PrimerBindingSite& site
)
{
    validate_canonical_sequence(
        site.primer_sequence,
        "primer_sequence"
    );

    validate_canonical_sequence(
        site.binding_sequence,
        "binding_sequence"
    );

    if (
        site.primer_sequence.size() !=
        site.binding_sequence.size()
    ) {
        throw std::invalid_argument(
            "PrimerBindingSite primer/binding "
            "sequence lengths disagree."
        );
    }

    if (
        site.genomic_end_exclusive <
        site.genomic_start
    ) {
        throw std::invalid_argument(
            "PrimerBindingSite genomic interval "
            "is invalid."
        );
    }

    const std::uint64_t genomic_length =
        site.genomic_end_exclusive -
        site.genomic_start;

    if (
        genomic_length !=
        site.primer_sequence.size()
    ) {
        throw std::invalid_argument(
            "PrimerBindingSite genomic interval "
            "length disagrees with primer length."
        );
    }

    const std::uint64_t expected =
        expected_mismatch_mask(
            site.primer_sequence,
            site.binding_sequence
        );

    if (
        expected !=
        site.mismatch_mask
    ) {
        throw std::logic_error(
            "PrimerBindingSite mismatch mask "
            "disagrees with primer/binding "
            "sequences."
        );
    }

    if (
        static_cast<std::size_t>(
            std::popcount(
                expected
            )
        ) !=
        site.mismatch_count()
    ) {
        throw std::logic_error(
            "PrimerBindingSite mismatch vector "
            "disagrees with mismatch mask."
        );
    }
}


double require_finite(
    const double value,
    const char* name
)
{
    if (!std::isfinite(value)) {
        throw std::runtime_error(
            std::string(
                "Non-finite thermodynamic value: "
            ) +
            name
        );
    }

    return value;
}


}  // namespace


BindingSiteThermodynamicProfile
profile_binding_site_thermodynamics(
    const PrimerBindingSite& site,
    const Primer3ThermodynamicBackend& backend,
    const ThermodynamicConditions& conditions
)
{
    validate_site(
        site
    );

    /*
     * The binding sequence is normalized into the
     * same biological 5' -> 3' orientation as the
     * primer by PrimerBindingSite.
     *
     * ntthal needs the antiparallel partner.
     */
    const std::string perfect_target =
        reverse_complement_binding_sequence(
            site.primer_sequence
        );

    const std::string observed_target =
        reverse_complement_binding_sequence(
            site.binding_sequence
        );


    const double perfect_tm =
        require_finite(
            backend.duplex_tm(
                site.primer_sequence,
                perfect_target,
                ThermodynamicAlignment::Any,
                conditions
            ),
            "perfect_match_tm"
        );


    const double observed_tm =
        require_finite(
            backend.duplex_tm(
                site.primer_sequence,
                observed_target,
                ThermodynamicAlignment::Any,
                conditions
            ),
            "observed_binding_tm"
        );


    BindingSiteThermodynamicProfile output;

    output.primer =
        site.primer;

    output.reverse_strand =
        site.reverse_strand;

    output.primer_length =
        site.primer_sequence.size();

    output.mismatch_count =
        site.mismatch_count();

    output.perfect_match_tm_celsius =
        perfect_tm;

    output.observed_binding_tm_celsius =
        observed_tm;

    output.delta_tm_celsius =
        perfect_tm -
        observed_tm;

    return output;
}


PrimerPairBindingThermodynamicProfile
profile_pair_binding_thermodynamics(
    const PrimerPairBindingSites& sites,
    const Primer3ThermodynamicBackend& backend,
    const ThermodynamicConditions& conditions
)
{
    PrimerPairBindingThermodynamicProfile output;

    output.left =
        profile_binding_site_thermodynamics(
            sites.left,
            backend,
            conditions
        );

    output.right =
        profile_binding_site_thermodynamics(
            sites.right,
            backend,
            conditions
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

    return output;
}


}  // namespace primerpair
