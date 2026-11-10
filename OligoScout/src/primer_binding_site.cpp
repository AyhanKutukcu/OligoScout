#include "primerpair/primer_binding_site.hpp"

#include <limits>
#include <stdexcept>
#include <string>

namespace primerpair {

namespace {

char complement(
    const char base
) {
    switch (base) {
        case 'A':
            return 'T';

        case 'C':
            return 'G';

        case 'G':
            return 'C';

        case 'T':
            return 'A';

        case 'N':
            return 'N';

        default:
            throw std::invalid_argument(
                "Unsupported nucleotide in "
                "reverse complement."
            );
    }
}


std::string_view primer_sequence_for(
    const PrimerIdentity identity,
    const std::string_view primer1,
    const std::string_view primer2
) {
    switch (identity) {
        case PrimerIdentity::Primer1:
            return primer1;

        case PrimerIdentity::Primer2:
            return primer2;
    }

    throw std::logic_error(
        "Unknown primer identity."
    );
}


PrimerBindingSite extract_one_site(
    const PackedReference& reference,
    const PrimerIdentity identity,
    const std::uint64_t genomic_start,
    const bool reverse_strand,
    const std::uint64_t expected_mismatch_mask,
    const std::size_t expected_mismatch_count,
    const std::string_view primer1,
    const std::string_view primer2
) {
    const std::string_view primer =
        primer_sequence_for(
            identity,
            primer1,
            primer2
        );

    if (primer.empty()) {
        throw std::invalid_argument(
            "Primer sequence cannot be empty."
        );
    }

    /*
     * The mismatch representation is currently
     * uint64_t, therefore this layer must respect
     * the same maximum length.
     */
    if (primer.size() > 64) {
        throw std::invalid_argument(
            "Primer binding-site extraction "
            "supports at most 64 bases."
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
        static_cast<std::uint64_t>(
            primer.size()
        );

    site.primer_sequence =
        std::string(
            primer
        );

    site.genomic_sequence =
        packed_reference_subsequence(
            reference,
            genomic_start,
            primer.size()
        );

    if (reverse_strand) {
        site.binding_sequence =
            reverse_complement_binding_sequence(
                site.genomic_sequence
            );
    } else {
        site.binding_sequence =
            site.genomic_sequence;
    }


    std::uint64_t observed_mask = 0;

    for (
        std::size_t i = 0;
        i < primer.size();
        ++i
    ) {
        const char expected =
            primer.at(i);

        const char observed =
            site.binding_sequence.at(i);

        if (observed == 'N') {
            throw std::runtime_error(
                "Primer binding site contains N."
            );
        }

        if (expected == observed) {
            continue;
        }

        observed_mask |=
            (
                std::uint64_t{1}
                <<
                i
            );

        site.mismatches.push_back(
            PrimerBindingMismatch{
                i,
                expected,
                observed,
                primer.size() -
                    1 -
                    i
            }
        );
    }

    site.mismatch_mask =
        observed_mask;


    /*
     * Strong internal consistency check:
     *
     * Search-result mismatch information and the
     * actual PackedReference sequence must agree.
     */
    if (
        observed_mask !=
        expected_mismatch_mask
    ) {
        throw std::runtime_error(
            "Binding-site mismatch mask does not "
            "match PrimerPairHit."
        );
    }

    if (
        site.mismatch_count() !=
        expected_mismatch_count
    ) {
        throw std::runtime_error(
            "Binding-site mismatch count does not "
            "match PrimerPairHit."
        );
    }

    return site;
}

}  // namespace


std::string packed_reference_subsequence(
    const PackedReference& reference,
    const std::uint64_t start,
    const std::size_t length
) {
    if (
        start >
        reference.size()
    ) {
        throw std::out_of_range(
            "Binding-site start exceeds reference."
        );
    }

    if (
        static_cast<std::uint64_t>(
            length
        )
        >
        reference.size() -
            start
    ) {
        throw std::out_of_range(
            "Binding-site interval exceeds reference."
        );
    }

    std::string result;

    result.reserve(
        length
    );

    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        result.push_back(
            reference.base_at(
                start +
                static_cast<std::uint64_t>(
                    i
                )
            )
        );
    }

    return result;
}


std::string reverse_complement_binding_sequence(
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
        result.push_back(
            complement(
                *it
            )
        );
    }

    return result;
}


PrimerPairBindingSites extract_pair_binding_sites(
    const PackedReference& reference,
    const PrimerPairHit& hit,
    const std::string_view primer1,
    const std::string_view primer2
) {
    PrimerPairBindingSites result;

    /*
     * PCR-compatible pair representation:
     *
     * left hit  = forward-facing primer site
     * right hit = reverse-facing primer site
     */
    result.left =
        extract_one_site(
            reference,
            hit.left_primer,
            hit.left_position,
            false,
            hit.left_mismatch_mask,
            hit.left_mismatches,
            primer1,
            primer2
        );

    result.right =
        extract_one_site(
            reference,
            hit.right_primer,
            hit.right_position,
            true,
            hit.right_mismatch_mask,
            hit.right_mismatches,
            primer1,
            primer2
        );

    return result;
}

}  // namespace primerpair
