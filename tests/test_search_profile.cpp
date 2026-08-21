#include "primerpair/search_profile.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>

#include "primerpair/sensitive_primer_search.hpp"

namespace {

void expect(
    const bool condition,
    const std::string& name
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " + name
        );
    }

    std::cout
        << "[PASS] "
        << name
        << '\n';
}

}  // namespace


int main() {
    try {
        const std::string primer =
            "ACGTTGCAAGTCCTGAACGA";

        const std::string rc =
            primerpair::
                reverse_complement(
                    primer
                );

        std::string reverse_mismatch =
            rc;

        reverse_mismatch.back() =
            reverse_mismatch.back() == 'A'
                ? 'C'
                : 'A';

        std::string reference =
            "NNNNNNNNNN";

        reference += primer;

        reference +=
            "NNNNNNNNNN";

        reference += rc;

        reference +=
            "NNNNNNNNNN";

        reference +=
            reverse_mismatch;

        reference +=
            "NNNNNNNNNN";

        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        /*
         * Existing validated STRICT engine.
         */
        const primerpair::
            StrandAwarePrimerSearchEngine
                legacy_engine(
                    index,
                    packed
                );

        /*
         * Exhaustive SENSITIVE correctness
         * reference.
         */
        const primerpair::
            SensitivePrimerSearchEngine
                sensitive_reference(
                    index
                );

        /*
         * Profile façade.
         */
        const primerpair::
            ProfiledPrimerSearchEngine
                profile_engine(
                    index,
                    packed
                );

        /*
         * --------------------------------------------------
         * STRICT regression
         * --------------------------------------------------
         */

        const auto legacy =
            legacy_engine.search(
                primer,
                12,
                1
            );

        const auto strict =
            profile_engine.search(
                primer,
                primerpair::
                    SearchProfile::
                        Strict,
                12,
                1
            );

        expect(
            strict.profile ==
                primerpair::
                    SearchProfile::
                        Strict,
            "Explicit STRICT profile stored"
        );

        expect(
            !strict
                .sensitive_result_available,
            "STRICT does not expose sensitive result"
        );

        expect(
            strict.search_result
                .primer_length ==
            legacy.primer_length,
            "STRICT primer length unchanged"
        );

        expect(
            strict.search_result.hits ==
                legacy.hits,
            "STRICT hits identical to legacy engine"
        );

        expect(
            strict.search_result.hit_count() ==
                legacy.hit_count(),
            "STRICT hit count identical"
        );

        expect(
            strict.hits() ==
                legacy.hits,
            "Generic STRICT hits accessor"
        );

        expect(
            strict.hit_count() ==
                legacy.hit_count(),
            "Generic STRICT hit count accessor"
        );

        expect(
            strict.primer_length() ==
                legacy.primer_length,
            "Generic STRICT primer length accessor"
        );

        /*
         * --------------------------------------------------
         * Profile names
         * --------------------------------------------------
         */

        expect(
            std::string(
                primerpair::to_string(
                    primerpair::
                        SearchProfile::
                            Strict
                )
            ) ==
            "STRICT",
            "STRICT profile name"
        );

        expect(
            std::string(
                primerpair::to_string(
                    primerpair::
                        SearchProfile::
                            Sensitive
                )
            ) ==
            "SENSITIVE",
            "SENSITIVE profile name"
        );

        /*
         * --------------------------------------------------
         * Basic SENSITIVE façade integration
         * --------------------------------------------------
         */

        const auto reference_sensitive =
            sensitive_reference.search(
                primer,
                1
            );

        const auto sensitive =
            profile_engine.search(
                primer,
                primerpair::
                    SearchProfile::
                        Sensitive,
                12,
                1
            );

        expect(
            sensitive.profile ==
                primerpair::
                    SearchProfile::
                        Sensitive,
            "Explicit SENSITIVE profile stored"
        );

        expect(
            sensitive
                .sensitive_result_available,
            "SENSITIVE result available"
        );

        expect(
            sensitive.hits() ==
                reference_sensitive.hits,
            "SENSITIVE façade equals exhaustive reference"
        );

        expect(
            sensitive.hit_count() ==
                reference_sensitive.hit_count(),
            "SENSITIVE façade hit count"
        );

        expect(
            sensitive.primer_length() ==
                reference_sensitive.primer_length,
            "SENSITIVE façade primer length"
        );

        std::cout
            << "All search-profile tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
