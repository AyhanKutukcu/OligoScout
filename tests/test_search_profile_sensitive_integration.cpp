#include "primerpair/search_profile.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
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

char mutate_base(
    const char base
) {
    return
        base == 'A'
            ? 'C'
            : 'A';
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

        /*
         * Forward hit with mismatch at ORIGINAL
         * primer position 19 = terminal 3-prime.
         */
        std::string forward_terminal =
            primer;

        forward_terminal.back() =
            mutate_base(
                forward_terminal.back()
            );

        /*
         * Reverse-orientation terminal 3-prime
         * mismatch.
         *
         * Original primer position 19 maps to
         * reverse_query position 0.
         */
        std::string reverse_terminal =
            rc;

        reverse_terminal.front() =
            mutate_base(
                reverse_terminal.front()
            );

        const std::string separator(
            20,
            'N'
        );

        std::string reference =
            separator;

        const std::uint64_t
            exact_forward_position =
                reference.size();

        reference += primer;
        reference += separator;

        const std::uint64_t
            exact_reverse_position =
                reference.size();

        reference += rc;
        reference += separator;

        const std::uint64_t
            forward_terminal_position =
                reference.size();

        reference +=
            forward_terminal;

        reference +=
            separator;

        const std::uint64_t
            reverse_terminal_position =
                reference.size();

        reference +=
            reverse_terminal;

        reference +=
            separator;

        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::
            SensitivePrimerSearchEngine
                exhaustive_reference(
                    index
                );

        const primerpair::
            ProfiledPrimerSearchEngine
                profile_engine(
                    index,
                    packed
                );

        /*
         * --------------------------------------------------
         * Differential equality k=0..3
         * --------------------------------------------------
         */

        for (
            std::size_t budget = 0;
            budget <= 3;
            ++budget
        ) {
            const auto reference_result =
                exhaustive_reference.search(
                    primer,
                    budget
                );

            const auto profiled_result =
                profile_engine.search(
                    primer,
                    primerpair::
                        SearchProfile::
                            Sensitive,
                    12,
                    budget
                );

            expect(
                profiled_result.profile ==
                    primerpair::
                        SearchProfile::
                            Sensitive,
                "SENSITIVE profile stored k=" +
                    std::to_string(
                        budget
                    )
            );

            expect(
                profiled_result
                    .sensitive_result_available,
                "SENSITIVE backend result available k=" +
                    std::to_string(
                        budget
                    )
            );

            expect(
                profiled_result.hits() ==
                    reference_result.hits,
                "SENSITIVE profile equals exhaustive reference k=" +
                    std::to_string(
                        budget
                    )
            );
        }

        /*
         * --------------------------------------------------
         * Biological profile difference at k=1
         * --------------------------------------------------
         */

        const auto strict =
            profile_engine.search(
                primer,
                primerpair::
                    SearchProfile::
                        Strict,
                12,
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

        /*
         * Exact forward and exact reverse must
         * remain present in SENSITIVE.
         */
        const auto exact_forward =
            std::find_if(
                sensitive.hits().begin(),
                sensitive.hits().end(),
                [
                    exact_forward_position
                ](
                    const auto& hit
                ) {
                    return
                        hit.position ==
                            exact_forward_position
                        &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                    Forward
                        &&
                        hit.mismatches == 0
                        &&
                        hit.mismatch_mask == 0;
                }
            );

        expect(
            exact_forward !=
                sensitive.hits().end(),
            "SENSITIVE preserves exact forward hit"
        );

        const auto exact_reverse =
            std::find_if(
                sensitive.hits().begin(),
                sensitive.hits().end(),
                [
                    exact_reverse_position
                ](
                    const auto& hit
                ) {
                    return
                        hit.position ==
                            exact_reverse_position
                        &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                    Reverse
                        &&
                        hit.mismatches == 0
                        &&
                        hit.mismatch_mask == 0;
                }
            );

        expect(
            exact_reverse !=
                sensitive.hits().end(),
            "SENSITIVE preserves exact reverse hit"
        );

        /*
         * Forward terminal 3-prime mismatch:
         *
         * bit 19 must be set.
         */
        const std::uint64_t
            terminal_mask =
                std::uint64_t{1}
                <<
                19;

        const auto forward_terminal_hit =
            std::find_if(
                sensitive.hits().begin(),
                sensitive.hits().end(),
                [
                    forward_terminal_position,
                    terminal_mask
                ](
                    const auto& hit
                ) {
                    return
                        hit.position ==
                            forward_terminal_position
                        &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                    Forward
                        &&
                        hit.mismatches == 1
                        &&
                        hit.mismatch_mask ==
                            terminal_mask;
                }
            );

        expect(
            forward_terminal_hit !=
                sensitive.hits().end(),
            "SENSITIVE recovers forward terminal 3-prime mismatch"
        );

        /*
         * Reverse terminal 3-prime mismatch must
         * also map to ORIGINAL position 19.
         */
        const auto reverse_terminal_hit =
            std::find_if(
                sensitive.hits().begin(),
                sensitive.hits().end(),
                [
                    reverse_terminal_position,
                    terminal_mask
                ](
                    const auto& hit
                ) {
                    return
                        hit.position ==
                            reverse_terminal_position
                        &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                    Reverse
                        &&
                        hit.mismatches == 1
                        &&
                        hit.mismatch_mask ==
                            terminal_mask;
                }
            );

        expect(
            reverse_terminal_hit !=
                sensitive.hits().end(),
            "SENSITIVE recovers reverse terminal 3-prime mismatch"
        );

        /*
         * STRICT exact 12-nt 3-prime anchor must
         * reject those same terminal mismatches.
         */
        const auto strict_forward_terminal =
            std::find_if(
                strict.hits().begin(),
                strict.hits().end(),
                [
                    forward_terminal_position
                ](
                    const auto& hit
                ) {
                    return
                        hit.position ==
                            forward_terminal_position
                        &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                    Forward;
                }
            );

        expect(
            strict_forward_terminal ==
                strict.hits().end(),
            "STRICT rejects forward terminal 3-prime mismatch"
        );

        const auto strict_reverse_terminal =
            std::find_if(
                strict.hits().begin(),
                strict.hits().end(),
                [
                    reverse_terminal_position
                ](
                    const auto& hit
                ) {
                    return
                        hit.position ==
                            reverse_terminal_position
                        &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                    Reverse;
                }
            );

        expect(
            strict_reverse_terminal ==
                strict.hits().end(),
            "STRICT rejects reverse terminal 3-prime mismatch"
        );

        /*
         * k<3 must stay on exhaustive backend.
         */
        expect(
            !sensitive
                .sensitive_result
                .estimator_used,
            "SENSITIVE k1 skips adaptive estimator"
        );

        /*
         * k=3 should invoke adaptive routing.
         */
        const auto sensitive_k3 =
            profile_engine.search(
                primer,
                primerpair::
                    SearchProfile::
                        Sensitive,
                12,
                3
            );

        expect(
            sensitive_k3
                .sensitive_result
                .estimator_used,
            "SENSITIVE k3 invokes adaptive estimator"
        );

        expect(
            sensitive_k3
                .sensitive_result
                .cost_estimate
                .max_seed_occurrences
            <=
            primerpair::
                kDefaultSensitiveK3MaxSeedThreshold,
            "Synthetic k3 lies below adaptive threshold"
        );

        expect(
            sensitive_k3
                .sensitive_result
                .backend
            ==
            primerpair::
                SensitiveAdaptiveBackend::
                    Candidate,
            "Synthetic k3 routed to candidate backend"
        );

        std::cout
            << "All SENSITIVE profile integration tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
