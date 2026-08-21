#include "primerpair/sensitive_primer_search.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

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

        std::string reference(
            10,
            'N'
        );

        /*
         * --------------------------------------------------
         * Exact forward
         * --------------------------------------------------
         */
        const std::uint64_t
            exact_forward_position =
                reference.size();

        reference += primer;

        reference.append(
            20,
            'N'
        );

        /*
         * --------------------------------------------------
         * Exact reverse
         * --------------------------------------------------
         */
        const std::uint64_t
            exact_reverse_position =
                reference.size();

        reference += rc;

        reference.append(
            20,
            'N'
        );

        /*
         * --------------------------------------------------
         * Forward biological 3-prime terminal mismatch
         *
         * Original primer position 19.
         * --------------------------------------------------
         */
        std::string forward_terminal =
            primer;

        forward_terminal.back() =
            mutate_base(
                forward_terminal.back()
            );

        const std::uint64_t
            forward_terminal_position =
                reference.size();

        reference +=
            forward_terminal;

        reference.append(
            20,
            'N'
        );

        /*
         * --------------------------------------------------
         * Reverse biological 3-prime terminal mismatch
         *
         * In reverse-query coordinates, original
         * position 19 corresponds to query position 0.
         * --------------------------------------------------
         */
        std::string reverse_terminal =
            rc;

        reverse_terminal.front() =
            mutate_base(
                reverse_terminal.front()
            );

        const std::uint64_t
            reverse_terminal_position =
                reference.size();

        reference +=
            reverse_terminal;

        reference.append(
            20,
            'N'
        );

        /*
         * --------------------------------------------------
         * Reverse biological 5-prime mismatch
         *
         * reverse-query.back()
         * -> original primer position 0
         * --------------------------------------------------
         */
        std::string reverse_5prime =
            rc;

        reverse_5prime.back() =
            mutate_base(
                reverse_5prime.back()
            );

        const std::uint64_t
            reverse_5prime_position =
                reference.size();

        reference +=
            reverse_5prime;

        reference.append(
            10,
            'N'
        );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::
            SensitivePrimerSearchEngine
                engine(
                    index
                );

        /*
         * k=0: only exact forward/reverse.
         */
        const auto exact =
            engine.search(
                primer,
                0
            );

        expect(
            exact.hit_count() == 2,
            "Sensitive k0 exact both strands"
        );

        /*
         * k=1 must recover terminal mismatches that
         * STRICT deliberately rejects.
         */
        const auto sensitive =
            engine.search(
                primer,
                1
            );

        const auto find_hit =
            [
                &sensitive
            ](
                const std::uint64_t position,
                const primerpair::
                    PrimerOrientation orientation,
                const std::uint64_t mask
            ) {
                return
                    std::find_if(
                        sensitive.hits.begin(),
                        sensitive.hits.end(),
                        [
                            position,
                            orientation,
                            mask
                        ](
                            const auto& hit
                        ) {
                            return
                                hit.position ==
                                    position &&
                                hit.orientation ==
                                    orientation &&
                                hit.mismatches ==
                                    1 &&
                                hit.mismatch_mask ==
                                    mask;
                        }
                    ) !=
                    sensitive.hits.end();
            };

        expect(
            find_hit(
                forward_terminal_position,
                primerpair::
                    PrimerOrientation::
                    Forward,
                std::uint64_t{1} << 19
            ),
            "Forward 3-prime terminal mismatch recovered"
        );

        expect(
            find_hit(
                reverse_terminal_position,
                primerpair::
                    PrimerOrientation::
                    Reverse,
                std::uint64_t{1} << 19
            ),
            "Reverse 3-prime terminal mismatch maps to original position 19"
        );

        expect(
            find_hit(
                reverse_5prime_position,
                primerpair::
                    PrimerOrientation::
                    Reverse,
                std::uint64_t{1} << 0
            ),
            "Reverse 5-prime mismatch maps to original position 0"
        );

        /*
         * Exact hits still remain.
         */
        const auto exact_forward_it =
            std::find_if(
                sensitive.hits.begin(),
                sensitive.hits.end(),
                [
                    exact_forward_position
                ](
                    const auto& hit
                ) {
                    return
                        hit.position ==
                            exact_forward_position &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                Forward &&
                        hit.mismatches == 0 &&
                        hit.mismatch_mask == 0;
                }
            );

        expect(
            exact_forward_it !=
                sensitive.hits.end(),
            "Exact forward preserved in SENSITIVE"
        );

        const auto exact_reverse_it =
            std::find_if(
                sensitive.hits.begin(),
                sensitive.hits.end(),
                [
                    exact_reverse_position
                ](
                    const auto& hit
                ) {
                    return
                        hit.position ==
                            exact_reverse_position &&
                        hit.orientation ==
                            primerpair::
                                PrimerOrientation::
                                Reverse &&
                        hit.mismatches == 0 &&
                        hit.mismatch_mask == 0;
                }
            );

        expect(
            exact_reverse_it !=
                sensitive.hits.end(),
            "Exact reverse preserved in SENSITIVE"
        );

        /*
         * Lowercase query must behave identically.
         */
        std::string lowercase =
            primer;

        std::transform(
            lowercase.begin(),
            lowercase.end(),
            lowercase.begin(),
            [](
                const unsigned char base
            ) {
                return
                    static_cast<char>(
                        std::tolower(
                            base
                        )
                    );
            }
        );

        const auto lowercase_result =
            engine.search(
                lowercase,
                1
            );

        expect(
            lowercase_result.hits ==
                sensitive.hits,
            "Sensitive lowercase search"
        );

        bool excessive_budget_rejected =
            false;

        try {
            static_cast<void>(
                engine.search(
                    primer,
                    4
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            excessive_budget_rejected =
                true;
        }

        expect(
            excessive_budget_rejected,
            "Sensitive mismatch budget >3 rejected"
        );

        std::cout
            << "sensitive_k0_hits\t"
            << exact.hit_count()
            << '\n';

        std::cout
            << "sensitive_k1_hits\t"
            << sensitive.hit_count()
            << '\n';

        std::cout
            << "All sensitive primer search tests passed.\n";

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
