#include "primerpair/bidirectional_fm_index.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

std::vector<std::uint64_t> sorted(
    std::vector<std::uint64_t> values
) {
    std::sort(
        values.begin(),
        values.end()
    );

    return values;
}

bool same_interval(
    const primerpair::Interval& lhs,
    const primerpair::Interval& rhs
) {
    return
        lhs.begin == rhs.begin &&
        lhs.end == rhs.end;
}

void expect_same_state(
    const primerpair::BidirectionalInterval& lhs,
    const primerpair::BidirectionalInterval& rhs,
    const std::string& name
) {
    /*
     * Empty FM intervals için begin/end endpoint'i
     * semantik olarak anlamlı değildir.
     *
     * Örneğin {0,0} ve {17,17} aynı şeyi ifade eder:
     * sıfır eşleşme.
     *
     * Bu nedenle absent-pattern durumunda yalnızca
     * iki state'in de gerçekten empty olduğunu test ederiz.
     */
    if (lhs.empty() || rhs.empty()) {
        expect(
            lhs.empty() && rhs.empty(),
            name + " both states empty"
        );

        expect(
            lhs.forward.empty() &&
            lhs.reverse.empty(),
            name + " synchronized empty state"
        );

        expect(
            lhs.size() == 0 &&
            rhs.size() == 0,
            name + " zero occurrence count"
        );

        return;
    }

    /*
     * Non-empty state'lerde endpoint eşitliği zorunludur.
     * Gerçek BiFM senkronizasyonunu kanıtlayan kısım budur.
     */
    expect(
        same_interval(
            lhs.forward,
            rhs.forward
        ),
        name + " forward interval"
    );

    expect(
        same_interval(
            lhs.reverse,
            rhs.reverse
        ),
        name + " reverse interval"
    );

    expect(
        lhs.size() == rhs.size(),
        name + " occurrence count"
    );
}

}  // namespace

int main() {
    using primerpair::BidirectionalFMIndex;

    try {
        const std::string reference =
            "AACGTTAACGTAACGTAA";

        const BidirectionalFMIndex index(
            reference
        );

        const auto seed =
            index.search(
                "ACGT"
            );

        expect(
            seed.matched_length() == 4,
            "Seed matched length"
        );

        expect(
            seed.size() == 3,
            "Seed occurrence count"
        );

        expect(
            sorted(
                index.locate(seed)
            ) ==
                std::vector<std::uint64_t>{
                    1,
                    7,
                    12
                },
            "Seed genomic positions"
        );

        /*
         * --------------------------------------------------
         * LEFT EXTENSION EQUIVALENCE
         * --------------------------------------------------
         *
         * extend_left(search(P), c)
         *
         * mutlaka
         *
         * search(cP)
         *
         * ile aynı forward + reverse intervali vermeli.
         */
        constexpr std::array<char, 5>
            bases{
                'A',
                'C',
                'G',
                'N',
                'T'
            };

        for (const char base : bases) {
            std::string expected_pattern;
            expected_pattern.push_back(
                base
            );
            expected_pattern +=
                "ACGT";

            const auto extended =
                index.extend_left(
                    seed,
                    base
                );

            const auto direct =
                index.search(
                    expected_pattern
                );

            expect_same_state(
                extended,
                direct,
                std::string(
                    "Left extension "
                ) +
                    base
            );
        }

        /*
         * --------------------------------------------------
         * RIGHT EXTENSION EQUIVALENCE
         * --------------------------------------------------
         */
        for (const char base : bases) {
            std::string expected_pattern =
                "ACGT";

            expected_pattern.push_back(
                base
            );

            const auto extended =
                index.extend_right(
                    seed,
                    base
                );

            const auto direct =
                index.search(
                    expected_pattern
                );

            expect_same_state(
                extended,
                direct,
                std::string(
                    "Right extension "
                ) +
                    base
            );
        }

        /*
         * AACGT
         */
        const auto left =
            index.extend_left(
                seed,
                'A'
            );

        expect(
            left.matched_length() == 5,
            "Left-extension matched length"
        );

        expect(
            sorted(
                index.locate(left)
            ) ==
                std::vector<std::uint64_t>{
                    0,
                    6,
                    11
                },
            "Left-extension genomic positions"
        );

        /*
         * AACGTA
         */
        const auto both =
            index.extend_right(
                left,
                'A'
            );

        const auto direct_both =
            index.search(
                "AACGTA"
            );

        expect_same_state(
            both,
            direct_both,
            "Alternating left/right extension"
        );

        expect(
            both.matched_length() == 6,
            "Bidirectional matched length"
        );

        expect(
            sorted(
                index.locate(both)
            ) ==
                std::vector<std::uint64_t>{
                    6,
                    11
                },
            "Bidirectional genomic positions"
        );

        /*
         * Sentinel-sensitive test.
         *
         * Pattern reference başlangıcında da bulunduğu
         * için predecessor dağılımında '$' vardır.
         *
         * '$' lexicographic offset hesabına katılmazsa
         * reverse interval endpointleri burada bozulur.
         */
        const auto boundary_seed =
            index.search(
                "A"
            );

        for (const char base : bases) {
            std::string expected;
            expected.push_back(
                base
            );
            expected +=
                "A";

            expect_same_state(
                index.extend_left(
                    boundary_seed,
                    base
                ),
                index.search(
                    expected
                ),
                std::string(
                    "Sentinel-sensitive left "
                ) +
                    base
            );
        }


        /*
         * --------------------------------------------------
         * RANDOMIZED DIFFERENTIAL STRESS TEST
         * --------------------------------------------------
         *
         * Senkron BiFM extension sonuçlarını doğrudan
         * search() sonuçlarıyla karşılaştırır.
         */
        std::string stress_reference(
            4096,
            'A'
        );

        constexpr std::array<char, 4>
            dna{
                'A',
                'C',
                'G',
                'T'
            };

        std::uint64_t rng =
            0x9E3779B97F4A7C15ULL;

        const auto next_random =
            [&rng]() -> std::uint64_t {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;

                return rng;
            };

        for (std::size_t i = 0;
             i < stress_reference.size();
             ++i) {

            stress_reference.at(i) =
                dna.at(
                    next_random() & 3ULL
                );
        }

        const BidirectionalFMIndex stress_index(
            stress_reference
        );

        const auto semantically_equal =
            [](
                const primerpair::BidirectionalInterval& lhs,
                const primerpair::BidirectionalInterval& rhs
            ) -> bool {

                if (lhs.empty() || rhs.empty()) {
                    return
                        lhs.empty() &&
                        rhs.empty() &&
                        lhs.forward.empty() &&
                        lhs.reverse.empty() &&
                        lhs.size() == 0 &&
                        rhs.size() == 0;
                }

                return
                    lhs.forward.begin ==
                        rhs.forward.begin &&
                    lhs.forward.end ==
                        rhs.forward.end &&
                    lhs.reverse.begin ==
                        rhs.reverse.begin &&
                    lhs.reverse.end ==
                        rhs.reverse.end &&
                    lhs.size() ==
                        rhs.size();
            };

        bool randomized_single_step_ok = true;

        for (std::size_t trial = 0;
             trial < 1000;
             ++trial) {

            const std::size_t seed_length =
                1 +
                static_cast<std::size_t>(
                    next_random() % 12
                );

            const std::size_t available =
                stress_reference.size() -
                seed_length +
                1;

            const std::size_t position =
                static_cast<std::size_t>(
                    next_random() %
                    available
                );

            const std::string seed_pattern =
                stress_reference.substr(
                    position,
                    seed_length
                );

            const auto state =
                stress_index.search(
                    seed_pattern
                );

            for (const char base : bases) {

                std::string left_pattern;
                left_pattern.reserve(
                    seed_pattern.size() + 1
                );

                left_pattern.push_back(
                    base
                );

                left_pattern +=
                    seed_pattern;

                const auto left_extended =
                    stress_index.extend_left(
                        state,
                        base
                    );

                const auto left_direct =
                    stress_index.search(
                        left_pattern
                    );

                if (!semantically_equal(
                        left_extended,
                        left_direct
                    )) {

                    randomized_single_step_ok = false;
                    break;
                }

                std::string right_pattern =
                    seed_pattern;

                right_pattern.push_back(
                    base
                );

                const auto right_extended =
                    stress_index.extend_right(
                        state,
                        base
                    );

                const auto right_direct =
                    stress_index.search(
                        right_pattern
                    );

                if (!semantically_equal(
                        right_extended,
                        right_direct
                    )) {

                    randomized_single_step_ok = false;
                    break;
                }
            }

            if (!randomized_single_step_ok) {
                break;
            }
        }

        expect(
            randomized_single_step_ok,
            "1000 randomized differential extension trials"
        );

        /*
         * --------------------------------------------------
         * MULTI-STEP / ALTERNATING EXTENSION TEST
         * --------------------------------------------------
         *
         * Ortadaki seed'den başlayıp önce sola,
         * sonra sağa doğru gerçek substring yeniden
         * oluşturulur.
         */
        bool chained_extension_ok = true;

        constexpr std::size_t full_length = 20;
        constexpr std::size_t seed_begin = 7;
        constexpr std::size_t seed_length = 6;

        for (std::size_t trial = 0;
             trial < 250;
             ++trial) {

            const std::size_t available =
                stress_reference.size() -
                full_length +
                1;

            const std::size_t position =
                static_cast<std::size_t>(
                    next_random() %
                    available
                );

            const std::string full_pattern =
                stress_reference.substr(
                    position,
                    full_length
                );

            std::string current_pattern =
                full_pattern.substr(
                    seed_begin,
                    seed_length
                );

            auto current_state =
                stress_index.search(
                    current_pattern
                );

            /*
             * Sol taraf:
             * positions 6,5,...,0
             */
            for (int i =
                     static_cast<int>(
                         seed_begin
                     ) - 1;
                 i >= 0;
                 --i) {

                const char base =
                    full_pattern.at(
                        static_cast<std::size_t>(
                            i
                        )
                    );

                current_state =
                    stress_index.extend_left(
                        current_state,
                        base
                    );

                current_pattern.insert(
                    current_pattern.begin(),
                    base
                );

                const auto direct =
                    stress_index.search(
                        current_pattern
                    );

                if (!semantically_equal(
                        current_state,
                        direct
                    )) {

                    chained_extension_ok = false;
                    break;
                }
            }

            if (!chained_extension_ok) {
                break;
            }

            /*
             * Sağ taraf:
             * seed'in sağından 19'a kadar.
             */
            const std::size_t right_begin =
                seed_begin +
                seed_length;

            for (std::size_t i = right_begin;
                 i < full_length;
                 ++i) {

                const char base =
                    full_pattern.at(i);

                current_state =
                    stress_index.extend_right(
                        current_state,
                        base
                    );

                current_pattern.push_back(
                    base
                );

                const auto direct =
                    stress_index.search(
                        current_pattern
                    );

                if (!semantically_equal(
                        current_state,
                        direct
                    )) {

                    chained_extension_ok = false;
                    break;
                }
            }

            if (!chained_extension_ok) {
                break;
            }
        }

        expect(
            chained_extension_ok,
            "250 chained bidirectional extension trials"
        );

        const auto lowercase =
            index.search(
                "acgt"
            );

        expect(
            same_interval(
                lowercase.forward,
                seed.forward
            ),
            "Lowercase forward normalization"
        );

        expect(
            same_interval(
                lowercase.reverse,
                seed.reverse
            ),
            "Lowercase reverse normalization"
        );

        bool invalid_rejected = false;

        try {
            static_cast<void>(
                index.search(
                    "ACGX"
                )
            );
        } catch (
            const std::invalid_argument&
        ) {
            invalid_rejected = true;
        }

        expect(
            invalid_rejected,
            "Invalid nucleotide rejection"
        );

        expect(
            index.forward_index()
                .suffix_array_sample_rate() ==
                primerpair::FMIndex::
                    kDefaultSuffixArraySampleRate,
            "Forward default SA rate"
        );

        expect(
            index.reverse_index()
                .suffix_array_sample_rate() ==
                primerpair::FMIndex::
                    kDefaultSuffixArraySampleRate,
            "Reverse default SA rate"
        );

        std::cout
            << "All synchronized bidirectional FM-index tests passed.\n";

        return 0;

    } catch (const std::exception& error) {
        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
