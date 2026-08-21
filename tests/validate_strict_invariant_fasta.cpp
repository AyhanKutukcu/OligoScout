#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/mismatch_features.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace {

std::string load_fasta(
    const std::string& path
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open FASTA."
        );
    }

    std::string sequence;
    std::string line;

    while (std::getline(input, line)) {

        if (
            line.empty() ||
            line.front() == '>'
        ) {
            continue;
        }

        for (const char raw : line) {

            if (
                std::isspace(
                    static_cast<unsigned char>(raw)
                )
            ) {
                continue;
            }

            const char base =
                static_cast<char>(
                    std::toupper(
                        static_cast<unsigned char>(raw)
                    )
                );

            switch (base) {

                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':
                    sequence.push_back(base);
                    break;

                case 'R':
                case 'Y':
                case 'S':
                case 'W':
                case 'K':
                case 'M':
                case 'B':
                case 'D':
                case 'H':
                case 'V':
                    sequence.push_back('N');
                    break;

                default:
                    throw std::runtime_error(
                        "Unsupported FASTA nucleotide."
                    );
            }
        }
    }

    if (sequence.empty()) {
        throw std::runtime_error(
            "Empty FASTA."
        );
    }

    return sequence;
}

bool all_acgt(
    const std::string_view sequence
) {
    return std::all_of(
        sequence.begin(),
        sequence.end(),
        [](
            const char base
        ) {
            return
                base == 'A' ||
                base == 'C' ||
                base == 'G' ||
                base == 'T';
        }
    );
}

}  // namespace

int main(
    int argc,
    char* argv[]
) {
    try {
        if (argc < 2 || argc > 3) {

            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " [primer_count]\n";

            return 2;
        }

        constexpr std::size_t
            primer_length = 20;

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            approximate_prefix_length =
                primer_length -
                anchor_length;

        constexpr std::uint64_t
            stride = 104729ULL;

        const std::size_t primer_count =
            argc == 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 2048;

        if (primer_count == 0) {
            throw std::invalid_argument(
                "Primer count must be > 0."
            );
        }

        const std::string reference =
            load_fasta(
                argv[1]
            );

        if (
            reference.size() <
            primer_length
        ) {
            throw std::runtime_error(
                "Reference shorter than primer."
            );
        }

        /*
         * Deterministically sample real 20-mers
         * from chr22.
         */
        std::vector<std::string>
            primers;

        primers.reserve(
            primer_count
        );

        const std::uint64_t available =
            static_cast<std::uint64_t>(
                reference.size() -
                primer_length +
                1
            );

        std::uint64_t cursor = 0;

        for (
            std::size_t attempt = 0;
            attempt <
                reference.size() &&
            primers.size() <
                primer_count;
            ++attempt
        ) {
            cursor =
                (
                    cursor +
                    stride
                ) %
                available;

            const std::string_view candidate(
                reference.data() + cursor,
                primer_length
            );

            if (!all_acgt(candidate)) {
                continue;
            }

            primers.emplace_back(
                candidate
            );
        }

        if (
            primers.size() !=
            primer_count
        ) {
            throw std::runtime_error(
                "Could not collect requested primers."
            );
        }

        const primerpair::PackedReference
            packed(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::
            StrandAwarePrimerSearchEngine
                engine(
                    index,
                    packed
                );

        /*
         * For a 20-mer with a 12-nt exact 3-prime
         * anchor, only bits 0..7 may ever be set.
         */
        const std::uint64_t
            allowed_prefix_mask =
                (
                    std::uint64_t{1}
                    <<
                    approximate_prefix_length
                ) -
                1;

        std::cout
            << "reference_bp\t"
            << reference.size()
            << '\n';

        std::cout
            << "primer_count\t"
            << primers.size()
            << '\n';

        std::cout
            << "primer_length\t"
            << primer_length
            << '\n';

        std::cout
            << "anchor_length\t"
            << anchor_length
            << '\n';

        std::cout
            << "allowed_mismatch_positions\t0-"
            << (
                approximate_prefix_length -
                1
            )
            << '\n';

        std::cout
            << "budget"
            << '\t'
            << "total_hits"
            << '\t'
            << "forward_hits"
            << '\t'
            << "reverse_hits"
            << '\t'
            << "exact_hits"
            << '\t'
            << "mismatch_hits"
            << '\t'
            << "min_nearest_3prime"
            << '\t'
            << "max_nearest_3prime"
            << '\t'
            << "anchor_overlap_violations"
            << '\t'
            << "prefix_mask_violations"
            << '\t'
            << "count_mask_disagreements"
            << '\t'
            << "terminal_3prime_violations"
            << '\t'
            << "last3_violations"
            << '\t'
            << "last5_violations"
            << '\t'
            << "last8_violations"
            << '\t'
            << "budget_violations"
            << '\t'
            << "all_invariants"
            << '\n';

        for (
            std::size_t budget = 1;
            budget <= 3;
            ++budget
        ) {
            std::uint64_t total_hits = 0;
            std::uint64_t forward_hits = 0;
            std::uint64_t reverse_hits = 0;

            std::uint64_t exact_hits = 0;
            std::uint64_t mismatch_hits = 0;

            std::uint64_t
                anchor_overlap_violations = 0;

            std::uint64_t
                prefix_mask_violations = 0;

            std::uint64_t
                count_mask_disagreements = 0;

            std::uint64_t
                terminal_3prime_violations = 0;

            std::uint64_t
                last3_violations = 0;

            std::uint64_t
                last5_violations = 0;

            std::uint64_t
                last8_violations = 0;

            std::uint64_t
                budget_violations = 0;

            std::size_t
                min_nearest =
                    std::numeric_limits<
                        std::size_t
                    >::max();

            std::size_t
                max_nearest = 0;

            bool have_mismatch_hit =
                false;

            for (const auto& primer : primers) {

                const auto result =
                    engine.search(
                        primer,
                        anchor_length,
                        budget
                    );

                for (const auto& hit : result.hits) {

                    ++total_hits;

                    if (
                        hit.orientation ==
                        primerpair::
                            PrimerOrientation::
                            Forward
                    ) {
                        ++forward_hits;

                    } else {
                        ++reverse_hits;
                    }

                    const auto features =
                        primerpair::
                            extract_mismatch_features(
                                hit.mismatch_mask,
                                primer_length
                            );

                    if (
                        features.mismatch_count !=
                        hit.mismatches
                    ) {
                        ++count_mask_disagreements;
                    }

                    if (
                        hit.mismatches >
                        budget
                    ) {
                        ++budget_violations;
                    }

                    if (
                        hit.mismatches == 0
                    ) {
                        ++exact_hits;

                        if (
                            hit.mismatch_mask != 0 ||
                            features
                                .nearest_mismatch_to_3prime
                                .has_value()
                        ) {
                            ++count_mask_disagreements;
                        }

                    } else {

                        ++mismatch_hits;

                        have_mismatch_hit = true;

                        if (
                            !features
                                .nearest_mismatch_to_3prime
                                .has_value()
                        ) {
                            ++count_mask_disagreements;

                        } else {

                            const std::size_t nearest =
                                features
                                    .nearest_mismatch_to_3prime
                                    .value();

                            min_nearest =
                                std::min(
                                    min_nearest,
                                    nearest
                                );

                            max_nearest =
                                std::max(
                                    max_nearest,
                                    nearest
                                );

                            /*
                             * Exact 12-nt 3-prime anchor:
                             *
                             * nearest legal mismatch
                             * must be >= 12 nt away.
                             */
                            if (
                                nearest <
                                anchor_length
                            ) {
                                ++anchor_overlap_violations;
                            }
                        }
                    }

                    if (
                        primerpair::
                            mismatch_mask_overlaps_3prime_anchor(
                                hit.mismatch_mask,
                                primer_length,
                                anchor_length
                            )
                    ) {
                        ++anchor_overlap_violations;
                    }

                    /*
                     * Strong bit-level invariant:
                     *
                     * bits 8..19 must always be zero.
                     */
                    if (
                        (
                            hit.mismatch_mask &
                            ~allowed_prefix_mask
                        ) != 0
                    ) {
                        ++prefix_mask_violations;
                    }

                    if (
                        features
                            .terminal_3prime_mismatch
                    ) {
                        ++terminal_3prime_violations;
                    }

                    if (
                        features.last_3_count != 0
                    ) {
                        ++last3_violations;
                    }

                    if (
                        features.last_5_count != 0
                    ) {
                        ++last5_violations;
                    }

                    if (
                        features.last_8_count != 0
                    ) {
                        ++last8_violations;
                    }
                }
            }

            const bool all_invariants =
                anchor_overlap_violations == 0 &&
                prefix_mask_violations == 0 &&
                count_mask_disagreements == 0 &&
                terminal_3prime_violations == 0 &&
                last3_violations == 0 &&
                last5_violations == 0 &&
                last8_violations == 0 &&
                budget_violations == 0 &&
                (
                    !have_mismatch_hit ||
                    min_nearest >=
                        anchor_length
                );

            std::cout
                << budget
                << '\t'
                << total_hits
                << '\t'
                << forward_hits
                << '\t'
                << reverse_hits
                << '\t'
                << exact_hits
                << '\t'
                << mismatch_hits
                << '\t';

            if (have_mismatch_hit) {

                std::cout
                    << min_nearest
                    << '\t'
                    << max_nearest;

            } else {

                std::cout
                    << "NA"
                    << '\t'
                    << "NA";
            }

            std::cout
                << '\t'
                << anchor_overlap_violations
                << '\t'
                << prefix_mask_violations
                << '\t'
                << count_mask_disagreements
                << '\t'
                << terminal_3prime_violations
                << '\t'
                << last3_violations
                << '\t'
                << last5_violations
                << '\t'
                << last8_violations
                << '\t'
                << budget_violations
                << '\t'
                << (
                    all_invariants
                        ? "YES"
                        : "NO"
                )
                << '\n';

            if (!all_invariants) {
                throw std::logic_error(
                    "STRICT mismatch invariant violated."
                );
            }
        }

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
