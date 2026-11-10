#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "primerpair/mismatch_features.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/sensitive_primer_search.hpp"
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
                    static_cast<unsigned char>(
                        raw
                    )
                )
            ) {
                continue;
            }

            const char base =
                static_cast<char>(
                    std::toupper(
                        static_cast<unsigned char>(
                            raw
                        )
                    )
                );

            switch (base) {

                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':

                    sequence.push_back(
                        base
                    );

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

                    sequence.push_back(
                        'N'
                    );

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

bool hit_less(
    const primerpair::
        OrientedPrimerSearchHit& lhs,
    const primerpair::
        OrientedPrimerSearchHit& rhs
) noexcept {
    if (
        lhs.position !=
        rhs.position
    ) {
        return
            lhs.position <
            rhs.position;
    }

    if (
        lhs.orientation !=
        rhs.orientation
    ) {
        return
            static_cast<int>(
                lhs.orientation
            ) <
            static_cast<int>(
                rhs.orientation
            );
    }

    if (
        lhs.mismatches !=
        rhs.mismatches
    ) {
        return
            lhs.mismatches <
            rhs.mismatches;
    }

    return
        lhs.mismatch_mask <
        rhs.mismatch_mask;
}

double percentile(
    std::vector<double> values,
    const double fraction
) {
    if (values.empty()) {
        return 0.0;
    }

    std::sort(
        values.begin(),
        values.end()
    );

    const double scaled =
        fraction *
        static_cast<double>(
            values.size() - 1
        );

    const std::size_t index =
        static_cast<std::size_t>(
            scaled
        );

    return
        values.at(
            index
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

        constexpr std::uint64_t
            stride = 104729ULL;

        const std::size_t primer_count =
            argc == 3
                ? static_cast<std::size_t>(
                      std::stoull(
                          argv[2]
                      )
                  )
                : 256;

        if (primer_count == 0) {
            throw std::invalid_argument(
                "Primer count must be > 0."
            );
        }

        const std::string reference =
            load_fasta(
                argv[1]
            );

        const std::uint64_t available =
            static_cast<std::uint64_t>(
                reference.size() -
                primer_length +
                1
            );

        /*
         * --------------------------------------------------
         * Collect deterministic UNIQUE real chr22 primers.
         * --------------------------------------------------
         */

        std::vector<std::string>
            primers;

        primers.reserve(
            primer_count
        );

        std::unordered_set<std::string>
            seen;

        seen.reserve(
            primer_count * 2
        );

        std::uint64_t cursor = 0;

        for (
            std::size_t attempt = 0;
            attempt < reference.size() &&
            primers.size() < primer_count;
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

            const std::string primer(
                candidate
            );

            if (
                !seen.insert(
                    primer
                ).second
            ) {
                continue;
            }

            primers.push_back(
                primer
            );
        }

        if (
            primers.size() !=
            primer_count
        ) {
            throw std::runtime_error(
                "Could not collect requested "
                "unique primers."
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
                strict_engine(
                    index,
                    packed
                );

        const primerpair::
            SensitivePrimerSearchEngine
                sensitive_engine(
                    index
                );

        std::cout
            << std::fixed
            << std::setprecision(3);

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
            << "strict_anchor_length\t"
            << anchor_length
            << '\n';

        std::cout
            << "budget"
            << '\t'
            << "strict_hits"
            << '\t'
            << "sensitive_hits"
            << '\t'
            << "additional_hits"
            << '\t'
            << "sensitive_vs_strict_hits"
            << '\t'
            << "additional_pct_sensitive"
            << '\t'
            << "additional_anchor_overlap"
            << '\t'
            << "additional_without_anchor_overlap"
            << '\t'
            << "additional_terminal_3prime"
            << '\t'
            << "additional_last3"
            << '\t'
            << "additional_last5"
            << '\t'
            << "additional_last8"
            << '\t'
            << "additional_min_nearest_3prime"
            << '\t'
            << "additional_max_nearest_3prime"
            << '\t'
            << "strict_subset_violations"
            << '\t'
            << "strict_median_us"
            << '\t'
            << "sensitive_median_us"
            << '\t'
            << "runtime_ratio"
            << '\t'
            << "strict_p95_us"
            << '\t'
            << "sensitive_p95_us"
            << '\t'
            << "all_checks"
            << '\n';

        for (
            std::size_t budget = 0;
            budget <= 3;
            ++budget
        ) {
            std::uint64_t
                strict_hits_total = 0;

            std::uint64_t
                sensitive_hits_total = 0;

            std::uint64_t
                additional_total = 0;

            std::uint64_t
                additional_anchor_overlap = 0;

            std::uint64_t
                additional_without_anchor_overlap = 0;

            std::uint64_t
                additional_terminal = 0;

            std::uint64_t
                additional_last3 = 0;

            std::uint64_t
                additional_last5 = 0;

            std::uint64_t
                additional_last8 = 0;

            std::uint64_t
                strict_subset_violations = 0;

            std::size_t
                additional_min_nearest =
                    std::numeric_limits<
                        std::size_t
                    >::max();

            std::size_t
                additional_max_nearest = 0;

            bool have_additional =
                false;

            std::vector<double>
                strict_times_us;

            std::vector<double>
                sensitive_times_us;

            strict_times_us.reserve(
                primers.size()
            );

            sensitive_times_us.reserve(
                primers.size()
            );

            for (
                std::size_t primer_index = 0;
                primer_index < primers.size();
                ++primer_index
            ) {
                const auto& primer =
                    primers.at(
                        primer_index
                    );

                primerpair::
                    StrandAwarePrimerSearchResult
                        strict;

                primerpair::
                    SensitivePrimerSearchResult
                        sensitive;

                double strict_us = 0.0;
                double sensitive_us = 0.0;

                /*
                 * Alternate engine order to reduce
                 * systematic cache/order bias.
                 */
                if (
                    (
                        primer_index +
                        budget
                    ) %
                    2 ==
                    0
                ) {
                    {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        strict =
                            strict_engine.search(
                                primer,
                                anchor_length,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        strict_us =
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0;
                    }

                    {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        sensitive =
                            sensitive_engine.search(
                                primer,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        sensitive_us =
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0;
                    }

                } else {

                    {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        sensitive =
                            sensitive_engine.search(
                                primer,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        sensitive_us =
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0;
                    }

                    {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        strict =
                            strict_engine.search(
                                primer,
                                anchor_length,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        strict_us =
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0;
                    }
                }

                strict_times_us.push_back(
                    strict_us
                );

                sensitive_times_us.push_back(
                    sensitive_us
                );

                strict_hits_total +=
                    strict.hit_count();

                sensitive_hits_total +=
                    sensitive.hit_count();

                /*
                 * Both engines are expected to use
                 * the same deterministic ordering.
                 */
                if (
                    !std::is_sorted(
                        strict.hits.begin(),
                        strict.hits.end(),
                        hit_less
                    ) ||
                    !std::is_sorted(
                        sensitive.hits.begin(),
                        sensitive.hits.end(),
                        hit_less
                    )
                ) {
                    throw std::logic_error(
                        "Hit vectors are not sorted."
                    );
                }

                /*
                 * Fundamental correctness condition:
                 *
                 * SENSITIVE must be a superset of STRICT.
                 */
                if (
                    !std::includes(
                        sensitive.hits.begin(),
                        sensitive.hits.end(),

                        strict.hits.begin(),
                        strict.hits.end(),

                        hit_less
                    )
                ) {
                    ++strict_subset_violations;

                    continue;
                }

                std::vector<
                    primerpair::
                        OrientedPrimerSearchHit
                > additional;

                additional.reserve(
                    sensitive.hits.size() -
                    strict.hits.size()
                );

                std::set_difference(
                    sensitive.hits.begin(),
                    sensitive.hits.end(),

                    strict.hits.begin(),
                    strict.hits.end(),

                    std::back_inserter(
                        additional
                    ),

                    hit_less
                );

                additional_total +=
                    additional.size();

                for (
                    const auto& hit :
                    additional
                ) {
                    const auto features =
                        primerpair::
                            extract_mismatch_features(
                                hit.mismatch_mask,
                                primer_length
                            );

                    /*
                     * Every additional SENSITIVE hit
                     * should differ from STRICT because
                     * at least one mismatch lies inside
                     * the exact 12-nt 3-prime anchor.
                     */
                    const bool anchor_overlap =
                        primerpair::
                            mismatch_mask_overlaps_3prime_anchor(
                                hit.mismatch_mask,
                                primer_length,
                                anchor_length
                            );

                    if (anchor_overlap) {

                        ++additional_anchor_overlap;

                    } else {

                        ++additional_without_anchor_overlap;
                    }

                    if (
                        features
                            .terminal_3prime_mismatch
                    ) {
                        ++additional_terminal;
                    }

                    if (
                        features.last_3_count > 0
                    ) {
                        ++additional_last3;
                    }

                    if (
                        features.last_5_count > 0
                    ) {
                        ++additional_last5;
                    }

                    if (
                        features.last_8_count > 0
                    ) {
                        ++additional_last8;
                    }

                    if (
                        features
                            .nearest_mismatch_to_3prime
                            .has_value()
                    ) {
                        have_additional = true;

                        const std::size_t nearest =
                            features
                                .nearest_mismatch_to_3prime
                                .value();

                        additional_min_nearest =
                            std::min(
                                additional_min_nearest,
                                nearest
                            );

                        additional_max_nearest =
                            std::max(
                                additional_max_nearest,
                                nearest
                            );
                    }
                }
            }

            const double strict_median =
                percentile(
                    strict_times_us,
                    0.50
                );

            const double sensitive_median =
                percentile(
                    sensitive_times_us,
                    0.50
                );

            const double strict_p95 =
                percentile(
                    strict_times_us,
                    0.95
                );

            const double sensitive_p95 =
                percentile(
                    sensitive_times_us,
                    0.95
                );

            const double hit_ratio =
                strict_hits_total > 0
                    ? static_cast<double>(
                          sensitive_hits_total
                      ) /
                      static_cast<double>(
                          strict_hits_total
                      )
                    : 0.0;

            const double additional_pct =
                sensitive_hits_total > 0
                    ? 100.0 *
                      static_cast<double>(
                          additional_total
                      ) /
                      static_cast<double>(
                          sensitive_hits_total
                      )
                    : 0.0;

            const double runtime_ratio =
                strict_median > 0.0
                    ? sensitive_median /
                      strict_median
                    : 0.0;

            /*
             * k=0 should be EXACTLY identical.
             *
             * k>0:
             * every extra hit must overlap the
             * STRICT 3-prime anchor.
             */
            const bool all_checks =
                strict_subset_violations == 0 &&
                additional_without_anchor_overlap == 0 &&
                (
                    budget != 0 ||
                    additional_total == 0
                ) &&
                sensitive_hits_total >=
                    strict_hits_total;

            std::cout
                << budget
                << '\t'
                << strict_hits_total
                << '\t'
                << sensitive_hits_total
                << '\t'
                << additional_total
                << '\t'
                << hit_ratio
                << '\t'
                << additional_pct
                << '\t'
                << additional_anchor_overlap
                << '\t'
                << additional_without_anchor_overlap
                << '\t'
                << additional_terminal
                << '\t'
                << additional_last3
                << '\t'
                << additional_last5
                << '\t'
                << additional_last8
                << '\t';

            if (have_additional) {

                std::cout
                    << additional_min_nearest
                    << '\t'
                    << additional_max_nearest;

            } else {

                std::cout
                    << "NA"
                    << '\t'
                    << "NA";
            }

            std::cout
                << '\t'
                << strict_subset_violations
                << '\t'
                << strict_median
                << '\t'
                << sensitive_median
                << '\t'
                << runtime_ratio
                << '\t'
                << strict_p95
                << '\t'
                << sensitive_p95
                << '\t'
                << (
                    all_checks
                        ? "YES"
                        : "NO"
                )
                << '\n';

            if (!all_checks) {
                throw std::logic_error(
                    "STRICT-vs-SENSITIVE "
                    "correctness check failed."
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
