#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"

namespace {

struct DesignedPair {
    std::string primer1;
    std::string primer2;

    std::uint64_t left_position{0};
    std::uint64_t right_position{0};
    std::uint64_t expected_amplicon_length{0};
};


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
                    sequence.push_back(base);
                    break;

                case 'N':
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


bool product_less(
    const primerpair::PrimerPairHit& lhs,
    const primerpair::PrimerPairHit& rhs
) {
    return std::tie(
        lhs.amplicon_start,
        lhs.amplicon_end_exclusive,
        lhs.left_primer,
        lhs.right_primer,
        lhs.left_position,
        lhs.right_position,
        lhs.left_mismatches,
        lhs.right_mismatches,
        lhs.left_mismatch_mask,
        lhs.right_mismatch_mask
    )
    <
    std::tie(
        rhs.amplicon_start,
        rhs.amplicon_end_exclusive,
        rhs.left_primer,
        rhs.right_primer,
        rhs.left_position,
        rhs.right_position,
        rhs.left_mismatches,
        rhs.right_mismatches,
        rhs.left_mismatch_mask,
        rhs.right_mismatch_mask
    );
}


std::vector<primerpair::PrimerPairHit>
canonical_products(
    const primerpair::PrimerPairSearchResult& result
) {
    auto products =
        result.amplicons;

    std::sort(
        products.begin(),
        products.end(),
        product_less
    );

    products.erase(
        std::unique(
            products.begin(),
            products.end()
        ),
        products.end()
    );

    return products;
}


bool contains_designed_target(
    const primerpair::PrimerPairSearchResult& result,
    const DesignedPair& pair
) {
    return std::any_of(
        result.amplicons.begin(),
        result.amplicons.end(),
        [
            &pair
        ](
            const auto& hit
        ) {
            return
                hit.left_primer ==
                    primerpair::
                        PrimerIdentity::
                            Primer1
                &&
                hit.right_primer ==
                    primerpair::
                        PrimerIdentity::
                            Primer2
                &&
                hit.left_position ==
                    pair.left_position
                &&
                hit.right_position ==
                    pair.right_position
                &&
                hit.amplicon_length ==
                    pair.expected_amplicon_length;
        }
    );
}


bool mask_has_last_n(
    const std::uint64_t mask,
    const std::size_t primer_length,
    const std::size_t n
) {
    if (
        primer_length == 0 ||
        n == 0
    ) {
        return false;
    }

    const std::size_t width =
        std::min(
            primer_length,
            n
        );

    const std::size_t start =
        primer_length - width;

    for (
        std::size_t position = start;
        position < primer_length;
        ++position
    ) {
        if (
            (
                mask &
                (
                    std::uint64_t{1}
                    <<
                    position
                )
            ) != 0
        ) {
            return true;
        }
    }

    return false;
}


bool product_has_last_n_mismatch(
    const primerpair::PrimerPairHit& hit,
    const std::size_t primer1_length,
    const std::size_t primer2_length,
    const std::size_t n
) {
    const std::size_t left_length =
        hit.left_primer ==
            primerpair::PrimerIdentity::Primer1
            ? primer1_length
            : primer2_length;

    const std::size_t right_length =
        hit.right_primer ==
            primerpair::PrimerIdentity::Primer1
            ? primer1_length
            : primer2_length;

    return
        mask_has_last_n(
            hit.left_mismatch_mask,
            left_length,
            n
        )
        ||
        mask_has_last_n(
            hit.right_mismatch_mask,
            right_length,
            n
        );
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

    const double position =
        fraction *
        static_cast<double>(
            values.size() - 1
        );

    const auto lower =
        static_cast<std::size_t>(
            position
        );

    const auto upper =
        std::min(
            lower + 1,
            values.size() - 1
        );

    const double weight =
        position -
        static_cast<double>(
            lower
        );

    return
        values.at(lower) *
            (1.0 - weight)
        +
        values.at(upper) *
            weight;
}


double mean(
    const std::vector<double>& values
) {
    if (values.empty()) {
        return 0.0;
    }

    long double total = 0.0L;

    for (const double value : values) {
        total += value;
    }

    return static_cast<double>(
        total /
        static_cast<long double>(
            values.size()
        )
    );
}

}  // namespace


int main(
    int argc,
    char* argv[]
) {
    try {
        if (
            argc < 2 ||
            argc > 4
        ) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " [pair_count]"
                << " [repeats]\n";

            return 2;
        }

        constexpr std::size_t
            primer_length = 20;

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::uint64_t
            min_amplicon = 50;

        constexpr std::uint64_t
            max_amplicon = 3000;

        constexpr std::uint64_t
            sampling_stride = 130363ULL;

        constexpr std::uint64_t
            initial_offset = 50021ULL;

        const std::size_t target_pairs =
            argc >= 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 256;

        const std::size_t repeats =
            argc >= 4
                ? static_cast<std::size_t>(
                      std::stoull(argv[3])
                  )
                : 3;

        if (
            target_pairs == 0 ||
            repeats == 0
        ) {
            throw std::invalid_argument(
                "pair_count and repeats must be > 0."
            );
        }

        const std::string reference =
            load_fasta(
                argv[1]
            );

        if (
            reference.size() <=
                max_amplicon +
                primer_length
        ) {
            throw std::runtime_error(
                "Reference is too short."
            );
        }

        /*
         * --------------------------------------------------
         * Deterministic real-reference pair sampling
         * --------------------------------------------------
         */

        std::vector<DesignedPair>
            pairs;

        pairs.reserve(
            target_pairs
        );

        std::unordered_set<std::string>
            seen_pairs;

        seen_pairs.reserve(
            target_pairs * 2
        );

        const std::uint64_t scan_span =
            static_cast<std::uint64_t>(
                reference.size()
            ) -
            max_amplicon -
            1;

        std::uint64_t cursor =
            initial_offset %
            scan_span;

        const std::size_t max_attempts =
            target_pairs *
            100000;

        for (
            std::size_t attempt = 0;
            attempt < max_attempts &&
            pairs.size() < target_pairs;
            ++attempt
        ) {
            cursor =
                (
                    cursor +
                    sampling_stride
                ) %
                scan_span;

            const std::uint64_t
                amplicon_length =
                    min_amplicon +
                    (
                        (
                            pairs.size() * 97ULL +
                            attempt * 13ULL
                        )
                        %
                        (
                            max_amplicon -
                            min_amplicon +
                            1
                        )
                    );

            const std::uint64_t
                right_position =
                    cursor +
                    amplicon_length -
                    primer_length;

            const std::string_view
                left_site(
                    reference.data() +
                        cursor,
                    primer_length
                );

            const std::string_view
                right_site(
                    reference.data() +
                        right_position,
                    primer_length
                );

            if (
                !all_acgt(
                    left_site
                )
                ||
                !all_acgt(
                    right_site
                )
            ) {
                continue;
            }

            const std::string primer1(
                left_site
            );

            const std::string primer2 =
                primerpair::
                    reverse_complement(
                        right_site
                    );

            if (primer1 == primer2) {
                continue;
            }

            const std::string key =
                primer1 +
                "\t" +
                primer2;

            if (
                !seen_pairs
                    .insert(key)
                    .second
            ) {
                continue;
            }

            pairs.push_back(
                DesignedPair{
                    primer1,
                    primer2,
                    cursor,
                    right_position,
                    amplicon_length
                }
            );
        }

        if (
            pairs.size() !=
            target_pairs
        ) {
            throw std::runtime_error(
                "Could not collect requested "
                "number of primer pairs."
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

        const primerpair::PrimerPairSearchEngine
            pair_engine(
                index,
                packed
            );

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "# reference_bp\t"
            << reference.size()
            << '\n';

        std::cout
            << "# designed_pairs\t"
            << pairs.size()
            << '\n';

        std::cout
            << "# primer_length\t"
            << primer_length
            << '\n';

        std::cout
            << "# anchor_length\t"
            << anchor_length
            << '\n';

        std::cout
            << "# amplicon_range\t"
            << min_amplicon
            << '-'
            << max_amplicon
            << '\n';

        std::cout
            << "# repeats\t"
            << repeats
            << '\n';

        std::cout
            << "budget"
            << '\t'
            << "strict_single_total"
            << '\t'
            << "sensitive_single_total"
            << '\t'
            << "single_hit_expansion"
            << '\t'
            << "strict_naive_cartesian"
            << '\t'
            << "sensitive_naive_cartesian"
            << '\t'
            << "strict_amplicons"
            << '\t'
            << "sensitive_amplicons"
            << '\t'
            << "additional_amplicons"
            << '\t'
            << "pair_amplicon_expansion"
            << '\t'
            << "expansion_compression"
            << '\t'
            << "strict_pair_reduction"
            << '\t'
            << "sensitive_pair_reduction"
            << '\t'
            << "add_anchor12"
            << '\t'
            << "add_last8"
            << '\t'
            << "add_last5"
            << '\t'
            << "add_last3"
            << '\t'
            << "add_terminal3"
            << '\t'
            << "add_mm0"
            << '\t'
            << "add_mm1"
            << '\t'
            << "add_mm2"
            << '\t'
            << "add_mm3"
            << '\t'
            << "add_mm4"
            << '\t'
            << "add_mm5"
            << '\t'
            << "add_mm6"
            << '\t'
            << "strict_targets"
            << '\t'
            << "sensitive_targets"
            << '\t'
            << "subset_violations"
            << '\t'
            << "k0_equality_violations"
            << '\t'
            << "anchor_explanation_violations"
            << '\t'
            << "strict_mean_us"
            << '\t'
            << "sensitive_mean_us"
            << '\t'
            << "runtime_ratio"
            << '\t'
            << "strict_median_us"
            << '\t'
            << "sensitive_median_us"
            << '\t'
            << "strict_p95_us"
            << '\t'
            << "sensitive_p95_us"
            << '\t'
            << "all_checks"
            << '\n';

        constexpr std::size_t budgets[] = {
            0,
            1,
            2,
            3
        };

        for (
            const std::size_t budget :
            budgets
        ) {
            long double
                strict_single_total = 0.0L;

            long double
                sensitive_single_total = 0.0L;

            long double
                strict_naive_total = 0.0L;

            long double
                sensitive_naive_total = 0.0L;

            std::uint64_t
                strict_amplicon_total = 0;

            std::uint64_t
                sensitive_amplicon_total = 0;

            std::uint64_t
                additional_total = 0;

            std::uint64_t
                add_anchor12 = 0;

            std::uint64_t
                add_last8 = 0;

            std::uint64_t
                add_last5 = 0;

            std::uint64_t
                add_last3 = 0;

            std::uint64_t
                add_terminal3 = 0;

            std::array<
                std::uint64_t,
                7
            > mismatch_distribution{};

            std::size_t
                strict_targets = 0;

            std::size_t
                sensitive_targets = 0;

            std::size_t
                subset_violations = 0;

            std::size_t
                equality_violations = 0;

            std::size_t
                anchor_explanation_violations = 0;

            std::vector<double>
                strict_pair_medians;

            std::vector<double>
                sensitive_pair_medians;

            strict_pair_medians.reserve(
                pairs.size()
            );

            sensitive_pair_medians.reserve(
                pairs.size()
            );

            for (
                std::size_t pair_index = 0;
                pair_index < pairs.size();
                ++pair_index
            ) {
                const auto& pair =
                    pairs.at(
                        pair_index
                    );

                std::vector<double>
                    strict_times;

                std::vector<double>
                    sensitive_times;

                strict_times.reserve(
                    repeats
                );

                sensitive_times.reserve(
                    repeats
                );

                primerpair::
                    PrimerPairSearchResult
                        strict_result;

                primerpair::
                    PrimerPairSearchResult
                        sensitive_result;

                for (
                    std::size_t repeat = 0;
                    repeat < repeats;
                    ++repeat
                ) {
                    auto run_strict =
                        [&]() {
                            const auto start =
                                std::chrono::
                                    steady_clock::now();

                            strict_result =
                                pair_engine.search(
                                    pair.primer1,
                                    pair.primer2,
                                    primerpair::
                                        SearchProfile::
                                            Strict,
                                    anchor_length,
                                    budget,
                                    min_amplicon,
                                    max_amplicon
                                );

                            const auto stop =
                                std::chrono::
                                    steady_clock::now();

                            strict_times.push_back(
                                static_cast<double>(
                                    std::chrono::
                                        duration_cast<
                                            std::chrono::
                                                nanoseconds
                                        >(
                                            stop - start
                                        ).count()
                                )
                                /
                                1000.0
                            );
                        };

                    auto run_sensitive =
                        [&]() {
                            const auto start =
                                std::chrono::
                                    steady_clock::now();

                            sensitive_result =
                                pair_engine.search(
                                    pair.primer1,
                                    pair.primer2,
                                    primerpair::
                                        SearchProfile::
                                            Sensitive,
                                    anchor_length,
                                    budget,
                                    min_amplicon,
                                    max_amplicon
                                );

                            const auto stop =
                                std::chrono::
                                    steady_clock::now();

                            sensitive_times.push_back(
                                static_cast<double>(
                                    std::chrono::
                                        duration_cast<
                                            std::chrono::
                                                nanoseconds
                                        >(
                                            stop - start
                                        ).count()
                                )
                                /
                                1000.0
                            );
                        };

                    /*
                     * Alternate execution order to reduce
                     * systematic cache/order bias.
                     */
                    if (
                        (
                            pair_index +
                            repeat
                        )
                        %
                        2
                        ==
                        0
                    ) {
                        run_strict();
                        run_sensitive();

                    } else {

                        run_sensitive();
                        run_strict();
                    }
                }

                strict_pair_medians.push_back(
                    percentile(
                        strict_times,
                        0.50
                    )
                );

                sensitive_pair_medians.push_back(
                    percentile(
                        sensitive_times,
                        0.50
                    )
                );

                if (
                    contains_designed_target(
                        strict_result,
                        pair
                    )
                ) {
                    ++strict_targets;
                }

                if (
                    contains_designed_target(
                        sensitive_result,
                        pair
                    )
                ) {
                    ++sensitive_targets;
                }

                const auto strict_products =
                    canonical_products(
                        strict_result
                    );

                const auto sensitive_products =
                    canonical_products(
                        sensitive_result
                    );

                if (
                    !std::includes(
                        sensitive_products.begin(),
                        sensitive_products.end(),

                        strict_products.begin(),
                        strict_products.end(),

                        product_less
                    )
                ) {
                    ++subset_violations;
                }

                if (
                    budget == 0 &&
                    strict_products !=
                        sensitive_products
                ) {
                    ++equality_violations;
                }

                std::vector<
                    primerpair::PrimerPairHit
                > additional;

                std::set_difference(
                    sensitive_products.begin(),
                    sensitive_products.end(),

                    strict_products.begin(),
                    strict_products.end(),

                    std::back_inserter(
                        additional
                    ),

                    product_less
                );

                additional_total +=
                    additional.size();

                for (
                    const auto& hit :
                    additional
                ) {
                    const auto
                        total_mismatches =
                            hit.total_mismatches();

                    if (
                        total_mismatches <
                        mismatch_distribution.size()
                    ) {
                        ++mismatch_distribution.at(
                            total_mismatches
                        );
                    }

                    const bool anchor12 =
                        product_has_last_n_mismatch(
                            hit,
                            primer_length,
                            primer_length,
                            12
                        );

                    if (anchor12) {
                        ++add_anchor12;

                    } else {

                        ++anchor_explanation_violations;
                    }

                    if (
                        product_has_last_n_mismatch(
                            hit,
                            primer_length,
                            primer_length,
                            8
                        )
                    ) {
                        ++add_last8;
                    }

                    if (
                        product_has_last_n_mismatch(
                            hit,
                            primer_length,
                            primer_length,
                            5
                        )
                    ) {
                        ++add_last5;
                    }

                    if (
                        product_has_last_n_mismatch(
                            hit,
                            primer_length,
                            primer_length,
                            3
                        )
                    ) {
                        ++add_last3;
                    }

                    if (
                        product_has_last_n_mismatch(
                            hit,
                            primer_length,
                            primer_length,
                            1
                        )
                    ) {
                        ++add_terminal3;
                    }
                }

                const long double
                    strict_p1_hits =
                        static_cast<long double>(
                            strict_result
                                .primer1_single_hit_count
                        );

                const long double
                    strict_p2_hits =
                        static_cast<long double>(
                            strict_result
                                .primer2_single_hit_count
                        );

                const long double
                    sensitive_p1_hits =
                        static_cast<long double>(
                            sensitive_result
                                .primer1_single_hit_count
                        );

                const long double
                    sensitive_p2_hits =
                        static_cast<long double>(
                            sensitive_result
                                .primer2_single_hit_count
                        );

                strict_single_total +=
                    strict_p1_hits +
                    strict_p2_hits;

                sensitive_single_total +=
                    sensitive_p1_hits +
                    sensitive_p2_hits;

                strict_naive_total +=
                    strict_p1_hits *
                    strict_p2_hits;

                sensitive_naive_total +=
                    sensitive_p1_hits *
                    sensitive_p2_hits;

                strict_amplicon_total +=
                    strict_result
                        .amplicon_count();

                sensitive_amplicon_total +=
                    sensitive_result
                        .amplicon_count();
            }

            const long double
                single_hit_expansion =
                    strict_single_total > 0.0L
                        ? sensitive_single_total /
                            strict_single_total
                        : 0.0L;

            const long double
                pair_amplicon_expansion =
                    strict_amplicon_total > 0
                        ? static_cast<long double>(
                              sensitive_amplicon_total
                          )
                          /
                          static_cast<long double>(
                              strict_amplicon_total
                          )
                        : 0.0L;

            /*
             * >1 means pair constraint attenuated the
             * SENSITIVE expansion relative to the raw
             * single-primer hit expansion.
             */
            const long double
                expansion_compression =
                    pair_amplicon_expansion > 0.0L
                        ? single_hit_expansion /
                            pair_amplicon_expansion
                        : 0.0L;

            const long double
                strict_pair_reduction =
                    strict_amplicon_total > 0
                        ? strict_naive_total /
                            static_cast<long double>(
                                strict_amplicon_total
                            )
                        : 0.0L;

            const long double
                sensitive_pair_reduction =
                    sensitive_amplicon_total > 0
                        ? sensitive_naive_total /
                            static_cast<long double>(
                                sensitive_amplicon_total
                            )
                        : 0.0L;

            const double strict_mean =
                mean(
                    strict_pair_medians
                );

            const double sensitive_mean =
                mean(
                    sensitive_pair_medians
                );

            const double runtime_ratio =
                strict_mean > 0.0
                    ? sensitive_mean /
                        strict_mean
                    : 0.0;

            const bool all_checks =
                strict_targets ==
                    pairs.size()
                &&
                sensitive_targets ==
                    pairs.size()
                &&
                subset_violations == 0
                &&
                (
                    budget != 0 ||
                    equality_violations == 0
                )
                &&
                anchor_explanation_violations == 0;

            std::cout
                << budget
                << '\t'
                << static_cast<double>(
                       strict_single_total
                   )
                << '\t'
                << static_cast<double>(
                       sensitive_single_total
                   )
                << '\t'
                << static_cast<double>(
                       single_hit_expansion
                   )
                << '\t'
                << static_cast<double>(
                       strict_naive_total
                   )
                << '\t'
                << static_cast<double>(
                       sensitive_naive_total
                   )
                << '\t'
                << strict_amplicon_total
                << '\t'
                << sensitive_amplicon_total
                << '\t'
                << additional_total
                << '\t'
                << static_cast<double>(
                       pair_amplicon_expansion
                   )
                << '\t'
                << static_cast<double>(
                       expansion_compression
                   )
                << '\t'
                << static_cast<double>(
                       strict_pair_reduction
                   )
                << '\t'
                << static_cast<double>(
                       sensitive_pair_reduction
                   )
                << '\t'
                << add_anchor12
                << '\t'
                << add_last8
                << '\t'
                << add_last5
                << '\t'
                << add_last3
                << '\t'
                << add_terminal3
                << '\t'
                << mismatch_distribution.at(0)
                << '\t'
                << mismatch_distribution.at(1)
                << '\t'
                << mismatch_distribution.at(2)
                << '\t'
                << mismatch_distribution.at(3)
                << '\t'
                << mismatch_distribution.at(4)
                << '\t'
                << mismatch_distribution.at(5)
                << '\t'
                << mismatch_distribution.at(6)
                << '\t'
                << strict_targets
                << '/'
                << pairs.size()
                << '\t'
                << sensitive_targets
                << '/'
                << pairs.size()
                << '\t'
                << subset_violations
                << '\t'
                << equality_violations
                << '\t'
                << anchor_explanation_violations
                << '\t'
                << strict_mean
                << '\t'
                << sensitive_mean
                << '\t'
                << runtime_ratio
                << '\t'
                << percentile(
                       strict_pair_medians,
                       0.50
                   )
                << '\t'
                << percentile(
                       sensitive_pair_medians,
                       0.50
                   )
                << '\t'
                << percentile(
                       strict_pair_medians,
                       0.95
                   )
                << '\t'
                << percentile(
                       sensitive_pair_medians,
                       0.95
                   )
                << '\t'
                << (
                    all_checks
                        ? "YES"
                        : "NO"
                )
                << '\n';

            if (!all_checks) {
                throw std::logic_error(
                    "STRICT/SENSITIVE pair invariant failed."
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
