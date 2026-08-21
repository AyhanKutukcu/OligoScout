#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

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

std::size_t count_orientation(
    const primerpair::
        StrandAwarePrimerSearchResult& result,
    const primerpair::PrimerOrientation orientation
) {
    return static_cast<std::size_t>(
        std::count_if(
            result.hits.begin(),
            result.hits.end(),
            [
                orientation
            ](
                const auto& hit
            ) {
                return
                    hit.orientation ==
                    orientation;
            }
        )
    );
}

void normalize_products(
    std::vector<
        primerpair::PrimerPairHit
    >& products
) {
    std::sort(
        products.begin(),
        products.end(),
        [](
            const auto& lhs,
            const auto& rhs
        ) {
            if (
                lhs.amplicon_start !=
                rhs.amplicon_start
            ) {
                return
                    lhs.amplicon_start <
                    rhs.amplicon_start;
            }

            if (
                lhs.amplicon_end_exclusive !=
                rhs.amplicon_end_exclusive
            ) {
                return
                    lhs.amplicon_end_exclusive <
                    rhs.amplicon_end_exclusive;
            }

            if (
                lhs.left_primer !=
                rhs.left_primer
            ) {
                return
                    static_cast<int>(
                        lhs.left_primer
                    ) <
                    static_cast<int>(
                        rhs.left_primer
                    );
            }

            if (
                lhs.right_primer !=
                rhs.right_primer
            ) {
                return
                    static_cast<int>(
                        lhs.right_primer
                    ) <
                    static_cast<int>(
                        rhs.right_primer
                    );
            }

            if (
                lhs.left_mismatches !=
                rhs.left_mismatches
            ) {
                return
                    lhs.left_mismatches <
                    rhs.left_mismatches;
            }

            return
                lhs.right_mismatches <
                rhs.right_mismatches;
        }
    );

    products.erase(
        std::unique(
            products.begin(),
            products.end()
        ),
        products.end()
    );
}

void brute_orientation(
    const std::vector<
        primerpair::OrientedPrimerSearchHit
    >& forward_source,

    const std::size_t forward_length,
    const primerpair::PrimerIdentity forward_identity,

    const std::vector<
        primerpair::OrientedPrimerSearchHit
    >& reverse_source,

    const std::size_t reverse_length,
    const primerpair::PrimerIdentity reverse_identity,

    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon,

    std::vector<
        primerpair::PrimerPairHit
    >& output
) {
    for (const auto& forward : forward_source) {

        if (
            forward.orientation !=
            primerpair::PrimerOrientation::Forward
        ) {
            continue;
        }

        for (const auto& reverse : reverse_source) {

            if (
                reverse.orientation !=
                primerpair::PrimerOrientation::Reverse
            ) {
                continue;
            }

            /*
             * Non-overlapping inward-facing primers.
             */
            if (
                forward.position >
                std::numeric_limits<
                    std::uint64_t
                >::max() -
                    forward_length
            ) {
                continue;
            }

            const std::uint64_t forward_end =
                forward.position +
                static_cast<std::uint64_t>(
                    forward_length
                );

            if (
                reverse.position <
                forward_end
            ) {
                continue;
            }

            if (
                reverse.position >
                std::numeric_limits<
                    std::uint64_t
                >::max() -
                    reverse_length
            ) {
                continue;
            }

            const std::uint64_t amplicon_end =
                reverse.position +
                static_cast<std::uint64_t>(
                    reverse_length
                );

            if (
                amplicon_end <=
                forward.position
            ) {
                continue;
            }

            const std::uint64_t length =
                amplicon_end -
                forward.position;

            if (
                length < min_amplicon ||
                length > max_amplicon
            ) {
                continue;
            }

            output.push_back(
                primerpair::PrimerPairHit{
                    forward_identity,
                    reverse_identity,

                    forward.position,
                    reverse.position,

                    forward.mismatches,
                    reverse.mismatches,

                    forward.position,
                    amplicon_end,
                    length,

                    forward.mismatch_mask,
                    reverse.mismatch_mask
                }
            );
        }
    }
}

std::vector<primerpair::PrimerPairHit>
brute_force_products(
    const primerpair::
        StrandAwarePrimerSearchResult& primer1,

    const std::size_t primer1_length,

    const primerpair::
        StrandAwarePrimerSearchResult& primer2,

    const std::size_t primer2_length,

    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon
) {
    std::vector<
        primerpair::PrimerPairHit
    > products;

    brute_orientation(
        primer1.hits,
        primer1_length,
        primerpair::PrimerIdentity::Primer1,

        primer2.hits,
        primer2_length,
        primerpair::PrimerIdentity::Primer2,

        min_amplicon,
        max_amplicon,

        products
    );

    brute_orientation(
        primer2.hits,
        primer2_length,
        primerpair::PrimerIdentity::Primer2,

        primer1.hits,
        primer1_length,
        primerpair::PrimerIdentity::Primer1,

        min_amplicon,
        max_amplicon,

        products
    );

    normalize_products(
        products
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
                        Primer1 &&
                hit.right_primer ==
                    primerpair::
                        PrimerIdentity::
                        Primer2 &&
                hit.left_position ==
                    pair.left_position &&
                hit.right_position ==
                    pair.right_position &&
                hit.amplicon_length ==
                    pair.expected_amplicon_length;
        }
    );
}

double median(
    std::vector<double> values
) {
    if (values.empty()) {
        return 0.0;
    }

    std::sort(
        values.begin(),
        values.end()
    );

    const std::size_t middle =
        values.size() / 2;

    if (
        values.size() %
        2 == 1
    ) {
        return values.at(
            middle
        );
    }

    return
        (
            values.at(
                middle - 1
            ) +
            values.at(
                middle
            )
        ) /
        2.0;
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
                << " [pair_count]\n";

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
            stride = 104729ULL;

        const std::size_t target_pairs =
            argc == 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 256;

        if (target_pairs == 0) {
            throw std::invalid_argument(
                "Pair count must be > 0."
            );
        }

        const std::string reference =
            load_fasta(
                argv[1]
            );

        if (
            reference.size() <
            max_amplicon +
                primer_length
        ) {
            throw std::runtime_error(
                "Reference too short."
            );
        }

        std::vector<DesignedPair> pairs;

        pairs.reserve(
            target_pairs
        );

        std::uint64_t cursor = 0;

        const std::uint64_t usable =
            static_cast<std::uint64_t>(
                reference.size()
            ) -
            max_amplicon -
            primer_length;

        /*
         * Construct real-genome primer pairs from
         * chr22 itself.
         *
         * Primer1 = reference sequence on left.
         *
         * Primer2 = reverse-complement of reference
         * sequence at the right binding site.
         *
         * Therefore every designed pair has at least
         * one known valid genomic amplicon.
         */
        for (
            std::size_t attempt = 0;
            attempt <
                reference.size() &&
            pairs.size() <
                target_pairs;
            ++attempt
        ) {
            cursor =
                (
                    cursor +
                    stride
                ) %
                usable;

            /*
             * Deterministic amplicon length
             * distributed through 100..3000 bp.
             */
            const std::uint64_t
                amplicon_length =
                    100ULL +
                    (
                        (
                            static_cast<
                                std::uint64_t
                            >(
                                attempt
                            ) *
                            2654435761ULL
                        ) %
                        2901ULL
                    );

            const std::uint64_t
                right_position =
                    cursor +
                    amplicon_length -
                    primer_length;

            if (
                right_position +
                    primer_length >
                reference.size()
            ) {
                continue;
            }

            const std::string_view left_site(
                reference.data() +
                    cursor,
                primer_length
            );

            const std::string_view right_site(
                reference.data() +
                    right_position,
                primer_length
            );

            if (
                !all_acgt(
                    left_site
                ) ||
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

            /*
             * Avoid an uninformative identical
             * primer pair.
             */
            if (primer1 == primer2) {
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

        const primerpair::
            StrandAwarePrimerSearchEngine
                strand_engine(
                    index,
                    packed
                );

        const primerpair::
            PrimerPairSearchEngine
                pair_engine(
                    index,
                    packed
                );

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "reference_bp\t"
            << reference.size()
            << '\n';

        std::cout
            << "designed_pairs\t"
            << pairs.size()
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
            << "amplicon_range\t"
            << min_amplicon
            << '-'
            << max_amplicon
            << '\n';

        std::cout
            << "budget"
            << '\t'
            << "mean_primer1_hits"
            << '\t'
            << "mean_primer2_hits"
            << '\t'
            << "total_naive_cartesian"
            << '\t'
            << "total_inward_pair_space"
            << '\t'
            << "total_valid_amplicons"
            << '\t'
            << "orientation_reduction"
            << '\t'
            << "distance_reduction"
            << '\t'
            << "total_pair_reduction"
            << '\t'
            << "target_recovery"
            << '\t'
            << "bruteforce_verified"
            << '\t'
            << "pair_search_median_us"
            << '\n';

        /*
         * Exact and maximum-MVP mismatch budgets.
         *
         * k=0 demonstrates pair constraint alone.
         * k=3 demonstrates the stressed off-target
         * candidate space.
         */
        constexpr std::size_t budgets[] = {
            0,
            3
        };

        for (
            const std::size_t budget :
            budgets
        ) {
            long double
                primer1_hits_total = 0.0L;

            long double
                primer2_hits_total = 0.0L;

            long double
                naive_total = 0.0L;

            long double
                inward_total = 0.0L;

            long double
                valid_total = 0.0L;

            std::size_t targets_recovered = 0;

            std::size_t
                brute_verified = 0;

            std::vector<double>
                pair_times_us;

            pair_times_us.reserve(
                pairs.size()
            );

            for (const auto& pair : pairs) {

                /*
                 * Independent strand results provide
                 * candidate-space counts.
                 */
                const auto p1 =
                    strand_engine.search(
                        pair.primer1,
                        anchor_length,
                        budget
                    );

                const auto p2 =
                    strand_engine.search(
                        pair.primer2,
                        anchor_length,
                        budget
                    );

                const std::size_t p1_forward =
                    count_orientation(
                        p1,
                        primerpair::
                            PrimerOrientation::
                            Forward
                    );

                const std::size_t p1_reverse =
                    count_orientation(
                        p1,
                        primerpair::
                            PrimerOrientation::
                            Reverse
                    );

                const std::size_t p2_forward =
                    count_orientation(
                        p2,
                        primerpair::
                            PrimerOrientation::
                            Forward
                    );

                const std::size_t p2_reverse =
                    count_orientation(
                        p2,
                        primerpair::
                            PrimerOrientation::
                            Reverse
                    );

                const long double p1_hits =
                    static_cast<long double>(
                        p1.hit_count()
                    );

                const long double p2_hits =
                    static_cast<long double>(
                        p2.hit_count()
                    );

                const long double
                    naive_cartesian =
                        p1_hits *
                        p2_hits;

                const long double
                    inward_space =
                        static_cast<long double>(
                            p1_forward
                        ) *
                        static_cast<long double>(
                            p2_reverse
                        ) +
                        static_cast<long double>(
                            p2_forward
                        ) *
                        static_cast<long double>(
                            p1_reverse
                        );

                const auto start =
                    std::chrono::
                        steady_clock::now();

                const auto result =
                    pair_engine.search(
                        pair.primer1,
                        pair.primer2,
                        anchor_length,
                        budget,
                        min_amplicon,
                        max_amplicon
                    );

                const auto stop =
                    std::chrono::
                        steady_clock::now();

                const double elapsed_us =
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

                pair_times_us.push_back(
                    elapsed_us
                );

                if (
                    contains_designed_target(
                        result,
                        pair
                    )
                ) {
                    ++targets_recovered;
                }

                /*
                 * Independent brute-force
                 * differential verification.
                 *
                 * Restrict it to manageable Cartesian
                 * spaces so intentionally repeat-rich
                 * cases cannot explode benchmark time.
                 */
                if (
                    brute_verified < 64 &&
                    naive_cartesian <=
                        1000000.0L
                ) {
                    const auto brute =
                        brute_force_products(
                            p1,
                            pair.primer1.size(),
                            p2,
                            pair.primer2.size(),
                            min_amplicon,
                            max_amplicon
                        );

                    if (
                        brute !=
                        result.amplicons
                    ) {
                        throw std::logic_error(
                            "Optimized pair search "
                            "differs from brute force."
                        );
                    }

                    ++brute_verified;
                }

                primer1_hits_total +=
                    p1_hits;

                primer2_hits_total +=
                    p2_hits;

                naive_total +=
                    naive_cartesian;

                inward_total +=
                    inward_space;

                valid_total +=
                    static_cast<long double>(
                        result.amplicon_count()
                    );
            }

            const long double denominator =
                static_cast<long double>(
                    pairs.size()
                );

            const long double
                orientation_reduction =
                    inward_total > 0.0L
                        ? naive_total /
                            inward_total
                        : 0.0L;

            const long double
                distance_reduction =
                    valid_total > 0.0L
                        ? inward_total /
                            valid_total
                        : 0.0L;

            const long double
                total_reduction =
                    valid_total > 0.0L
                        ? naive_total /
                            valid_total
                        : 0.0L;

            std::cout
                << budget
                << '\t'
                << static_cast<double>(
                       primer1_hits_total /
                       denominator
                   )
                << '\t'
                << static_cast<double>(
                       primer2_hits_total /
                       denominator
                   )
                << '\t'
                << static_cast<double>(
                       naive_total
                   )
                << '\t'
                << static_cast<double>(
                       inward_total
                   )
                << '\t'
                << static_cast<double>(
                       valid_total
                   )
                << '\t'
                << static_cast<double>(
                       orientation_reduction
                   )
                << '\t'
                << static_cast<double>(
                       distance_reduction
                   )
                << '\t'
                << static_cast<double>(
                       total_reduction
                   )
                << '\t'
                << targets_recovered
                << '/'
                << pairs.size()
                << '\t'
                << brute_verified
                << '\t'
                << median(
                       pair_times_us
                   )
                << '\n';

            if (
                targets_recovered !=
                pairs.size()
            ) {
                throw std::logic_error(
                    "Designed target recovery "
                    "was not 100%."
                );
            }

            if (brute_verified == 0) {
                throw std::logic_error(
                    "No brute-force differential "
                    "cases were verified."
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
