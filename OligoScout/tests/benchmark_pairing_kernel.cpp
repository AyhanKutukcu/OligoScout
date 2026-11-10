#include <algorithm>
#include <array>
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

struct PairCase {
    std::vector<primerpair::OrientedPrimerSearchHit> p1_forward;
    std::vector<primerpair::OrientedPrimerSearchHit> p1_reverse;

    std::vector<primerpair::OrientedPrimerSearchHit> p2_forward;
    std::vector<primerpair::OrientedPrimerSearchHit> p2_reverse;
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

    return sequence;
}

bool all_acgt(
    const std::string_view sequence
) {
    return std::all_of(
        sequence.begin(),
        sequence.end(),
        [](const char base) {
            return
                base == 'A' ||
                base == 'C' ||
                base == 'G' ||
                base == 'T';
        }
    );
}

std::vector<primerpair::OrientedPrimerSearchHit>
extract(
    const primerpair::StrandAwarePrimerSearchResult& result,
    const primerpair::PrimerOrientation orientation
) {
    std::vector<
        primerpair::OrientedPrimerSearchHit
    > output;

    for (const auto& hit : result.hits) {

        if (hit.orientation == orientation) {
            output.push_back(hit);
        }
    }

    return output;
}

std::uint64_t saturating_add(
    const std::uint64_t lhs,
    const std::uint64_t rhs
) noexcept {
    if (
        lhs >
        std::numeric_limits<
            std::uint64_t
        >::max() - rhs
    ) {
        return
            std::numeric_limits<
                std::uint64_t
            >::max();
    }

    return lhs + rhs;
}

void normalize(
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

void emit(
    const primerpair::OrientedPrimerSearchHit& forward,
    const primerpair::PrimerIdentity forward_identity,

    const primerpair::OrientedPrimerSearchHit& reverse,
    const primerpair::PrimerIdentity reverse_identity,

    const std::uint64_t reverse_length,

    std::vector<
        primerpair::PrimerPairHit
    >& output
) {
    const std::uint64_t end =
        reverse.position +
        reverse_length;

    output.push_back(
        primerpair::PrimerPairHit{
            forward_identity,
            reverse_identity,

            forward.position,
            reverse.position,

            forward.mismatches,
            reverse.mismatches,

            forward.position,
            end,
            end - forward.position
        }
    );
}

/*
 * --------------------------------------------------
 * OLD BASELINE:
 *
 * O(F log R + K)
 * --------------------------------------------------
 */
void append_lower_bound(
    const std::vector<
        primerpair::OrientedPrimerSearchHit
    >& forward_hits,

    const primerpair::PrimerIdentity forward_identity,

    const std::vector<
        primerpair::OrientedPrimerSearchHit
    >& reverse_hits,

    const primerpair::PrimerIdentity reverse_identity,

    const std::uint64_t primer_length,
    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon,

    std::vector<
        primerpair::PrimerPairHit
    >& output
) {
    for (const auto& forward : forward_hits) {

        std::uint64_t minimum =
            saturating_add(
                forward.position,
                primer_length
            );

        if (
            min_amplicon >
            primer_length
        ) {
            minimum =
                std::max(
                    minimum,
                    saturating_add(
                        forward.position,
                        min_amplicon -
                            primer_length
                    )
                );
        }

        const std::uint64_t maximum =
            saturating_add(
                forward.position,
                max_amplicon -
                    primer_length
            );

        auto it =
            std::lower_bound(
                reverse_hits.begin(),
                reverse_hits.end(),
                minimum,
                [](
                    const auto& hit,
                    const std::uint64_t value
                ) {
                    return
                        hit.position <
                        value;
                }
            );

        for (
            ;
            it != reverse_hits.end() &&
            it->position <= maximum;
            ++it
        ) {
            emit(
                forward,
                forward_identity,

                *it,
                reverse_identity,

                primer_length,
                output
            );
        }
    }
}

/*
 * --------------------------------------------------
 * NEW SWEEP-LINE:
 *
 * O(F + R + K)
 * --------------------------------------------------
 */
void append_sweep(
    const std::vector<
        primerpair::OrientedPrimerSearchHit
    >& forward_hits,

    const primerpair::PrimerIdentity forward_identity,

    const std::vector<
        primerpair::OrientedPrimerSearchHit
    >& reverse_hits,

    const primerpair::PrimerIdentity reverse_identity,

    const std::uint64_t primer_length,
    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon,

    std::vector<
        primerpair::PrimerPairHit
    >& output
) {
    std::size_t lower = 0;
    std::size_t upper = 0;

    for (const auto& forward : forward_hits) {

        std::uint64_t minimum =
            saturating_add(
                forward.position,
                primer_length
            );

        if (
            min_amplicon >
            primer_length
        ) {
            minimum =
                std::max(
                    minimum,
                    saturating_add(
                        forward.position,
                        min_amplicon -
                            primer_length
                    )
                );
        }

        const std::uint64_t maximum =
            saturating_add(
                forward.position,
                max_amplicon -
                    primer_length
            );

        while (
            lower <
                reverse_hits.size() &&
            reverse_hits.at(lower)
                .position <
                minimum
        ) {
            ++lower;
        }

        if (upper < lower) {
            upper = lower;
        }

        while (
            upper <
                reverse_hits.size() &&
            reverse_hits.at(upper)
                .position <=
                maximum
        ) {
            ++upper;
        }

        for (
            std::size_t i = lower;
            i < upper;
            ++i
        ) {
            emit(
                forward,
                forward_identity,

                reverse_hits.at(i),
                reverse_identity,

                primer_length,
                output
            );
        }
    }
}

std::vector<primerpair::PrimerPairHit>
pair_lower_bound(
    const PairCase& data,
    const std::uint64_t primer_length,
    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon
) {
    std::vector<
        primerpair::PrimerPairHit
    > output;

    append_lower_bound(
        data.p1_forward,
        primerpair::PrimerIdentity::Primer1,

        data.p2_reverse,
        primerpair::PrimerIdentity::Primer2,

        primer_length,
        min_amplicon,
        max_amplicon,

        output
    );

    append_lower_bound(
        data.p2_forward,
        primerpair::PrimerIdentity::Primer2,

        data.p1_reverse,
        primerpair::PrimerIdentity::Primer1,

        primer_length,
        min_amplicon,
        max_amplicon,

        output
    );

    normalize(output);

    return output;
}

std::vector<primerpair::PrimerPairHit>
pair_sweep(
    const PairCase& data,
    const std::uint64_t primer_length,
    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon
) {
    std::vector<
        primerpair::PrimerPairHit
    > output;

    append_sweep(
        data.p1_forward,
        primerpair::PrimerIdentity::Primer1,

        data.p2_reverse,
        primerpair::PrimerIdentity::Primer2,

        primer_length,
        min_amplicon,
        max_amplicon,

        output
    );

    append_sweep(
        data.p2_forward,
        primerpair::PrimerIdentity::Primer2,

        data.p1_reverse,
        primerpair::PrimerIdentity::Primer1,

        primer_length,
        min_amplicon,
        max_amplicon,

        output
    );

    normalize(output);

    return output;
}

double median(
    std::vector<double> values
) {
    std::sort(
        values.begin(),
        values.end()
    );

    const std::size_t middle =
        values.size() / 2;

    if (values.size() % 2 == 1) {
        return values.at(middle);
    }

    return
        (
            values.at(middle - 1) +
            values.at(middle)
        ) /
        2.0;
}

std::uint64_t checksum(
    const std::vector<
        primerpair::PrimerPairHit
    >& products
) {
    std::uint64_t value =
        static_cast<std::uint64_t>(
            products.size()
        );

    for (const auto& product : products) {

        value ^=
            product.amplicon_start +
            0x9E3779B97F4A7C15ULL;

        value ^=
            product.amplicon_end_exclusive;
    }

    return value;
}

}  // namespace

int main(
    int argc,
    char* argv[]
) {
    try {
        if (argc < 2 || argc > 4) {

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
            stride = 104729ULL;

        const std::size_t pair_count =
            argc >= 3
                ? std::stoull(argv[2])
                : 256;

        const std::size_t repeats =
            argc >= 4
                ? std::stoull(argv[3])
                : 9;

        const std::string reference =
            load_fasta(
                argv[1]
            );

        const primerpair::PackedReference
            packed(reference);

        const primerpair::BidirectionalFMIndex
            index(reference);

        const primerpair::
            StrandAwarePrimerSearchEngine
                searcher(
                    index,
                    packed
                );

        constexpr std::array<
            std::size_t,
            2
        > budgets{
            0,
            3
        };

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "reference_bp\t"
            << reference.size()
            << '\n';

        std::cout
            << "pair_count\t"
            << pair_count
            << '\n';

        std::cout
            << "repeats\t"
            << repeats
            << '\n';

        std::cout
            << "budget"
            << '\t'
            << "mean_total_input_hits"
            << '\t'
            << "total_products"
            << '\t'
            << "lower_bound_median_us"
            << '\t'
            << "sweep_median_us"
            << '\t'
            << "speedup"
            << '\t'
            << "sweep_wins"
            << '\t'
            << "equivalent"
            << '\n';

        for (const auto budget : budgets) {

            std::vector<PairCase> cases;

            cases.reserve(
                pair_count
            );

            std::uint64_t cursor = 0;

            const std::uint64_t usable =
                reference.size() -
                max_amplicon -
                primer_length;

            for (
                std::size_t attempt = 0;
                attempt < reference.size() &&
                cases.size() < pair_count;
                ++attempt
            ) {
                cursor =
                    (
                        cursor +
                        stride
                    ) %
                    usable;

                const std::uint64_t length =
                    100ULL +
                    (
                        (
                            static_cast<
                                std::uint64_t
                            >(attempt) *
                            2654435761ULL
                        ) %
                        2901ULL
                    );

                const std::uint64_t right =
                    cursor +
                    length -
                    primer_length;

                const std::string_view left_site(
                    reference.data() + cursor,
                    primer_length
                );

                const std::string_view right_site(
                    reference.data() + right,
                    primer_length
                );

                if (
                    !all_acgt(left_site) ||
                    !all_acgt(right_site)
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

                const auto p1 =
                    searcher.search(
                        primer1,
                        anchor_length,
                        budget
                    );

                const auto p2 =
                    searcher.search(
                        primer2,
                        anchor_length,
                        budget
                    );

                cases.push_back(
                    PairCase{
                        extract(
                            p1,
                            primerpair::
                                PrimerOrientation::
                                Forward
                        ),

                        extract(
                            p1,
                            primerpair::
                                PrimerOrientation::
                                Reverse
                        ),

                        extract(
                            p2,
                            primerpair::
                                PrimerOrientation::
                                Forward
                        ),

                        extract(
                            p2,
                            primerpair::
                                PrimerOrientation::
                                Reverse
                        )
                    }
                );
            }

            if (
                cases.size() !=
                pair_count
            ) {
                throw std::runtime_error(
                    "Could not collect pair cases."
                );
            }

            long double input_hits = 0.0L;
            std::uint64_t total_products = 0;

            /*
             * Differential correctness before timing.
             */
            for (const auto& data : cases) {

                const auto old_result =
                    pair_lower_bound(
                        data,
                        primer_length,
                        min_amplicon,
                        max_amplicon
                    );

                const auto new_result =
                    pair_sweep(
                        data,
                        primer_length,
                        min_amplicon,
                        max_amplicon
                    );

                if (old_result != new_result) {
                    throw std::logic_error(
                        "Lower-bound and sweep-line "
                        "pairing differ."
                    );
                }

                input_hits +=
                    data.p1_forward.size() +
                    data.p1_reverse.size() +
                    data.p2_forward.size() +
                    data.p2_reverse.size();

                total_products +=
                    new_result.size();
            }

            std::vector<double>
                lower_times;

            std::vector<double>
                sweep_times;

            std::size_t sweep_wins = 0;

            std::uint64_t lower_checksum_guard = 0;
            std::uint64_t sweep_checksum_guard = 0;

            for (
                std::size_t repeat = 0;
                repeat < repeats;
                ++repeat
            ) {
                auto run_lower =
                    [&]() {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        std::uint64_t sum = 0;

                        for (const auto& data : cases) {

                            sum ^=
                                checksum(
                                    pair_lower_bound(
                                        data,
                                        primer_length,
                                        min_amplicon,
                                        max_amplicon
                                    )
                                );
                        }

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        lower_checksum_guard = sum;

                        return
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0 /
                            cases.size();
                    };

                auto run_sweep =
                    [&]() {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        std::uint64_t sum = 0;

                        for (const auto& data : cases) {

                            sum ^=
                                checksum(
                                    pair_sweep(
                                        data,
                                        primer_length,
                                        min_amplicon,
                                        max_amplicon
                                    )
                                );
                        }

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        sweep_checksum_guard = sum;

                        return
                            static_cast<double>(
                                std::chrono::
                                    duration_cast<
                                        std::chrono::
                                            nanoseconds
                                    >(
                                        stop - start
                                    ).count()
                            ) /
                            1000.0 /
                            cases.size();
                    };

                double lower_us = 0.0;
                double sweep_us = 0.0;

                /*
                 * Alternate timing order.
                 */
                if (repeat % 2 == 0) {

                    lower_us =
                        run_lower();

                    sweep_us =
                        run_sweep();

                } else {

                    sweep_us =
                        run_sweep();

                    lower_us =
                        run_lower();
                }

                lower_times.push_back(
                    lower_us
                );

                sweep_times.push_back(
                    sweep_us
                );

                if (sweep_us < lower_us) {
                    ++sweep_wins;
                }
            }

            const double lower_median =
                median(
                    lower_times
                );

            const double sweep_median =
                median(
                    sweep_times
                );

            std::cout
                << budget
                << '\t'
                << static_cast<double>(
                       input_hits /
                       cases.size()
                   )
                << '\t'
                << total_products
                << '\t'
                << lower_median
                << '\t'
                << sweep_median
                << '\t'
                << (
                    lower_median /
                    sweep_median
                )
                << '\t'
                << sweep_wins
                << '/'
                << repeats
                << '\t'
                << "YES"
                << '\n';

            if (
                lower_checksum_guard !=
                sweep_checksum_guard
            ) {
                throw std::logic_error(
                    "Timed lower-bound and sweep-line "
                    "checksums differ."
                );
            }

            std::cerr
                << "checksum_lower("
                << budget
                << ")\t"
                << lower_checksum_guard
                << '\n';

            std::cerr
                << "checksum_sweep("
                << budget
                << ")\t"
                << sweep_checksum_guard
                << '\n';
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
