#include <algorithm>
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
#include <unordered_set>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/sensitive_candidate_search.hpp"
#include "primerpair/sensitive_primer_search.hpp"

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

    return values.at(
        static_cast<std::size_t>(
            scaled
        )
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
                << " <fasta> [primer_count]\n";

            return 2;
        }

        constexpr std::size_t
            primer_length = 20;

        constexpr std::uint64_t
            stride = 104729ULL;

        const std::size_t primer_count =
            argc == 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 64;

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

        const std::uint64_t available =
            static_cast<std::uint64_t>(
                reference.size() -
                primer_length +
                1
            );

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
                !seen.insert(primer).second
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
                "Could not collect requested primers."
            );
        }

        const primerpair::PackedReference
            packed(reference);

        const primerpair::BidirectionalFMIndex
            index(reference);

        const primerpair::
            SensitivePrimerSearchEngine
                reference_engine(index);

        const primerpair::
            SensitiveCandidateSearchEngine
                candidate_engine(
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
            << "primer_count\t"
            << primers.size()
            << '\n';

        std::cout
            << "primer_length\t"
            << primer_length
            << '\n';

        std::cout
            << "budget"
            << '\t'
            << "reference_hits"
            << '\t'
            << "candidate_hits"
            << '\t'
            << "equality_violations"
            << '\t'
            << "reference_median_us"
            << '\t'
            << "candidate_median_us"
            << '\t'
            << "speedup"
            << '\t'
            << "reference_p95_us"
            << '\t'
            << "candidate_p95_us"
            << '\t'
            << "all_checks"
            << '\n';

        for (
            std::size_t budget = 0;
            budget <= 3;
            ++budget
        ) {
            std::uint64_t
                reference_hits_total = 0;

            std::uint64_t
                candidate_hits_total = 0;

            std::uint64_t
                equality_violations = 0;

            std::vector<double>
                reference_times;

            std::vector<double>
                candidate_times;

            reference_times.reserve(
                primers.size()
            );

            candidate_times.reserve(
                primers.size()
            );

            for (
                std::size_t i = 0;
                i < primers.size();
                ++i
            ) {
                const auto& primer =
                    primers.at(i);

                primerpair::
                    SensitivePrimerSearchResult
                        reference_result;

                primerpair::
                    SensitivePrimerSearchResult
                        candidate_result;

                double reference_us = 0.0;
                double candidate_us = 0.0;

                /*
                 * Alternate timing order to reduce
                 * systematic cache/order bias.
                 */
                if (
                    (i + budget) % 2 == 0
                ) {
                    {
                        const auto start =
                            std::chrono::
                                steady_clock::now();

                        reference_result =
                            reference_engine.search(
                                primer,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        reference_us =
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

                        candidate_result =
                            candidate_engine.search(
                                primer,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        candidate_us =
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

                        candidate_result =
                            candidate_engine.search(
                                primer,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        candidate_us =
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

                        reference_result =
                            reference_engine.search(
                                primer,
                                budget
                            );

                        const auto stop =
                            std::chrono::
                                steady_clock::now();

                        reference_us =
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

                reference_times.push_back(
                    reference_us
                );

                candidate_times.push_back(
                    candidate_us
                );

                reference_hits_total +=
                    reference_result.hit_count();

                candidate_hits_total +=
                    candidate_result.hit_count();

                if (
                    candidate_result.hits !=
                    reference_result.hits
                ) {
                    ++equality_violations;
                }
            }

            const double reference_median =
                percentile(
                    reference_times,
                    0.50
                );

            const double candidate_median =
                percentile(
                    candidate_times,
                    0.50
                );

            const double reference_p95 =
                percentile(
                    reference_times,
                    0.95
                );

            const double candidate_p95 =
                percentile(
                    candidate_times,
                    0.95
                );

            const double speedup =
                candidate_median > 0.0
                    ? reference_median /
                      candidate_median
                    : 0.0;

            const bool all_checks =
                equality_violations == 0 &&
                reference_hits_total ==
                    candidate_hits_total;

            std::cout
                << budget
                << '\t'
                << reference_hits_total
                << '\t'
                << candidate_hits_total
                << '\t'
                << equality_violations
                << '\t'
                << reference_median
                << '\t'
                << candidate_median
                << '\t'
                << speedup
                << '\t'
                << reference_p95
                << '\t'
                << candidate_p95
                << '\t'
                << (
                    all_checks
                        ? "YES"
                        : "NO"
                )
                << '\n';

            if (!all_checks) {
                throw std::logic_error(
                    "Candidate search differs from "
                    "exhaustive reference."
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
