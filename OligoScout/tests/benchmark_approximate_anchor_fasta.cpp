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
#include <utility>
#include <vector>

#include "primerpair/approximate_anchor_search.hpp"

namespace {

std::string load_single_fasta(
    const std::string& path
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open FASTA: " + path
        );
    }

    std::string sequence;
    std::string line;

    std::uint64_t ambiguous_to_n = 0;

    while (std::getline(input, line)) {

        if (line.empty()) {
            continue;
        }

        if (line.front() == '>') {
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

                /*
                 * Standard IUPAC ambiguous DNA symbols.
                 * Exact referans bazı olarak kullanmıyoruz;
                 * bunları N olarak indeksliyoruz.
                 */
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
                    ++ambiguous_to_n;
                    break;

                default:
                    throw std::runtime_error(
                        std::string(
                            "Unsupported FASTA nucleotide: "
                        ) +
                        base
                    );
            }
        }
    }

    if (sequence.empty()) {
        throw std::runtime_error(
            "FASTA contains no sequence."
        );
    }

    std::cerr
        << "FASTA ambiguous IUPAC -> N: "
        << ambiguous_to_n
        << '\n';

    return sequence;
}

bool all_acgt(
    const std::string_view sequence
) {
    for (const char base : sequence) {

        if (
            base != 'A' &&
            base != 'C' &&
            base != 'G' &&
            base != 'T'
        ) {
            return false;
        }
    }

    return true;
}

}  // namespace

int main(
    int argc,
    char* argv[]
) {
    try {
        if (argc < 2 || argc > 5) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " [iterations]"
                << " [primer_length]"
                << " [anchor_length]\n";

            return 2;
        }

        const std::string fasta_path =
            argv[1];

        const std::size_t iterations =
            argc >= 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 5000;

        const std::size_t primer_length =
            argc >= 4
                ? static_cast<std::size_t>(
                      std::stoull(argv[3])
                  )
                : 20;

        const std::size_t anchor_length =
            argc >= 5
                ? static_cast<std::size_t>(
                      std::stoull(argv[4])
                  )
                : 12;

        if (iterations == 0) {
            throw std::invalid_argument(
                "Iterations must be greater than zero."
            );
        }

        if (primer_length == 0) {
            throw std::invalid_argument(
                "Primer length must be greater than zero."
            );
        }

        if (
            anchor_length == 0 ||
            anchor_length > primer_length
        ) {
            throw std::invalid_argument(
                "Invalid anchor length."
            );
        }

        std::string reference =
            load_single_fasta(
                fasta_path
            );

        if (reference.size() < primer_length) {
            throw std::runtime_error(
                "Reference shorter than primer."
            );
        }

        constexpr std::size_t target_pool =
            512;

        constexpr std::uint64_t stride =
            104729ULL;

        const std::size_t available =
            reference.size() -
            primer_length +
            1;

        std::vector<std::string> primers;

        primers.reserve(
            target_pool
        );

        std::uint64_t candidate = 0;

        /*
         * Yalnız tamamen A/C/G/T olan primer
         * pencerelerini benchmark havuzuna al.
         *
         * N veya ambiguity içeren bölgeler dışlanır.
         */
        for (
            std::size_t attempts = 0;
            attempts < available &&
            primers.size() < target_pool;
            ++attempts
        ) {
            candidate =
                (
                    candidate +
                    stride
                ) %
                available;

            const std::string_view window(
                reference.data() + candidate,
                primer_length
            );

            if (!all_acgt(window)) {
                continue;
            }

            primers.emplace_back(
                window
            );
        }

        if (primers.size() < 100) {
            throw std::runtime_error(
                "Unable to collect enough ACGT primers."
            );
        }

        primerpair::BidirectionalFMIndex index(
            std::move(reference)
        );

        const primerpair::ApproximateAnchorSearcher
            searcher(index);

        /*
         * --------------------------------------------------
         * 12-mer anchor tekrar sıklığı
         * --------------------------------------------------
         */

        double anchor_occurrence_sum =
            0.0;

        std::uint64_t anchor_occurrence_max =
            0;

        std::size_t repeated_anchor_count =
            0;

        for (const auto& primer : primers) {

            const std::string_view anchor(
                primer.data() +
                    primer.size() -
                    anchor_length,
                anchor_length
            );

            const auto state =
                index.search(anchor);

            anchor_occurrence_sum +=
                static_cast<double>(
                    state.size()
                );

            anchor_occurrence_max =
                std::max(
                    anchor_occurrence_max,
                    state.size()
                );

            if (state.size() > 1) {
                ++repeated_anchor_count;
            }
        }

        const double mean_anchor_occurrences =
            anchor_occurrence_sum /
            static_cast<double>(
                primers.size()
            );

        const double repeated_anchor_fraction =
            static_cast<double>(
                repeated_anchor_count
            ) /
            static_cast<double>(
                primers.size()
            );

        std::uint64_t checksum = 0;

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "reference_length_bp\t"
            << index.forward_index()
                    .indexed_text_length()
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
            << "pattern_pool_size\t"
            << primers.size()
            << '\n';

        std::cout
            << "mean_anchor_occurrences\t"
            << mean_anchor_occurrences
            << '\n';

        std::cout
            << "max_anchor_occurrences\t"
            << anchor_occurrence_max
            << '\n';

        std::cout
            << "repeated_anchor_fraction\t"
            << repeated_anchor_fraction
            << '\n';

        std::cout
            << "mismatch_budget"
            << '\t'
            << "ns_per_query"
            << '\t'
            << "queries_per_second"
            << '\t'
            << "mean_final_branches"
            << '\t'
            << "mean_match_count"
            << '\n';

        for (std::size_t budget = 0;
             budget <= 3;
             ++budget) {

            std::uint64_t total_branches =
                0;

            std::uint64_t total_matches =
                0;

            const auto start =
                std::chrono::steady_clock::now();

            for (std::size_t i = 0;
                 i < iterations;
                 ++i) {

                const auto result =
                    searcher.search_5prime_mismatches(
                        primers.at(
                            i % primers.size()
                        ),
                        anchor_length,
                        budget
                    );

                total_branches +=
                    static_cast<std::uint64_t>(
                        result.hits.size()
                    );

                total_matches +=
                    result.total_match_count();

                checksum +=
                    static_cast<std::uint64_t>(
                        result.hits.size()
                    );

                checksum +=
                    result.total_match_count();
            }

            const auto stop =
                std::chrono::steady_clock::now();

            const double elapsed_ns =
                static_cast<double>(
                    std::chrono::duration_cast<
                        std::chrono::nanoseconds
                    >(
                        stop - start
                    ).count()
                );

            const double count =
                static_cast<double>(
                    iterations
                );

            const double ns_per_query =
                elapsed_ns /
                count;

            std::cout
                << budget
                << '\t'
                << ns_per_query
                << '\t'
                << 1e9 / ns_per_query
                << '\t'
                << static_cast<double>(
                       total_branches
                   ) / count
                << '\t'
                << static_cast<double>(
                       total_matches
                   ) / count
                << '\n';
        }

        std::cout
            << "checksum\t"
            << checksum
            << '\n';

        return 0;

    } catch (const std::exception& error) {
        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
