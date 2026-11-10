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
#include <vector>

#include "primerpair/approximate_anchor_search.hpp"

namespace {

struct PrimerCase {
    std::string primer;
    std::uint64_t anchor_count{0};
    std::uint64_t full_count{0};
};

enum class Difficulty {
    Easy,
    Moderate,
    Hard,
    RepeatRich
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
                        "Unsupported FASTA symbol."
                    );
            }
        }
    }

    if (sequence.empty()) {
        throw std::runtime_error(
            "Empty FASTA sequence."
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
        [](const char base) {
            return
                base == 'A' ||
                base == 'C' ||
                base == 'G' ||
                base == 'T';
        }
    );
}

Difficulty classify(
    const std::uint64_t anchor_count
) {
    if (anchor_count <= 10) {
        return Difficulty::Easy;
    }

    if (anchor_count <= 100) {
        return Difficulty::Moderate;
    }

    if (anchor_count <= 1000) {
        return Difficulty::Hard;
    }

    return Difficulty::RepeatRich;
}

const char* difficulty_name(
    const Difficulty difficulty
) {
    switch (difficulty) {
        case Difficulty::Easy:
            return "EASY";

        case Difficulty::Moderate:
            return "MODERATE";

        case Difficulty::Hard:
            return "HARD";

        case Difficulty::RepeatRich:
            return "REPEAT_RICH";
    }

    return "UNKNOWN";
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
                << " [queries_per_bin_budget]\n";

            return 2;
        }

        constexpr std::size_t
            primer_length = 20;

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            target_pool = 4096;

        constexpr std::uint64_t
            stride = 104729ULL;

        const std::size_t queries =
            argc == 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 5000;

        if (queries == 0) {
            throw std::invalid_argument(
                "Query count must be > 0."
            );
        }

        std::string reference =
            load_fasta(argv[1]);

        const std::size_t available =
            reference.size() -
            primer_length +
            1;

        std::vector<std::string> primers;

        primers.reserve(
            target_pool
        );

        std::uint64_t candidate = 0;

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

            primers.emplace_back(window);
        }

        if (primers.size() != target_pool) {
            throw std::runtime_error(
                "Could not collect 4096 primers."
            );
        }

        primerpair::BidirectionalFMIndex index(
            std::move(reference)
        );

        const primerpair::ApproximateAnchorSearcher
            searcher(index);

        std::array<
            std::vector<PrimerCase>,
            4
        > bins;

        for (const auto& primer : primers) {

            const std::string_view anchor(
                primer.data() +
                    primer.size() -
                    anchor_length,
                anchor_length
            );

            const std::uint64_t anchor_count =
                index.search(anchor).size();

            const std::uint64_t full_count =
                index.search(primer).size();

            PrimerCase current{
                primer,
                anchor_count,
                full_count
            };

            switch (classify(anchor_count)) {
                case Difficulty::Easy:
                    bins.at(0).push_back(
                        std::move(current)
                    );
                    break;

                case Difficulty::Moderate:
                    bins.at(1).push_back(
                        std::move(current)
                    );
                    break;

                case Difficulty::Hard:
                    bins.at(2).push_back(
                        std::move(current)
                    );
                    break;

                case Difficulty::RepeatRich:
                    bins.at(3).push_back(
                        std::move(current)
                    );
                    break;
            }
        }

        constexpr std::array<Difficulty, 4>
            difficulties{
                Difficulty::Easy,
                Difficulty::Moderate,
                Difficulty::Hard,
                Difficulty::RepeatRich
            };

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "primer_pool\t"
            << primers.size()
            << '\n';

        std::cout
            << "queries_per_bin_budget\t"
            << queries
            << '\n';

        std::cout
            << "difficulty"
            << '\t'
            << "bin_size"
            << '\t'
            << "mean_anchor_count"
            << '\t'
            << "mean_full20_count"
            << '\t'
            << "budget"
            << '\t'
            << "ns_per_query"
            << '\t'
            << "qps"
            << '\t'
            << "mean_final_branches"
            << '\t'
            << "mean_matches"
            << '\n';

        std::uint64_t checksum = 0;

        for (std::size_t b = 0;
             b < bins.size();
             ++b) {

            const auto& cases =
                bins.at(b);

            if (cases.empty()) {
                continue;
            }

            long double anchor_sum = 0.0;
            long double full_sum = 0.0;

            for (const auto& item : cases) {
                anchor_sum +=
                    item.anchor_count;

                full_sum +=
                    item.full_count;
            }

            const double mean_anchor =
                static_cast<double>(
                    anchor_sum /
                    cases.size()
                );

            const double mean_full =
                static_cast<double>(
                    full_sum /
                    cases.size()
                );

            /*
             * Warm-up.
             */
            for (std::size_t i = 0;
                 i < 100;
                 ++i) {

                const auto result =
                    searcher
                        .search_5prime_mismatches(
                            cases.at(
                                i %
                                cases.size()
                            ).primer,
                            anchor_length,
                            3
                        );

                checksum +=
                    result.hits.size();
            }

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
                     i < queries;
                     ++i) {

                    const auto& item =
                        cases.at(
                            i %
                            cases.size()
                        );

                    const auto result =
                        searcher
                            .search_5prime_mismatches(
                                item.primer,
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

                const double denominator =
                    static_cast<double>(
                        queries
                    );

                const double ns_per_query =
                    elapsed_ns /
                    denominator;

                std::cout
                    << difficulty_name(
                           difficulties.at(b)
                       )
                    << '\t'
                    << cases.size()
                    << '\t'
                    << mean_anchor
                    << '\t'
                    << mean_full
                    << '\t'
                    << budget
                    << '\t'
                    << ns_per_query
                    << '\t'
                    << 1e9 / ns_per_query
                    << '\t'
                    << static_cast<double>(
                           total_branches
                       ) /
                       denominator
                    << '\t'
                    << static_cast<double>(
                           total_matches
                       ) /
                       denominator
                    << '\n';
            }
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
