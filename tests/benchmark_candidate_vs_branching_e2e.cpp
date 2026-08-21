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
#include <utility>
#include <vector>

#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/approximate_anchor_search.hpp"
#include "primerpair/packed_reference.hpp"

namespace {

using PositionMismatch =
    std::pair<
        std::uint64_t,
        std::size_t
    >;

struct PrimerCase {
    std::string primer;
    std::uint64_t anchor_count{0};
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
    const std::uint64_t count
) {
    if (count <= 10) {
        return Difficulty::Easy;
    }

    if (count <= 100) {
        return Difficulty::Moderate;
    }

    if (count <= 1000) {
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

std::vector<PositionMismatch>
materialize_branching(
    const primerpair::BidirectionalFMIndex& index,
    const primerpair::ApproximateAnchorSearchResult& result
) {
    std::vector<PositionMismatch> output;

    /*
     * total_match_count is a useful lower-bound
     * reservation estimate here.
     */
    output.reserve(
        static_cast<std::size_t>(
            result.total_match_count()
        )
    );

    for (const auto& hit : result.hits) {

        const auto positions =
            index.locate(
                hit.state
            );

        for (const auto position : positions) {

            output.emplace_back(
                position,
                hit.mismatches
            );
        }
    }

    /*
     * A usable genomic hit list should be
     * coordinate ordered and deduplicated.
     */
    std::sort(
        output.begin(),
        output.end()
    );

    output.erase(
        std::unique(
            output.begin(),
            output.end()
        ),
        output.end()
    );

    return output;
}

std::vector<PositionMismatch>
materialize_candidate(
    const primerpair::AnchorCandidateSearchResult& result
) {
    std::vector<PositionMismatch> output;

    output.reserve(
        result.hits.size()
    );

    for (const auto& hit : result.hits) {

        output.emplace_back(
            hit.position,
            hit.mismatches
        );
    }

    /*
     * AnchorCandidateSearcher should already
     * naturally produce coordinate-ordered hits,
     * but normalize here for differential checking.
     */
    std::sort(
        output.begin(),
        output.end()
    );

    output.erase(
        std::unique(
            output.begin(),
            output.end()
        ),
        output.end()
    );

    return output;
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
                : 500;

        if (queries == 0) {
            throw std::invalid_argument(
                "Query count must be > 0."
            );
        }

        std::string reference =
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

            primers.emplace_back(
                window
            );
        }

        if (
            primers.size() !=
            target_pool
        ) {
            throw std::runtime_error(
                "Could not collect 4096 primers."
            );
        }

        /*
         * One packed reference copy.
         */
        const primerpair::PackedReference
            packed_reference(
                reference
            );

        primerpair::BidirectionalFMIndex
            index(
                std::move(reference)
            );

        const primerpair::ApproximateAnchorSearcher
            branching(
                index
            );

        const primerpair::AnchorCandidateSearcher
            candidate_searcher(
                index,
                packed_reference
            );

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

            const std::uint64_t
                anchor_count =
                    index.search(
                        anchor
                    ).size();

            PrimerCase item{
                primer,
                anchor_count
            };

            switch (
                classify(
                    anchor_count
                )
            ) {

                case Difficulty::Easy:
                    bins.at(0).push_back(
                        std::move(item)
                    );
                    break;

                case Difficulty::Moderate:
                    bins.at(1).push_back(
                        std::move(item)
                    );
                    break;

                case Difficulty::Hard:
                    bins.at(2).push_back(
                        std::move(item)
                    );
                    break;

                case Difficulty::RepeatRich:
                    bins.at(3).push_back(
                        std::move(item)
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
            << "packed_reference_bytes\t"
            << packed_reference.memory_bytes()
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
            << "budget"
            << '\t'
            << "branch_search_only_ns"
            << '\t'
            << "branch_e2e_ns"
            << '\t'
            << "candidate_e2e_ns"
            << '\t'
            << "branch_e2e_over_candidate"
            << '\t'
            << "mean_output_hits"
            << '\t'
            << "candidate_mean_verified"
            << '\t'
            << "equivalent"
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

            long double anchor_sum =
                0.0;

            for (const auto& item : cases) {
                anchor_sum +=
                    item.anchor_count;
            }

            const double mean_anchor =
                static_cast<double>(
                    anchor_sum /
                    cases.size()
                );

            /*
             * Warm-up.
             */
            for (std::size_t i = 0;
                 i < 50;
                 ++i) {

                const auto& item =
                    cases.at(
                        i %
                        cases.size()
                    );

                const auto branch_result =
                    branching
                        .search_5prime_mismatches(
                            item.primer,
                            anchor_length,
                            3
                        );

                const auto branch_hits =
                    materialize_branching(
                        index,
                        branch_result
                    );

                const auto candidate_result =
                    candidate_searcher.search(
                        item.primer,
                        anchor_length,
                        3
                    );

                checksum +=
                    branch_hits.size();

                checksum +=
                    candidate_result.hits.size();
            }

            for (std::size_t budget = 0;
                 budget <= 3;
                 ++budget) {

                /*
                 * ------------------------------------------
                 * 1. Branching search-only
                 * ------------------------------------------
                 */

                std::uint64_t
                    search_only_checksum = 0;

                const auto search_start =
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
                        branching
                            .search_5prime_mismatches(
                                item.primer,
                                anchor_length,
                                budget
                            );

                    search_only_checksum +=
                        result.total_match_count();

                    search_only_checksum +=
                        static_cast<std::uint64_t>(
                            result.hits.size()
                        );
                }

                const auto search_stop =
                    std::chrono::steady_clock::now();

                /*
                 * ------------------------------------------
                 * 2. Branching end-to-end:
                 *
                 * search + locate + sort + deduplicate
                 * ------------------------------------------
                 */

                std::uint64_t
                    branch_output_hits = 0;

                const auto branch_e2e_start =
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
                        branching
                            .search_5prime_mismatches(
                                item.primer,
                                anchor_length,
                                budget
                            );

                    const auto hits =
                        materialize_branching(
                            index,
                            result
                        );

                    branch_output_hits +=
                        static_cast<std::uint64_t>(
                            hits.size()
                        );
                }

                const auto branch_e2e_stop =
                    std::chrono::steady_clock::now();

                /*
                 * ------------------------------------------
                 * 3. Candidate end-to-end
                 * ------------------------------------------
                 */

                std::uint64_t
                    candidate_output_hits = 0;

                std::uint64_t
                    candidate_verified = 0;

                const auto candidate_start =
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
                        candidate_searcher.search(
                            item.primer,
                            anchor_length,
                            budget
                        );

                    candidate_output_hits +=
                        static_cast<std::uint64_t>(
                            result.hits.size()
                        );

                    candidate_verified +=
                        result.candidates_verified;
                }

                const auto candidate_stop =
                    std::chrono::steady_clock::now();

                const double denominator =
                    static_cast<double>(
                        queries
                    );

                const double search_only_ns =
                    static_cast<double>(
                        std::chrono::duration_cast<
                            std::chrono::nanoseconds
                        >(
                            search_stop -
                            search_start
                        ).count()
                    ) /
                    denominator;

                const double branch_e2e_ns =
                    static_cast<double>(
                        std::chrono::duration_cast<
                            std::chrono::nanoseconds
                        >(
                            branch_e2e_stop -
                            branch_e2e_start
                        ).count()
                    ) /
                    denominator;

                const double candidate_e2e_ns =
                    static_cast<double>(
                        std::chrono::duration_cast<
                            std::chrono::nanoseconds
                        >(
                            candidate_stop -
                            candidate_start
                        ).count()
                    ) /
                    denominator;

                const bool equivalent =
                    branch_output_hits ==
                    candidate_output_hits;

                checksum +=
                    search_only_checksum;

                checksum +=
                    branch_output_hits;

                checksum +=
                    candidate_output_hits;

                checksum +=
                    candidate_verified;

                std::cout
                    << difficulty_name(
                           difficulties.at(b)
                       )
                    << '\t'
                    << cases.size()
                    << '\t'
                    << mean_anchor
                    << '\t'
                    << budget
                    << '\t'
                    << search_only_ns
                    << '\t'
                    << branch_e2e_ns
                    << '\t'
                    << candidate_e2e_ns
                    << '\t'
                    << (
                        branch_e2e_ns /
                        candidate_e2e_ns
                    )
                    << '\t'
                    << (
                        static_cast<double>(
                            branch_output_hits
                        ) /
                        denominator
                    )
                    << '\t'
                    << (
                        static_cast<double>(
                            candidate_verified
                        ) /
                        denominator
                    )
                    << '\t'
                    << (
                        equivalent
                            ? "YES"
                            : "NO"
                    )
                    << '\n';

                if (!equivalent) {
                    throw std::logic_error(
                        "Branching and candidate "
                        "output counts differ."
                    );
                }
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
