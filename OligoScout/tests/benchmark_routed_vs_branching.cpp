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
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"

namespace {

struct PrimerCase {
    std::string primer;
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

std::size_t difficulty_index(
    const primerpair::SearchDifficulty difficulty
) {
    switch (difficulty) {

        case primerpair::SearchDifficulty::Easy:
            return 0;

        case primerpair::SearchDifficulty::Moderate:
            return 1;

        case primerpair::SearchDifficulty::Hard:
            return 2;

        case primerpair::SearchDifficulty::RepeatRich:
            return 3;
    }

    throw std::logic_error(
        "Unknown difficulty."
    );
}

const char* difficulty_name(
    const std::size_t index
) {
    constexpr std::array<const char*, 4>
        names{
            "EASY",
            "MODERATE",
            "HARD",
            "REPEAT_RICH"
        };

    return names.at(index);
}

std::vector<primerpair::PrimerSearchHit>
materialize_branching(
    const primerpair::BidirectionalFMIndex& index,
    const primerpair::ApproximateAnchorSearchResult& result
) {
    std::vector<primerpair::PrimerSearchHit> hits;

    hits.reserve(
        static_cast<std::size_t>(
            result.total_match_count()
        )
    );

    for (const auto& branch : result.hits) {

        const auto positions =
            index.locate(
                branch.state
            );

        for (const auto position : positions) {

            hits.push_back(
                {
                    position,
                    branch.mismatches
                }
            );
        }
    }

    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const auto& lhs,
            const auto& rhs
        ) {
            if (
                lhs.position !=
                rhs.position
            ) {
                return
                    lhs.position <
                    rhs.position;
            }

            return
                lhs.mismatches <
                rhs.mismatches;
        }
    );

    hits.erase(
        std::unique(
            hits.begin(),
            hits.end()
        ),
        hits.end()
    );

    return hits;
}

void update_checksum(
    std::uint64_t& checksum,
    const std::vector<
        primerpair::PrimerSearchHit
    >& hits
) {
    constexpr std::uint64_t prime =
        1099511628211ULL;

    checksum ^=
        static_cast<std::uint64_t>(
            hits.size()
        );

    checksum *= prime;

    for (const auto& hit : hits) {

        checksum ^= hit.position;
        checksum *= prime;

        checksum ^=
            static_cast<std::uint64_t>(
                hit.mismatches + 1
            );

        checksum *= prime;
    }
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

        constexpr std::size_t primer_length = 20;
        constexpr std::size_t anchor_length = 12;
        constexpr std::size_t target_pool = 4096;

        constexpr std::uint64_t stride =
            104729ULL;

        const std::size_t queries =
            argc == 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 500;

        if (queries == 0) {
            throw std::invalid_argument(
                "Queries must be > 0."
            );
        }

        std::string reference =
            load_fasta(argv[1]);

        const std::size_t available =
            reference.size() -
            primer_length +
            1;

        std::vector<std::string> primer_pool;

        primer_pool.reserve(
            target_pool
        );

        std::uint64_t candidate = 0;

        for (
            std::size_t attempts = 0;
            attempts < available &&
            primer_pool.size() < target_pool;
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

            primer_pool.emplace_back(
                window
            );
        }

        if (
            primer_pool.size() !=
            target_pool
        ) {
            throw std::runtime_error(
                "Could not collect 4096 primers."
            );
        }

        const primerpair::PackedReference
            packed_reference(
                reference
            );

        const primerpair::BidirectionalFMIndex
            index(
                reference
            );

        const primerpair::ApproximateAnchorSearcher
            branching(
                index
            );

        const primerpair::SinglePrimerSearchEngine
            routed_engine(
                index,
                packed_reference
            );

        const primerpair::SearchDifficultyEstimator
            estimator(
                index
            );

        std::array<
            std::vector<PrimerCase>,
            4
        > bins;

        for (const auto& primer : primer_pool) {

            const auto profile =
                estimator.estimate(
                    primer,
                    anchor_length
                );

            bins.at(
                difficulty_index(
                    profile.difficulty
                )
            ).push_back(
                PrimerCase{
                    primer
                }
            );
        }

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "primer_pool\t"
            << primer_pool.size()
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
            << "budget"
            << '\t'
            << "always_branching_ns"
            << '\t'
            << "routed_ns"
            << '\t'
            << "branch_over_routed_speedup"
            << '\t'
            << "direct_routes"
            << '\t'
            << "candidate_routes"
            << '\t'
            << "mean_hits"
            << '\t'
            << "equivalent"
            << '\n';

        for (std::size_t b = 0;
             b < bins.size();
             ++b) {

            const auto& cases =
                bins.at(b);

            if (cases.empty()) {
                continue;
            }

            /*
             * Warm-up.
             */
            for (std::size_t i = 0;
                 i < 50;
                 ++i) {

                const auto& primer =
                    cases.at(
                        i %
                        cases.size()
                    ).primer;

                const auto branch =
                    branching
                        .search_5prime_mismatches(
                            primer,
                            anchor_length,
                            3
                        );

                const auto branch_hits =
                    materialize_branching(
                        index,
                        branch
                    );

                const auto routed =
                    routed_engine.search(
                        primer,
                        anchor_length,
                        3
                    );

                static_cast<void>(
                    branch_hits
                );

                static_cast<void>(
                    routed
                );
            }

            for (std::size_t budget = 0;
                 budget <= 3;
                 ++budget) {

                std::uint64_t
                    branching_checksum =
                        1469598103934665603ULL;

                std::uint64_t
                    routed_checksum =
                        1469598103934665603ULL;

                std::uint64_t
                    branching_total_hits = 0;

                std::uint64_t
                    routed_total_hits = 0;

                /*
                 * --------------------------------------
                 * Always-branching E2E
                 * --------------------------------------
                 */

                const auto branch_start =
                    std::chrono::steady_clock::now();

                for (std::size_t i = 0;
                     i < queries;
                     ++i) {

                    const auto& primer =
                        cases.at(
                            i %
                            cases.size()
                        ).primer;

                    const auto result =
                        branching
                            .search_5prime_mismatches(
                                primer,
                                anchor_length,
                                budget
                            );

                    const auto hits =
                        materialize_branching(
                            index,
                            result
                        );

                    branching_total_hits +=
                        static_cast<std::uint64_t>(
                            hits.size()
                        );

                    update_checksum(
                        branching_checksum,
                        hits
                    );
                }

                const auto branch_stop =
                    std::chrono::steady_clock::now();

                /*
                 * --------------------------------------
                 * Routed E2E
                 * --------------------------------------
                 */

                std::uint64_t direct_routes = 0;
                std::uint64_t candidate_routes = 0;

                const auto routed_start =
                    std::chrono::steady_clock::now();

                for (std::size_t i = 0;
                     i < queries;
                     ++i) {

                    const auto& primer =
                        cases.at(
                            i %
                            cases.size()
                        ).primer;

                    const auto result =
                        routed_engine.search(
                            primer,
                            anchor_length,
                            budget
                        );

                    routed_total_hits +=
                        static_cast<std::uint64_t>(
                            result.hits.size()
                        );

                    update_checksum(
                        routed_checksum,
                        result.hits
                    );

                    switch (
                        result.decision
                            .recommended_strategy
                    ) {

                        case primerpair::SearchStrategy::
                            DirectBranching:

                            ++direct_routes;
                            break;

                        case primerpair::SearchStrategy::
                            AnchorCandidateVerification:

                            ++candidate_routes;
                            break;

                        case primerpair::SearchStrategy::
                            SplitSeedCandidate:

                            throw std::logic_error(
                                "Unexpected split-seed route."
                            );
                    }
                }

                const auto routed_stop =
                    std::chrono::steady_clock::now();

                const double denominator =
                    static_cast<double>(
                        queries
                    );

                const double branching_ns =
                    static_cast<double>(
                        std::chrono::duration_cast<
                            std::chrono::nanoseconds
                        >(
                            branch_stop -
                            branch_start
                        ).count()
                    ) /
                    denominator;

                const double routed_ns =
                    static_cast<double>(
                        std::chrono::duration_cast<
                            std::chrono::nanoseconds
                        >(
                            routed_stop -
                            routed_start
                        ).count()
                    ) /
                    denominator;

                const bool equivalent =
                    branching_total_hits ==
                        routed_total_hits &&
                    branching_checksum ==
                        routed_checksum;

                std::cout
                    << difficulty_name(b)
                    << '\t'
                    << cases.size()
                    << '\t'
                    << budget
                    << '\t'
                    << branching_ns
                    << '\t'
                    << routed_ns
                    << '\t'
                    << (
                        branching_ns /
                        routed_ns
                    )
                    << '\t'
                    << direct_routes
                    << '\t'
                    << candidate_routes
                    << '\t'
                    << (
                        static_cast<double>(
                            routed_total_hits
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
                        "Routed output differs from "
                        "always-branching baseline."
                    );
                }
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
