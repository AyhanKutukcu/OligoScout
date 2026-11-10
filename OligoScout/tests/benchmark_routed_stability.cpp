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

#include "primerpair/approximate_anchor_search.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"

namespace {

struct PrimerCase {
    std::string primer;
};

struct BatchResult {
    double ns_per_query{0.0};

    std::uint64_t checksum{0};
    std::uint64_t total_hits{0};

    std::uint64_t direct_routes{0};
    std::uint64_t candidate_routes{0};
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
    std::vector<
        primerpair::PrimerSearchHit
    > hits;

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

BatchResult run_branching(
    const primerpair::BidirectionalFMIndex& index,
    const primerpair::ApproximateAnchorSearcher& searcher,
    const std::vector<PrimerCase>& cases,
    const std::size_t queries,
    const std::size_t anchor_length,
    const std::size_t budget
) {
    BatchResult output;

    output.checksum =
        1469598103934665603ULL;

    const auto start =
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
            searcher
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

        output.total_hits +=
            static_cast<std::uint64_t>(
                hits.size()
            );

        update_checksum(
            output.checksum,
            hits
        );
    }

    const auto stop =
        std::chrono::steady_clock::now();

    output.ns_per_query =
        static_cast<double>(
            std::chrono::duration_cast<
                std::chrono::nanoseconds
            >(
                stop - start
            ).count()
        ) /
        static_cast<double>(
            queries
        );

    return output;
}

BatchResult run_routed(
    const primerpair::SinglePrimerSearchEngine& engine,
    const std::vector<PrimerCase>& cases,
    const std::size_t queries,
    const std::size_t anchor_length,
    const std::size_t budget
) {
    BatchResult output;

    output.checksum =
        1469598103934665603ULL;

    const auto start =
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
            engine.search(
                primer,
                anchor_length,
                budget
            );

        output.total_hits +=
            static_cast<std::uint64_t>(
                result.hits.size()
            );

        update_checksum(
            output.checksum,
            result.hits
        );

        switch (
            result.decision
                .recommended_strategy
        ) {

            case primerpair::SearchStrategy::
                DirectBranching:

                ++output.direct_routes;
                break;

            case primerpair::SearchStrategy::
                AnchorCandidateVerification:

                ++output.candidate_routes;
                break;

            case primerpair::SearchStrategy::
                SplitSeedCandidate:

                throw std::logic_error(
                    "Unexpected split-seed route."
                );
        }
    }

    const auto stop =
        std::chrono::steady_clock::now();

    output.ns_per_query =
        static_cast<double>(
            std::chrono::duration_cast<
                std::chrono::nanoseconds
            >(
                stop - start
            ).count()
        ) /
        static_cast<double>(
            queries
        );

    return output;
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
        if (argc < 2 || argc > 4) {

            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>"
                << " [queries_per_batch]"
                << " [repeats]\n";

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
            argc >= 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : 500;

        const std::size_t repeats =
            argc >= 4
                ? static_cast<std::size_t>(
                      std::stoull(argv[3])
                  )
                : 9;

        if (queries == 0) {
            throw std::invalid_argument(
                "Queries must be > 0."
            );
        }

        if (repeats < 3) {
            throw std::invalid_argument(
                "Use at least 3 repeats."
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

        std::vector<std::string>
            primer_pool;

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

        const primerpair::
            ApproximateAnchorSearcher
                branching(
                    index
                );

        const primerpair::
            SinglePrimerSearchEngine
                routed(
                    index,
                    packed_reference
                );

        const primerpair::
            SearchDifficultyEstimator
                estimator(
                    index
                );

        std::array<
            std::vector<PrimerCase>,
            4
        > bins;

        for (
            const auto& primer :
            primer_pool
        ) {

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
            << "queries_per_batch\t"
            << queries
            << '\n';

        std::cout
            << "repeats\t"
            << repeats
            << '\n';

        std::cout
            << "difficulty"
            << '\t'
            << "bin_size"
            << '\t'
            << "budget"
            << '\t'
            << "branch_median_ns"
            << '\t'
            << "routed_median_ns"
            << '\t'
            << "median_speedup"
            << '\t'
            << "routed_wins"
            << '\t'
            << "branch_min_ns"
            << '\t'
            << "branch_max_ns"
            << '\t'
            << "routed_min_ns"
            << '\t'
            << "routed_max_ns"
            << '\t'
            << "direct_routes"
            << '\t'
            << "candidate_routes"
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
             * Initial warm-up before timed repeats.
             */
            const std::size_t warm_queries =
                std::min<std::size_t>(
                    queries,
                    100
                );

            static_cast<void>(
                run_branching(
                    index,
                    branching,
                    cases,
                    warm_queries,
                    anchor_length,
                    3
                )
            );

            static_cast<void>(
                run_routed(
                    routed,
                    cases,
                    warm_queries,
                    anchor_length,
                    3
                )
            );

            for (std::size_t budget = 0;
                 budget <= 3;
                 ++budget) {

                std::vector<double>
                    branch_times;

                std::vector<double>
                    routed_times;

                branch_times.reserve(
                    repeats
                );

                routed_times.reserve(
                    repeats
                );

                std::size_t routed_wins =
                    0;

                bool equivalent =
                    true;

                std::uint64_t
                    direct_routes = 0;

                std::uint64_t
                    candidate_routes = 0;

                for (std::size_t repeat = 0;
                     repeat < repeats;
                     ++repeat) {

                    BatchResult branch_result;
                    BatchResult routed_result;

                    /*
                     * Alternate measurement order to
                     * reduce systematic cache/order bias.
                     */
                    if (
                        repeat %
                        2 == 0
                    ) {

                        branch_result =
                            run_branching(
                                index,
                                branching,
                                cases,
                                queries,
                                anchor_length,
                                budget
                            );

                        routed_result =
                            run_routed(
                                routed,
                                cases,
                                queries,
                                anchor_length,
                                budget
                            );

                    } else {

                        routed_result =
                            run_routed(
                                routed,
                                cases,
                                queries,
                                anchor_length,
                                budget
                            );

                        branch_result =
                            run_branching(
                                index,
                                branching,
                                cases,
                                queries,
                                anchor_length,
                                budget
                            );
                    }

                    branch_times.push_back(
                        branch_result.ns_per_query
                    );

                    routed_times.push_back(
                        routed_result.ns_per_query
                    );

                    if (
                        routed_result.ns_per_query <
                        branch_result.ns_per_query
                    ) {
                        ++routed_wins;
                    }

                    if (
                        branch_result.checksum !=
                            routed_result.checksum ||
                        branch_result.total_hits !=
                            routed_result.total_hits
                    ) {
                        equivalent = false;
                    }

                    direct_routes =
                        routed_result.direct_routes;

                    candidate_routes =
                        routed_result.candidate_routes;
                }

                const double branch_median =
                    median(
                        branch_times
                    );

                const double routed_median =
                    median(
                        routed_times
                    );

                const auto [
                    branch_min_it,
                    branch_max_it
                ] =
                    std::minmax_element(
                        branch_times.begin(),
                        branch_times.end()
                    );

                const auto [
                    routed_min_it,
                    routed_max_it
                ] =
                    std::minmax_element(
                        routed_times.begin(),
                        routed_times.end()
                    );

                std::cout
                    << difficulty_name(b)
                    << '\t'
                    << cases.size()
                    << '\t'
                    << budget
                    << '\t'
                    << branch_median
                    << '\t'
                    << routed_median
                    << '\t'
                    << (
                        branch_median /
                        routed_median
                    )
                    << '\t'
                    << routed_wins
                    << '/'
                    << repeats
                    << '\t'
                    << *branch_min_it
                    << '\t'
                    << *branch_max_it
                    << '\t'
                    << *routed_min_it
                    << '\t'
                    << *routed_max_it
                    << '\t'
                    << direct_routes
                    << '\t'
                    << candidate_routes
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
                        "always-branching output."
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
