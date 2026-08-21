#include "primerpair/fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/ip_bwt_rmi.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void expect(
    const bool condition,
    const char* message
) {
    if (!condition) {
        throw std::runtime_error(
            message
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}


std::string make_reference() {
    constexpr char alphabet[] = {
        'A',
        'C',
        'G',
        'T'
    };

    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;

    std::string reference;

    reference.reserve(
        5000
    );


    for (
        std::size_t i = 0;
        i < 4096;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL
            +
            1442695040888963407ULL;

        const std::size_t index =
            static_cast<std::size_t>(
                (
                    state >>
                    32
                )
                &
                3ULL
            );

        reference.push_back(
            alphabet[index]
        );
    }


    const std::string motif =
        "ACGTACGTACGTACGTACGTACGT"
        "TTTTAAAACCCCGGGG"
        "ACGTACGTACGTACGTACGTACGT";

    for (
        std::size_t repeat = 0;
        repeat < 12;
        ++repeat
    ) {
        reference +=
            motif;
    }


    return reference;
}


char mutate_base(
    const char base
) {
    switch (base) {
        case 'A':
            return 'C';

        case 'C':
            return 'G';

        case 'G':
            return 'T';

        case 'T':
            return 'A';

        default:
            throw std::runtime_error(
                "Unexpected nucleotide."
            );
    }
}


void compare_query(
    const primerpair::FMIndex& fm,
    const primerpair::IPBWTIndex& ip,
    const primerpair::IPBWTRMI& rmi,
    const std::string& query,
    const std::size_t primer_length,
    const std::size_t start
) {
    const auto fm_interval =
        fm.backward_search(
            query
        );

    const auto binary_interval =
        ip.exact_search(
            query
        );

    const auto learned_interval =
        rmi.exact_search(
            query
        );


    /*
     * Binary IP-BWT and learned IP-BWT operate on
     * the same IP-BWT row space, therefore their
     * raw intervals must always be identical.
     */
    const bool binary_rmi_interval_equal =
        binary_interval.begin ==
            learned_interval.begin
        &&
        binary_interval.end ==
            learned_interval.end;


    /*
     * FM and IP-BWT can produce different insertion
     * coordinates for an EMPTY result because their
     * search paths/chunk boundaries differ.
     *
     * For non-empty results, however, both index
     * structures use the same suffix/BW row order,
     * so interval coordinates must agree.
     */
    const bool all_empty =
        fm_interval.empty()
        &&
        binary_interval.empty()
        &&
        learned_interval.empty();


    const bool nonempty_intervals_equal =
        fm_interval.begin ==
            binary_interval.begin
        &&
        fm_interval.end ==
            binary_interval.end
        &&
        fm_interval.begin ==
            learned_interval.begin
        &&
        fm_interval.end ==
            learned_interval.end;


    /*
     * Ultimately exact-search correctness is defined
     * by genomic hit coordinates.
     */
    const auto fm_positions =
        fm.locate(
            fm_interval
        );

    const auto binary_positions =
        ip.locate(
            binary_interval
        );

    const auto learned_positions =
        ip.locate(
            learned_interval
        );


    const bool hit_sets_equal =
        fm_positions ==
            binary_positions
        &&
        fm_positions ==
            learned_positions;


    const bool equal =
        binary_rmi_interval_equal
        &&
        hit_sets_equal
        &&
        (
            all_empty
            ||
            nonempty_intervals_equal
        );


    if (!equal) {
        std::cerr
            << "SEARCH_MISMATCH\t"
            << "length="
            << primer_length
            << '\t'
            << "start="
            << start
            << '\t'
            << "query="
            << query
            << '\n';

        std::cerr
            << "FM\t"
            << fm_interval.begin
            << '\t'
            << fm_interval.end
            << '\t'
            << "hits="
            << fm_positions.size()
            << '\n';

        std::cerr
            << "IPBWT_BINARY\t"
            << binary_interval.begin
            << '\t'
            << binary_interval.end
            << '\t'
            << "hits="
            << binary_positions.size()
            << '\n';

        std::cerr
            << "IPBWT_RMI\t"
            << learned_interval.begin
            << '\t'
            << learned_interval.end
            << '\t'
            << "hits="
            << learned_positions.size()
            << '\n';

        throw std::runtime_error(
            "Primer-length exact-search "
            "semantic mismatch."
        );
    }
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string reference =
            make_reference();


        expect(
            reference.size() == 4864,
            "Deterministic reference length is stable"
        );

        expect(
            reference.substr(
                0,
                18
            ) ==
                "CCGGCGCCCTTCAGTGTC",
            "Deterministic reference fingerprint is stable"
        );


        FMIndex fm(
            reference
        );


        /*
         * Target LISA/IP-BWT configuration.
         */
        IPBWTIndex ip(
            reference,
            21
        );


        /*
         * More leaves than the tiny correctness
         * test because this reference is larger.
         */
        IPBWTRMI rmi(
            ip,
            64,
            8
        );


        expect(
            ip.chunk_length() == 21,
            "K=21 IP-BWT configured"
        );

        expect(
            rmi.leaf_count() == 64,
            "Sixty-four RMI leaves created"
        );


        std::size_t present_checks = 0;
        std::size_t mutated_checks = 0;


        /*
         * OligoScout MVP supports
         * primer lengths 18 through 35.
         *
         * Step seven samples many different
         * coordinates and chunk phases while
         * keeping the sanitizer test lightweight.
         */
        for (
            std::size_t length = 18;
            length <= 35;
            ++length
        ) {
            for (
                std::size_t start = 0;
                start + length <=
                    reference.size();
                start += 7
            ) {
                const std::string query =
                    reference.substr(
                        start,
                        length
                    );


                compare_query(
                    fm,
                    ip,
                    rmi,
                    query,
                    length,
                    start
                );

                ++present_checks;


                /*
                 * Also test a one-base-mutated query.
                 *
                 * Whether it exists elsewhere in the
                 * reference is irrelevant:
                 *
                 * all three backends must return the
                 * same exact interval.
                 */
                std::string mutated =
                    query;

                const std::size_t mutation_position =
                    length /
                    2;

                mutated.at(
                    mutation_position
                ) =
                    mutate_base(
                        mutated.at(
                            mutation_position
                        )
                    );


                compare_query(
                    fm,
                    ip,
                    rmi,
                    mutated,
                    length,
                    start
                );

                ++mutated_checks;
            }
        }


        expect(
            present_checks > 1000,
            "Thousands of present primer queries checked"
        );

        expect(
            mutated_checks ==
                present_checks,
            "Mutated-query matrix completed"
        );


        /*
         * Explicitly cover the most important
         * K=21 boundary relationships.
         */
        const std::size_t boundary_lengths[] = {
            18,
            20,
            21,
            22,
            30,
            35
        };


        for (
            const std::size_t length :
            boundary_lengths
        ) {
            const std::string query =
                reference.substr(
                    123,
                    length
                );


            compare_query(
                fm,
                ip,
                rmi,
                query,
                length,
                123
            );
        }


        expect(
            true,
            "K=21 boundary primer lengths agree"
        );


        std::cout
            << "reference_length\t"
            << reference.size()
            << '\n';

        std::cout
            << "primer_length_min\t18\n";

        std::cout
            << "primer_length_max\t35\n";

        std::cout
            << "present_queries\t"
            << present_checks
            << '\n';

        std::cout
            << "mutated_queries\t"
            << mutated_checks
            << '\n';

        std::cout
            << "total_query_checks\t"
            << (
                present_checks +
                mutated_checks +
                6
            )
            << '\n';

        std::cout
            << "rmi_mean_leaf_training_error\t"
            << static_cast<double>(
                rmi.mean_leaf_training_error()
            )
            << '\n';

        std::cout
            << "rmi_max_leaf_training_error\t"
            << rmi.maximum_leaf_training_error()
            << '\n';

        std::cout
            << "ALL_CHECKS\tYES\n";


        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';

        return 1;
    }
}
