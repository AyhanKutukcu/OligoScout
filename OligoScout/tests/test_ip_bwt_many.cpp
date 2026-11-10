#include "primerpair/ip_bwt_index.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>


namespace {

std::string
make_reference(
    const std::size_t length
) {
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
        length
    );


    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL
            +
            1442695040888963407ULL;


        reference.push_back(
            alphabet[
                static_cast<std::size_t>(
                    (
                        state >>
                        32
                    )
                    &
                    3ULL
                )
            ]
        );
    }


    /*
     * Add deliberately repetitive regions so the
     * many-query API also exercises intervals with
     * multiple matching rows.
     */
    const std::string motif =
        "ACGTACGTACGTACGT"
        "AAAAAAAACCCCCCCC"
        "GGGGGGGGTTTTTTTT"
        "ACACACACGTGTGTGT";


    for (
        std::size_t start = 4096;
        start + 512 < reference.size();
        start += 8192
    ) {
        for (
            std::size_t offset = 0;
            offset < 512;
            ++offset
        ) {
            reference[
                start +
                offset
            ] =
                motif[
                    offset %
                    motif.size()
                ];
        }
    }


    return reference;
}


void
check_equal(
    const primerpair::Interval expected,
    const primerpair::Interval observed,
    const std::size_t query_index
) {
    if (
        expected.begin !=
            observed.begin
        ||
        expected.end !=
            observed.end
    ) {
        std::cerr
            << "MISMATCH"
            << '\t'
            << "query_index="
            << query_index
            << '\t'
            << "expected="
            << expected.begin
            << ':'
            << expected.end
            << '\t'
            << "observed="
            << observed.begin
            << ':'
            << observed.end
            << '\n';


        throw std::runtime_error(
            "Many-query interval mismatch."
        );
    }
}


void
run_size_case(
    const primerpair::IPBWTIndex& index,
    const std::string& reference,
    const std::size_t query_count,
    std::size_t& checks
) {
    std::vector<std::string>
        owned;

    owned.reserve(
        query_count
    );


    for (
        std::size_t i = 0;
        i < query_count;
        ++i
    ) {
        /*
         * Mix ordinary genomic positions with
         * repetitive-region positions.
         */
        std::size_t position = 0;


        if (
            i %
            3 ==
            0
        ) {
            position =
                (
                    i *
                    997
                    +
                    4096
                )
                %
                (
                    reference.size()
                    -
                    21
                );

        } else {

            position =
                (
                    i *
                    7919
                    +
                    123
                )
                %
                (
                    reference.size()
                    -
                    21
                );
        }


        owned.push_back(
            reference.substr(
                position,
                21
            )
        );


        /*
         * Periodically make a deterministic query
         * unlikely to occur in the reference.
         */
        if (
            i %
            7 ==
            6
        ) {
            owned.back()[10] =
                (
                    owned.back()[10]
                    ==
                    'A'
                )
                ?
                'T'
                :
                'A';
        }
    }


    std::vector<std::string_view>
        queries;

    queries.reserve(
        owned.size()
    );


    for (
        const auto& query :
        owned
    ) {
        queries.emplace_back(
            query
        );
    }


    const auto observed =
        index.exact_prefix_search_many(
            queries
        );


    std::vector<primerpair::Interval>
        observed_preallocated(
            queries.size()
        );


    index.exact_prefix_search_many(
        std::span<const std::string_view>(
            queries.data(),
            queries.size()
        ),
        std::span<primerpair::Interval>(
            observed_preallocated.data(),
            observed_preallocated.size()
        )
    );


    if (
        observed.size() !=
        queries.size()
    ) {
        throw std::runtime_error(
            "Many-query result size mismatch."
        );
    }


    for (
        std::size_t i = 0;
        i < queries.size();
        ++i
    ) {
        const auto expected =
            index.exact_search(
                queries[i]
            );


        check_equal(
            expected,
            observed[i],
            i
        );


        check_equal(
            expected,
            observed_preallocated[i],
            i
        );


        ++checks;
    }
}


}  // namespace


int
main() {
    using namespace primerpair;


    try {
        const std::string reference =
            make_reference(
                200000
            );


        IPBWTIndex index(
            reference,
            21
        );


        std::size_t checks = 0;


        /*
         * Exercise all important batch/tail
         * boundaries.
         */
        for (
            const std::size_t query_count :
            {
                std::size_t{0},
                std::size_t{1},
                std::size_t{2},
                std::size_t{7},
                std::size_t{8},
                std::size_t{9},
                std::size_t{15},
                std::size_t{16},
                std::size_t{17},
                std::size_t{31},
                std::size_t{32},
                std::size_t{33},
                std::size_t{64},
                std::size_t{65}
            }
        ) {
            run_size_case(
                index,
                reference,
                query_count,
                checks
            );
        }


        /*
         * Invalid length in a scalar tail.
         */
        {
            std::string good =
                reference.substr(
                    100,
                    21
                );

            std::string bad =
                reference.substr(
                    200,
                    20
                );


            std::vector<std::string_view>
            queries{
                good,
                bad
            };


            bool threw = false;


            try {
                static_cast<void>(
                    index
                        .exact_prefix_search_many(
                            queries
                        )
                );

            } catch (
                const std::invalid_argument&
            ) {
                threw = true;
            }


            if (!threw) {
                throw std::runtime_error(
                    "Invalid scalar-tail length "
                    "was not rejected."
                );
            }


            ++checks;
        }


        /*
         * Invalid length inside a full Batch8.
         */
        {
            std::vector<std::string>
                owned;


            for (
                std::size_t i = 0;
                i < 8;
                ++i
            ) {
                owned.push_back(
                    reference.substr(
                        1000 +
                        i *
                        100,
                        21
                    )
                );
            }


            owned[4].resize(
                20
            );


            std::vector<std::string_view>
                queries;


            for (
                const auto& query :
                owned
            ) {
                queries.emplace_back(
                    query
                );
            }


            bool threw = false;


            try {
                static_cast<void>(
                    index
                        .exact_prefix_search_many(
                            queries
                        )
                );

            } catch (
                const std::invalid_argument&
            ) {
                threw = true;
            }


            if (!threw) {
                throw std::runtime_error(
                    "Invalid Batch8 length "
                    "was not rejected."
                );
            }


            ++checks;
        }


        /*
         * Allocation-free API must reject
         * mismatched output sizes.
         */
        {
            std::string query =
                reference.substr(
                    500,
                    21
                );


            std::vector<std::string_view>
            queries{
                query
            };


            std::vector<Interval>
                results;


            bool threw = false;


            try {
                index.exact_prefix_search_many(
                    std::span<const std::string_view>(
                        queries.data(),
                        queries.size()
                    ),
                    std::span<Interval>(
                        results.data(),
                        results.size()
                    )
                );

            } catch (
                const std::invalid_argument&
            ) {
                threw = true;
            }


            if (!threw) {
                throw std::runtime_error(
                    "Mismatched result span "
                    "was not rejected."
                );
            }


            ++checks;
        }


        std::cout
            << "checks\t"
            << checks
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
