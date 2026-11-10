#include "primerpair/fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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


void compare_query(
    const primerpair::FMIndex& fm,
    const primerpair::IPBWTIndex& ip,
    const std::string& query
) {
    const auto fm_interval =
        fm.backward_search(
            query
        );

    const auto ip_interval =
        ip.exact_search(
            query
        );

    if (
        fm_interval.begin !=
            ip_interval.begin ||
        fm_interval.end !=
            ip_interval.end
    ) {
        std::cerr
            << "QUERY\t"
            << query
            << '\n'
            << "FM\t"
            << fm_interval.begin
            << '\t'
            << fm_interval.end
            << '\n'
            << "IPBWT\t"
            << ip_interval.begin
            << '\t'
            << ip_interval.end
            << '\n';

        throw std::runtime_error(
            "FM and IP-BWT intervals differ."
        );
    }


    auto fm_positions =
        fm.locate(
            fm_interval
        );

    auto ip_positions =
        ip.locate(
            ip_interval
        );

    std::sort(
        fm_positions.begin(),
        fm_positions.end()
    );

    std::sort(
        ip_positions.begin(),
        ip_positions.end()
    );

    if (
        fm_positions !=
        ip_positions
    ) {
        std::cerr
            << "QUERY\t"
            << query
            << '\n';

        throw std::runtime_error(
            "FM and IP-BWT locations differ."
        );
    }
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string reference =
            "CATTATTAGGA";

        FMIndex fm(
            reference
        );

        IPBWTIndex ip(
            reference,
            3
        );


        expect(
            ip.chunk_length() == 3,
            "IP-BWT chunk length is three"
        );

        expect(
            ip.row_count() ==
                reference.size() + 1,
            "IP-BWT contains sentinel row"
        );


        const std::vector<std::string>
            hand_queries{
                "A",
                "AT",
                "ATT",
                "ATTA",
                "TTA",
                "GGA",
                "CATT",
                "TAG",
                "GG",
                "CCC"
            };


        for (
            const auto& query :
            hand_queries
        ) {
            compare_query(
                fm,
                ip,
                query
            );
        }

        expect(
            true,
            "Hand-selected FM/IP-BWT queries agree"
        );


        /*
         * Exhaustively compare every substring
         * of the reference up to length six.
         */
        std::size_t compared = 0;

        for (
            std::size_t start = 0;
            start < reference.size();
            ++start
        ) {
            for (
                std::size_t length = 1;
                length <= 6 &&
                start + length <=
                    reference.size();
                ++length
            ) {
                compare_query(
                    fm,
                    ip,
                    reference.substr(
                        start,
                        length
                    )
                );

                ++compared;
            }
        }


        expect(
            compared > 0,
            "Exhaustive substring comparisons completed"
        );


        const auto atta =
            ip.exact_search(
                "ATTA"
            );

        const auto atta_positions =
            ip.locate(
                atta
            );

        expect(
            atta.size() == 2,
            "ATTA has two matches"
        );

        expect(
            atta_positions ==
                std::vector<std::uint64_t>{
                    1,
                    4
                },
            "ATTA genomic positions correct"
        );


        {
            const std::string long_reference =
                "ACGTACGTACGTACGTACGTACGT"
                "TGCATGCATGCATGCATGCATGCA";

            FMIndex fm_k21(
                long_reference
            );

            IPBWTIndex ip_k21(
                long_reference,
                21
            );

            compare_query(
                fm_k21,
                ip_k21,
                "ACGTACGTACGTACGTACGTA"
            );

            expect(
                ip_k21.chunk_length() == 21,
                "Twenty-one-symbol packed prefix supported"
            );

            expect(
                ip_k21.compact_storage_bytes() > 0,
                "Compact numeric storage allocated"
            );
        }


        bool oversized_chunk_rejected =
            false;

        try {
            IPBWTIndex invalid(
                "ACGTACGTACGTACGTACGTACGT",
                22
            );

            (void)invalid;

        } catch (
            const std::invalid_argument&
        ) {
            oversized_chunk_rejected =
                true;
        }

        expect(
            oversized_chunk_rejected,
            "Twenty-two-symbol packed prefix rejected"
        );


        bool zero_chunk_rejected =
            false;

        try {
            IPBWTIndex invalid(
                reference,
                0
            );

            (void)invalid;

        } catch (
            const std::invalid_argument&
        ) {
            zero_chunk_rejected =
                true;
        }

        expect(
            zero_chunk_rejected,
            "Zero chunk length rejected"
        );


        {
            const std::string n_reference =
                "ACGTNNNNACGT";

            FMIndex fm_with_n(
                n_reference
            );

            IPBWTIndex ip_with_n(
                n_reference,
                3
            );

            compare_query(
                fm_with_n,
                ip_with_n,
                "ACG"
            );

            compare_query(
                fm_with_n,
                ip_with_n,
                "CGT"
            );

            expect(
                ip_with_n.row_count() ==
                    n_reference.size() + 1,
                "Reference N symbols supported"
            );
        }


        bool unsupported_symbol_rejected =
            false;

        try {
            IPBWTIndex invalid(
                "ACGTRACGT",
                3
            );

            (void)invalid;

        } catch (
            const std::invalid_argument&
        ) {
            unsupported_symbol_rejected =
                true;
        }

        expect(
            unsupported_symbol_rejected,
            "Unsupported reference symbol rejected"
        );


        std::cout
            << "queries_compared\t"
            << (
                compared +
                hand_queries.size()
            )
            << '\n';

        std::cout
            << "chunk_length\t"
            << ip.chunk_length()
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
