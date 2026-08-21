#include "primerpair/fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/ip_bwt_rmi.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
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


void generate_kmers(
    const std::size_t remaining,
    std::string& current,
    std::vector<std::string>& output
) {
    if (
        remaining == 0
    ) {
        output.push_back(
            current
        );

        return;
    }


    constexpr char alphabet[] = {
        'A',
        'C',
        'G',
        'T'
    };


    for (
        const char base :
        alphabet
    ) {
        current.push_back(
            base
        );

        generate_kmers(
            remaining - 1,
            current,
            output
        );

        current.pop_back();
    }
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string reference =
            "CATTATTAGGA"
            "ACGTACGTACGT"
            "GATTACAGATTACA"
            "TTTTACGTAAAA"
            "ACACACACACAC";


        FMIndex fm(
            reference
        );


        IPBWTIndex ip(
            reference,
            3
        );


        IPBWTRMI rmi(
            ip,
            8,
            4
        );


        expect(
            rmi.leaf_count() == 8,
            "Eight RMI leaf models created"
        );


        expect(
            rmi.correction_margin() == 4,
            "RMI correction margin preserved"
        );


        /*
         * ----------------------------------------------------
         * Every A/C/G/T 3-mer × representative IP-BWT rows.
         * ----------------------------------------------------
         */

        std::vector<std::string>
            kmers;

        std::string current;

        generate_kmers(
            3,
            current,
            kmers
        );


        const std::uint64_t n =
            ip.row_count();


        std::vector<std::uint64_t>
            rows{
                0,
                1,
                n / 4,
                n / 2,
                (3 * n) / 4,
                n > 0
                    ? n - 1
                    : 0,
                n
            };


        std::size_t lower_bound_checks =
            0;

        std::size_t local_correction_checks =
            0;

        std::size_t global_fallback_checks =
            0;


        for (
            const auto& kmer :
            kmers
        ) {
            for (
                const std::uint64_t row :
                rows
            ) {
                const auto binary =
                    ip.lower_bound(
                        kmer,
                        row
                    );


                const auto diagnostic =
                    rmi.lower_bound_diagnostics(
                        kmer,
                        row
                    );

                const auto learned =
                    diagnostic.position;


                const auto galloping =
                    rmi.lower_bound_galloping(
                        kmer,
                        row
                    );


                if (
                    diagnostic.used_global_fallback
                ) {
                    ++global_fallback_checks;
                } else {
                    ++local_correction_checks;
                }


                if (
                    binary !=
                        learned
                    ||
                    binary !=
                        galloping
                ) {
                    std::cerr
                        << "LOWER_BOUND_MISMATCH\t"
                        << kmer
                        << '\t'
                        << row
                        << '\t'
                        << binary
                        << '\t'
                        << learned
                        << '\n';

                    throw std::runtime_error(
                        "RMI corrected lower_bound "
                        "differs from binary search."
                    );
                }


                ++lower_bound_checks;
            }
        }


        expect(
            lower_bound_checks ==
                kmers.size() *
                rows.size(),
            "RMI lower-bound matrix completed"
        );


        /*
         * Sentinel-bearing keys are required by
         * partial final chunks in IP-BWT search.
         */
        const std::vector<std::string>
            sentinel_chunks{
                "$AA",
                "A$A",
                "AC$",
                "T$A"
            };


        for (
            const auto& chunk :
            sentinel_chunks
        ) {
            for (
                const std::uint64_t row :
                rows
            ) {
                expect(
                    ip.lower_bound(
                        chunk,
                        row
                    )
                    ==
                    rmi.lower_bound(
                        chunk,
                        row
                    ),
                    "Sentinel compound-key lookup agrees"
                );
            }
        }


        /*
         * ----------------------------------------------------
         * End-to-end exact-search equivalence:
         *
         * FM
         * IP-BWT binary
         * IP-BWT RMI
         * ----------------------------------------------------
         */

        std::size_t query_checks = 0;


        for (
            std::size_t start = 0;
            start < reference.size();
            ++start
        ) {
            for (
                std::size_t length = 1;
                length <= 12 &&
                start + length <=
                    reference.size();
                ++length
            ) {
                const std::string query =
                    reference.substr(
                        start,
                        length
                    );


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


                const auto galloping_interval =
                    rmi.exact_search_galloping(
                        query
                    );


                if (
                    fm_interval.begin !=
                        binary_interval.begin
                    ||
                    fm_interval.end !=
                        binary_interval.end
                    ||
                    fm_interval.begin !=
                        learned_interval.begin
                    ||
                    fm_interval.end !=
                        learned_interval.end
                    ||
                    binary_interval.begin !=
                        galloping_interval.begin
                    ||
                    binary_interval.end !=
                        galloping_interval.end
                ) {
                    std::cerr
                        << "QUERY_INTERVAL_MISMATCH\t"
                        << query
                        << '\n';

                    throw std::runtime_error(
                        "FM / binary IP-BWT / "
                        "RMI IP-BWT intervals differ."
                    );
                }


                const auto fm_positions =
                    fm.locate(
                        fm_interval
                    );


                const auto learned_positions =
                    ip.locate(
                        learned_interval
                    );


                if (
                    fm_positions !=
                    learned_positions
                ) {
                    std::cerr
                        << "QUERY_LOCATE_MISMATCH\t"
                        << query
                        << '\n';

                    throw std::runtime_error(
                        "FM and RMI locations differ."
                    );
                }


                ++query_checks;
            }
        }


        expect(
            query_checks > 0,
            "End-to-end RMI exact-search checks completed"
        );


        /*
         * K=21 is the target packed-prefix configuration.
         */
        {
            const std::string long_reference =
                "ACGTACGTACGTACGTACGTACGT"
                "TGCATGCATGCATGCATGCATGCA"
                "GATTACAGATTACAGATTACAGATTACA";


            IPBWTIndex ip21(
                long_reference,
                21
            );


            IPBWTRMI rmi21(
                ip21,
                8,
                4
            );


            const std::string query =
                "ACGTACGTACGTACGTACGTA";


            expect(
                ip21.exact_search(
                    query
                ).begin
                ==
                rmi21.exact_search(
                    query
                ).begin,
                "K=21 RMI interval begin agrees"
            );


            expect(
                ip21.exact_search(
                    query
                ).end
                ==
                rmi21.exact_search(
                    query
                ).end,
                "K=21 RMI interval end agrees"
            );
        }


        bool zero_leaf_rejected =
            false;

        try {
            IPBWTRMI invalid(
                ip,
                0,
                4
            );

            (void)invalid;

        } catch (
            const std::invalid_argument&
        ) {
            zero_leaf_rejected =
                true;
        }


        expect(
            zero_leaf_rejected,
            "Zero RMI leaf count rejected"
        );


        std::cout
            << "lower_bound_checks\t"
            << lower_bound_checks
            << '\n';


        std::cout
            << "query_checks\t"
            << query_checks
            << '\n';


        std::cout
            << "local_correction_checks\t"
            << local_correction_checks
            << '\n';


        std::cout
            << "global_fallback_checks\t"
            << global_fallback_checks
            << '\n';


        const double fallback_rate =
            lower_bound_checks == 0
            ?
            0.0
            :
            static_cast<double>(
                global_fallback_checks
            )
            /
            static_cast<double>(
                lower_bound_checks
            );


        std::cout
            << "global_fallback_rate\t"
            << fallback_rate
            << '\n';




        std::cout
            << "rmi_leaf_count\t"
            << rmi.leaf_count()
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
