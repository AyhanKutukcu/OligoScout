#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/multiplex_primer_search.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


namespace {


using Clock =
    std::chrono::steady_clock;


std::string make_reference(
    const std::size_t length
) {
    constexpr char bases[] = {
        'A', 'C', 'G', 'T'
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
            6364136223846793005ULL +
            1442695040888963407ULL;

        reference.push_back(
            bases[
                static_cast<std::size_t>(
                    (state >> 32U) &
                    3ULL
                )
            ]
        );
    }


    return reference;
}


double ns_between(
    const Clock::time_point begin,
    const Clock::time_point end
) {
    return
        static_cast<double>(
            std::chrono::duration_cast<
                std::chrono::nanoseconds
            >(
                end - begin
            ).count()
        );
}


}  // namespace


int main(
    int argc,
    char** argv
) {
    try {
        using namespace primerpair;


        std::size_t repetitions = 20;

        if (argc >= 2) {
            repetitions =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[1]
                    )
                );
        }


        if (repetitions == 0) {
            throw std::invalid_argument(
                "Repetitions must be > 0."
            );
        }


        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            primer_length = 24;


        const std::string reference =
            make_reference(
                150000
            );


        PackedReference packed(
            reference
        );

        BidirectionalFMIndex bifm(
            reference
        );

        IPBWTIndex ipbwt(
            reference,
            anchor_length
        );

        SearchDifficultyEstimator estimator(
            bifm
        );

        BatchedAnchorLookup anchor_lookup(
            ipbwt,
            estimator
        );

        AnchorCandidateSearcher verifier(
            bifm,
            packed
        );

        BatchedCandidateSearchEngine
            candidate_engine(
                ipbwt,
                anchor_lookup,
                verifier
            );

        SinglePrimerSearchEngine
            single_engine(
                bifm,
                packed
            );

        HybridBatchedPrimerSearchEngine
            forward_hybrid(
                candidate_engine,
                single_engine
            );

        StrandAwarePrimerSearchEngine
            legacy_strand(
                bifm,
                packed
            );

        HybridStrandAwarePrimerSearchEngine
            hybrid_strand(
                forward_hybrid,
                anchor_lookup,
                ipbwt,
                packed,
                legacy_strand
            );

        MultiplexPrimerSearchEngine
            multiplex(
                hybrid_strand
            );

        PrimerPairSearchEngine
            legacy_pair(
                bifm,
                packed
            );


        /*
         * ============================================
         * Scenario A
         *
         * Shared-primer intended panel.
         * ============================================
         */

        constexpr std::size_t
            intended_pair_count = 60;


        const std::string shared =
            reference.substr(
                1000,
                primer_length
            );


        std::vector<std::string>
            right_storage;

        right_storage.reserve(
            intended_pair_count
        );


        std::vector<
            MultiplexPrimerPairRequest
        > intended_requests;

        intended_requests.reserve(
            intended_pair_count
        );


        for (
            std::size_t i = 0;
            i < intended_pair_count;
            ++i
        ) {
            const std::size_t position =
                1300 +
                i * 1700;


            right_storage.push_back(
                reverse_complement(
                    reference.substr(
                        position,
                        primer_length
                    )
                )
            );


            intended_requests.push_back(
                MultiplexPrimerPairRequest{
                    shared,
                    right_storage.back(),
                    3,
                    50,
                    3000
                }
            );
        }


        /*
         * Correctness checksum.
         */
        std::uint64_t
            legacy_intended_checksum = 0;


        for (
            const auto& request :
            intended_requests
        ) {
            const auto result =
                legacy_pair.search(
                    request.primer1,
                    request.primer2,
                    anchor_length,
                    3,
                    50,
                    3000
                );

            legacy_intended_checksum +=
                result.amplicons.size();
        }


        const auto multiplex_once =
            multiplex.search(
                intended_requests,
                anchor_length,
                false
            );


        std::uint64_t
            multiplex_intended_checksum = 0;


        for (
            const auto& pair :
            multiplex_once.intended_pairs
        ) {
            multiplex_intended_checksum +=
                pair.amplicons.size();
        }


        if (
            legacy_intended_checksum !=
            multiplex_intended_checksum
        ) {
            throw std::runtime_error(
                "Intended benchmark checksum mismatch."
            );
        }


        std::uint64_t sink = 0;


        const auto legacy_begin =
            Clock::now();


        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            for (
                const auto& request :
                intended_requests
            ) {
                const auto result =
                    legacy_pair.search(
                        request.primer1,
                        request.primer2,
                        anchor_length,
                        3,
                        50,
                        3000
                    );

                sink +=
                    result.amplicons.size();
            }
        }


        const auto legacy_end =
            Clock::now();


        const auto multiplex_begin =
            Clock::now();


        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            const auto result =
                multiplex.search(
                    intended_requests,
                    anchor_length,
                    false
                );

            for (
                const auto& pair :
                result.intended_pairs
            ) {
                sink +=
                    pair.amplicons.size();
            }
        }


        const auto multiplex_end =
            Clock::now();


        const double
            legacy_intended_ns =
                ns_between(
                    legacy_begin,
                    legacy_end
                );


        const double
            multiplex_intended_ns =
                ns_between(
                    multiplex_begin,
                    multiplex_end
                );


        /*
         * ============================================
         * Scenario B
         *
         * Cross-join reuse.
         * ============================================
         */

        constexpr std::size_t
            cross_pair_count = 12;


        const std::string primer_a =
            reference.substr(
                10000,
                primer_length
            );


        const std::string primer_b =
            reverse_complement(
                reference.substr(
                    10400,
                    primer_length
                )
            );


        std::vector<
            MultiplexPrimerPairRequest
        > cross_requests;


        for (
            std::size_t i = 0;
            i < cross_pair_count;
            ++i
        ) {
            cross_requests.push_back(
                MultiplexPrimerPairRequest{
                    primer_a,
                    primer_b,
                    3,
                    50,
                    1000
                }
            );
        }


        constexpr std::size_t
            cross_slot_requests =
                (
                    cross_pair_count *
                    (
                        cross_pair_count - 1
                    ) /
                    2
                ) *
                4;


        const std::array<
            std::string_view,
            2
        > cross_sequences{
            primer_a,
            primer_b
        };


        std::uint64_t
            legacy_cross_checksum = 0;


        for (
            std::size_t pair_a = 0;
            pair_a < cross_pair_count;
            ++pair_a
        ) {
            for (
                std::size_t pair_b =
                    pair_a + 1;

                pair_b < cross_pair_count;
                ++pair_b
            ) {
                for (
                    std::size_t slot_a = 0;
                    slot_a < 2;
                    ++slot_a
                ) {
                    for (
                        std::size_t slot_b = 0;
                        slot_b < 2;
                        ++slot_b
                    ) {
                        const auto product =
                            legacy_pair.search(
                                cross_sequences.at(
                                    slot_a
                                ),
                                cross_sequences.at(
                                    slot_b
                                ),
                                anchor_length,
                                3,
                                50,
                                1000
                            );

                        legacy_cross_checksum +=
                            product
                                .amplicons
                                .size();
                    }
                }
            }
        }


        const auto cross_once =
            multiplex.search(
                cross_requests,
                anchor_length,
                true,
                50,
                1000
            );


        if (
            legacy_cross_checksum !=
            cross_once
                .cross_amplicons
                .size()
        ) {
            throw std::runtime_error(
                "Cross benchmark checksum mismatch."
            );
        }


        const auto
            legacy_cross_begin =
                Clock::now();


        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            for (
                std::size_t pair_a = 0;
                pair_a < cross_pair_count;
                ++pair_a
            ) {
                for (
                    std::size_t pair_b =
                        pair_a + 1;

                    pair_b <
                        cross_pair_count;
                    ++pair_b
                ) {
                    for (
                        std::size_t slot_a = 0;
                        slot_a < 2;
                        ++slot_a
                    ) {
                        for (
                            std::size_t slot_b = 0;
                            slot_b < 2;
                            ++slot_b
                        ) {
                            const auto product =
                                legacy_pair.search(
                                    cross_sequences.at(
                                        slot_a
                                    ),
                                    cross_sequences.at(
                                        slot_b
                                    ),
                                    anchor_length,
                                    3,
                                    50,
                                    1000
                                );

                            sink +=
                                product
                                    .amplicons
                                    .size();
                        }
                    }
                }
            }
        }


        const auto legacy_cross_end =
            Clock::now();


        const auto
            multiplex_cross_begin =
                Clock::now();


        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            const auto result =
                multiplex.search(
                    cross_requests,
                    anchor_length,
                    true,
                    50,
                    1000
                );

            sink +=
                result
                    .cross_amplicons
                    .size();
        }


        const auto multiplex_cross_end =
            Clock::now();


        const double legacy_cross_ns =
            ns_between(
                legacy_cross_begin,
                legacy_cross_end
            );


        const double multiplex_cross_ns =
            ns_between(
                multiplex_cross_begin,
                multiplex_cross_end
            );


        std::cout
            << "repetitions\t"
            << repetitions
            << '\n';


        std::cout
            << "intended_pairs\t"
            << intended_pair_count
            << '\n';


        std::cout
            << "intended_primer_slots\t"
            << intended_pair_count * 2
            << '\n';


        std::cout
            << "intended_unique_primers\t"
            << multiplex_once
                .stats
                .unique_primer_queries
            << '\n';


        std::cout
            << "legacy_intended_ns_per_pair\t"
            << (
                legacy_intended_ns /
                (
                    repetitions *
                    intended_pair_count
                )
            )
            << '\n';


        std::cout
            << "multiplex_intended_ns_per_pair\t"
            << (
                multiplex_intended_ns /
                (
                    repetitions *
                    intended_pair_count
                )
            )
            << '\n';


        std::cout
            << "intended_speedup\t"
            << (
                legacy_intended_ns /
                multiplex_intended_ns
            )
            << '\n';


        std::cout
            << "cross_pairs\t"
            << cross_pair_count
            << '\n';


        std::cout
            << "cross_slot_requests\t"
            << cross_slot_requests
            << '\n';


        std::cout
            << "cross_unique_join_computations\t"
            << cross_once
                .stats
                .unique_cross_join_computations
            << '\n';


        std::cout
            << "cross_reused_join_requests\t"
            << cross_once
                .stats
                .reused_cross_join_requests
            << '\n';


        std::cout
            << "legacy_cross_ns_per_slot_request\t"
            << (
                legacy_cross_ns /
                (
                    repetitions *
                    cross_slot_requests
                )
            )
            << '\n';


        std::cout
            << "multiplex_cross_ns_per_slot_request\t"
            << (
                multiplex_cross_ns /
                (
                    repetitions *
                    cross_slot_requests
                )
            )
            << '\n';


        std::cout
            << "cross_speedup\t"
            << (
                legacy_cross_ns /
                multiplex_cross_ns
            )
            << '\n';


        std::cout
            << "VERIFY_INTENDED_CHECKSUM\tYES\n";

        std::cout
            << "VERIFY_CROSS_CHECKSUM\tYES\n";

        std::cerr
            << "sink\t"
            << sink
            << '\n';


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
