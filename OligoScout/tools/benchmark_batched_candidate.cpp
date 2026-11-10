#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/search_strategy.hpp"

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


std::string
make_reference(
    const std::size_t length
) {
    constexpr char alphabet[] = {
        'A', 'C', 'G', 'T'
    };

    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;

    std::string reference;
    reference.reserve(length);

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
                    (state >> 32U) & 3ULL
                )
            ]
        );
    }

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
            reference[start + offset] =
                motif[
                    offset %
                    motif.size()
                ];
        }
    }

    return reference;
}


char
mutated_base(
    const char base
) {
    return
        base == 'A'
        ? 'T'
        : 'A';
}


std::uint64_t
checksum_result(
    const primerpair::
        AnchorCandidateSearchResult& result
) {
    std::uint64_t checksum =
        result.anchor_occurrences
        +
        result.candidates_verified
        +
        result.hits.size();

    for (
        const auto& hit :
        result.hits
    ) {
        checksum ^=
            (
                hit.position *
                0x9E3779B97F4A7C15ULL
            )
            ^
            static_cast<std::uint64_t>(
                hit.mismatches + 1
            );
    }

    return checksum;
}

}  // namespace


int
main(
    const int argc,
    char** argv
) {
    try {
        using namespace primerpair;

        std::size_t repetitions = 100;

        if (argc >= 2) {
            repetitions =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[1]
                    )
                );
        }

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            primer_length = 24;


        const std::string reference =
            make_reference(
                200000
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
            batch_engine(
                ipbwt,
                anchor_lookup,
                verifier
            );


        std::vector<std::string>
            owned_primers;

        owned_primers.reserve(
            4096
        );


        /*
         * Same type of workload as the
         * correctness integration test, but
         * larger.
         */
        for (
            std::size_t i = 0;
            i < 1000;
            ++i
        ) {
            const std::size_t position =
                (
                    i *
                    3571
                    +
                    123
                )
                %
                (
                    reference.size()
                    -
                    primer_length
                );

            const std::string original =
                reference.substr(
                    position,
                    primer_length
                );

            owned_primers.push_back(
                original
            );

            for (
                std::size_t introduced = 1;
                introduced <= 3;
                ++introduced
            ) {
                std::string primer =
                    original;

                for (
                    std::size_t m = 0;
                    m < introduced;
                    ++m
                ) {
                    const std::size_t p =
                        1 + m * 3;

                    primer[p] =
                        mutated_base(
                            primer[p]
                        );
                }

                owned_primers.push_back(
                    std::move(
                        primer
                    )
                );
            }
        }


        std::vector<
            BatchedCandidateRequest
        > all_requests;

        all_requests.reserve(
            owned_primers.size()
        );

        for (
            const auto& primer :
            owned_primers
        ) {
            all_requests.push_back(
                BatchedCandidateRequest{
                    primer,
                    3
                }
            );
        }


        /*
         * Build a candidate-only panel using
         * the legacy router.
         *
         * This makes the timing comparison fair:
         * both paths execute the same candidate
         * work.
         */
        SearchStrategyRouter router(
            bifm
        );


        std::vector<
            BatchedCandidateRequest
        > candidate_requests;

        candidate_requests.reserve(
            all_requests.size()
        );


        for (
            const auto& request :
            all_requests
        ) {
            const auto decision =
                router.decide(
                    request.primer,
                    anchor_length,
                    request.max_mismatches
                );

            if (
                decision.recommended_strategy ==
                SearchStrategy::
                    AnchorCandidateVerification
            ) {
                candidate_requests.push_back(
                    request
                );
            }
        }


        if (
            candidate_requests.empty()
        ) {
            throw std::runtime_error(
                "No candidate-routed primers."
            );
        }


        /*
         * Correctness checksum before timing.
         */
        std::uint64_t legacy_checksum = 0;

        for (
            const auto& request :
            candidate_requests
        ) {
            legacy_checksum ^=
                checksum_result(
                    verifier.search(
                        request.primer,
                        anchor_length,
                        request.max_mismatches
                    )
                );
        }


        const auto batch_once =
            batch_engine.search(
                candidate_requests,
                anchor_length
            );


        std::uint64_t batch_checksum = 0;

        for (
            const auto& result :
            batch_once
        ) {
            if (
                !result.candidate_executed
            ) {
                throw std::runtime_error(
                    "Candidate-only panel "
                    "was not executed."
                );
            }

            batch_checksum ^=
                checksum_result(
                    result.candidate_result
                );
        }


        if (
            legacy_checksum !=
            batch_checksum
        ) {
            throw std::runtime_error(
                "Benchmark checksum mismatch."
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
                candidate_requests
            ) {
                sink ^=
                    checksum_result(
                        verifier.search(
                            request.primer,
                            anchor_length,
                            request.max_mismatches
                        )
                    );
            }
        }

        const auto legacy_end =
            Clock::now();


        const auto batch_begin =
            Clock::now();

        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            const auto results =
                batch_engine.search(
                    candidate_requests,
                    anchor_length
                );

            for (
                const auto& result :
                results
            ) {
                sink ^=
                    checksum_result(
                        result.candidate_result
                    );
            }
        }

        const auto batch_end =
            Clock::now();


        const double query_count =
            static_cast<double>(
                candidate_requests.size()
                *
                repetitions
            );


        const double legacy_ns =
            std::chrono::duration<double, std::nano>(
                legacy_end -
                legacy_begin
            ).count()
            /
            query_count;


        const double batch_ns =
            std::chrono::duration<double, std::nano>(
                batch_end -
                batch_begin
            ).count()
            /
            query_count;


        std::cout
            << "panel_size\t"
            << candidate_requests.size()
            << '\n';

        std::cout
            << "repetitions\t"
            << repetitions
            << '\n';

        std::cout
            << "legacy_ns_per_primer\t"
            << legacy_ns
            << '\n';

        std::cout
            << "batch_ns_per_primer\t"
            << batch_ns
            << '\n';

        std::cout
            << "speedup\t"
            << legacy_ns /
               batch_ns
            << '\n';

        std::cout
            << "checksum\t"
            << legacy_checksum
            << '\n';

        std::cout
            << "VERIFY_CHECKSUM\tYES\n";

        /*
         * Prevent optimization from deleting
         * benchmark work.
         */
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
