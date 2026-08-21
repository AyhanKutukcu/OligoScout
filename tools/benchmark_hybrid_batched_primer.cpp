#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
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
checksum_hits(
    const std::vector<
        primerpair::PrimerSearchHit
    >& hits
) {
    std::uint64_t checksum =
        hits.size();

    for (
        const auto& hit :
        hits
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


        std::size_t sa_sample_rate =
            FMIndex::kDefaultSuffixArraySampleRate;

        if (argc >= 3) {
            sa_sample_rate =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[2]
                    )
                );

            if (
                sa_sample_rate != 4 &&
                sa_sample_rate != 8 &&
                sa_sample_rate != 16 &&
                sa_sample_rate != 32
            ) {
                throw std::invalid_argument(
                    "SA sample rate must be "
                    "4, 8, 16, or 32."
                );
            }
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
            reference,
            sa_sample_rate
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
            legacy_engine(
                bifm,
                packed
            );

        HybridBatchedPrimerSearchEngine
            hybrid_engine(
                candidate_engine,
                legacy_engine
            );


        std::vector<std::string>
            owned_primers;

        owned_primers.reserve(
            5000
        );


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


        /*
         * Add repeat-rich primers to ensure
         * DirectBranching remains represented.
         */
        for (
            std::size_t i = 0;
            i < 100;
            ++i
        ) {
            owned_primers.push_back(
                "TTTTTTTTTTTT"
                "ACGTACGTACGT"
            );

            owned_primers.push_back(
                "AAAAAAAAAAAA"
                "ACACACACACAC"
            );

            owned_primers.push_back(
                "CCCCCCCCCCCC"
                "GTGTGTGTGTGT"
            );
        }


        std::vector<
            HybridBatchedPrimerRequest
        > requests;

        requests.reserve(
            owned_primers.size()
        );

        for (
            const auto& primer :
            owned_primers
        ) {
            requests.push_back(
                HybridBatchedPrimerRequest{
                    primer,
                    3
                }
            );
        }


        /*
         * Correctness pass.
         */
        const auto hybrid_once =
            hybrid_engine.search(
                requests,
                anchor_length
            );


        if (
            hybrid_once.size() !=
            requests.size()
        ) {
            throw std::runtime_error(
                "Hybrid result-count mismatch."
            );
        }


        std::uint64_t legacy_checksum = 0;
        std::uint64_t hybrid_checksum = 0;

        std::size_t candidate_routes = 0;
        std::size_t branching_routes = 0;


        for (
            std::size_t i = 0;
            i < requests.size();
            ++i
        ) {
            const auto legacy =
                legacy_engine.search(
                    requests[i].primer,
                    anchor_length,
                    requests[i]
                        .max_mismatches
                );


            legacy_checksum ^=
                checksum_hits(
                    legacy.hits
                );


            hybrid_checksum ^=
                checksum_hits(
                    hybrid_once[i].hits
                );


            if (
                hybrid_once[i]
                    .used_candidate_backend
            ) {
                ++candidate_routes;
            } else {
                ++branching_routes;
            }
        }


        if (
            legacy_checksum !=
            hybrid_checksum
        ) {
            throw std::runtime_error(
                "Hybrid benchmark checksum mismatch."
            );
        }


        std::uint64_t sink = 0;


        /*
         * Legacy full mixed-panel path.
         */
        const auto legacy_begin =
            Clock::now();

        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            for (
                const auto& request :
                requests
            ) {
                const auto result =
                    legacy_engine.search(
                        request.primer,
                        anchor_length,
                        request.max_mismatches
                    );

                sink ^=
                    checksum_hits(
                        result.hits
                    );
            }
        }

        const auto legacy_end =
            Clock::now();


        /*
         * Hybrid full mixed-panel path.
         */
        const auto hybrid_begin =
            Clock::now();

        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            const auto results =
                hybrid_engine.search(
                    requests,
                    anchor_length
                );

            for (
                const auto& result :
                results
            ) {
                sink ^=
                    checksum_hits(
                        result.hits
                    );
            }
        }

        const auto hybrid_end =
            Clock::now();


        const double query_count =
            static_cast<double>(
                requests.size()
                *
                repetitions
            );


        const double legacy_ns =
            std::chrono::duration<
                double,
                std::nano
            >(
                legacy_end -
                legacy_begin
            ).count()
            /
            query_count;


        const double hybrid_ns =
            std::chrono::duration<
                double,
                std::nano
            >(
                hybrid_end -
                hybrid_begin
            ).count()
            /
            query_count;


        std::cout
            << "sa_sample_rate\t"
            << sa_sample_rate
            << '\n';


        const std::uint64_t
            sampled_sa_bytes =
                static_cast<std::uint64_t>(
                    bifm
                        .forward_index()
                        .sampled_sa_memory_bytes()
                )
                +
                static_cast<std::uint64_t>(
                    bifm
                        .reverse_index()
                        .sampled_sa_memory_bytes()
                );


        std::cout
            << "sampled_sa_bytes\t"
            << sampled_sa_bytes
            << '\n';


        std::cout
            << "panel_size\t"
            << requests.size()
            << '\n';

        std::cout
            << "candidate_routes\t"
            << candidate_routes
            << '\n';

        std::cout
            << "branching_routes\t"
            << branching_routes
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
            << "hybrid_ns_per_primer\t"
            << hybrid_ns
            << '\n';

        std::cout
            << "speedup\t"
            << legacy_ns /
               hybrid_ns
            << '\n';

        std::cout
            << "checksum\t"
            << legacy_checksum
            << '\n';

        std::cout
            << "VERIFY_CHECKSUM\tYES\n";

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
