#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/single_primer_search.hpp"
#include "primerpair/search_strategy.hpp"

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
checksum(
    const primerpair::SinglePrimerSearchResult& result
) {
    std::uint64_t value =
        result.hits.size();

    for (
        const auto& hit :
        result.hits
    ) {
        value ^=
            (
                hit.position *
                0x9E3779B97F4A7C15ULL
            )
            ^
            static_cast<std::uint64_t>(
                hit.mismatches + 1
            );
    }

    return value;
}

}  // namespace


int
main(
    int argc,
    char** argv
) {
    try {
        using namespace primerpair;

        std::size_t repetitions = 100;

        if (argc >= 2) {
            repetitions =
                static_cast<std::size_t>(
                    std::stoull(argv[1])
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

        SinglePrimerSearchEngine engine(
            bifm,
            packed
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
                    i * 3571 + 123
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


        std::vector<std::string_view>
            candidate;

        std::vector<std::string_view>
            branching;


        candidate.reserve(
            owned_primers.size()
        );

        branching.reserve(
            owned_primers.size()
        );


        for (
            const auto& primer :
            owned_primers
        ) {
            const auto decision =
                engine.router().decide(
                    primer,
                    anchor_length,
                    3
                );


            if (
                decision.recommended_strategy ==
                SearchStrategy::
                    AnchorCandidateVerification
            ) {
                candidate.push_back(
                    primer
                );

            } else {
                branching.push_back(
                    primer
                );
            }
        }


        std::uint64_t sink = 0;


        const auto candidate_begin =
            Clock::now();

        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            for (
                const auto primer :
                candidate
            ) {
                sink ^=
                    checksum(
                        engine.search(
                            primer,
                            anchor_length,
                            3
                        )
                    );
            }
        }

        const auto candidate_end =
            Clock::now();


        const auto branching_begin =
            Clock::now();

        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            for (
                const auto primer :
                branching
            ) {
                sink ^=
                    checksum(
                        engine.search(
                            primer,
                            anchor_length,
                            3
                        )
                    );
            }
        }

        const auto branching_end =
            Clock::now();


        const double candidate_total_ns =
            std::chrono::duration<
                double,
                std::nano
            >(
                candidate_end -
                candidate_begin
            ).count();


        const double branching_total_ns =
            std::chrono::duration<
                double,
                std::nano
            >(
                branching_end -
                branching_begin
            ).count();


        const double candidate_ns =
            candidate_total_ns
            /
            static_cast<double>(
                candidate.size()
                *
                repetitions
            );


        const double branching_ns =
            branching_total_ns
            /
            static_cast<double>(
                branching.size()
                *
                repetitions
            );


        const double total =
            candidate_total_ns
            +
            branching_total_ns;


        std::cout
            << "candidate_count\t"
            << candidate.size()
            << '\n';

        std::cout
            << "branching_count\t"
            << branching.size()
            << '\n';

        std::cout
            << "candidate_ns_per_primer\t"
            << candidate_ns
            << '\n';

        std::cout
            << "branching_ns_per_primer\t"
            << branching_ns
            << '\n';

        std::cout
            << "candidate_time_fraction\t"
            << candidate_total_ns / total
            << '\n';

        std::cout
            << "branching_time_fraction\t"
            << branching_total_ns / total
            << '\n';

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
