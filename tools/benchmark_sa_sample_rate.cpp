#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_strategy.hpp"
#include "primerpair/single_primer_search.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using Clock =
    std::chrono::steady_clock;


std::string make_reference(
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


char mutated_base(
    const char base
) {
    return
        base == 'A'
            ? 'T'
            : 'A';
}


std::uint64_t checksum(
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


struct Result {
    std::size_t sample_rate{0};

    std::size_t branching_count{0};

    double ns_per_primer{0.0};

    std::uint64_t sampled_sa_bytes{0};

    std::uint64_t sampled_sa_count{0};

    std::uint64_t checksum{0};
};


Result run_rate(
    const std::string& reference,
    const std::vector<std::string>& primers,
    const std::size_t sample_rate,
    const std::size_t repetitions
) {
    using namespace primerpair;

    constexpr std::size_t
        anchor_length = 12;


    PackedReference packed(
        reference
    );


    BidirectionalFMIndex bifm(
        reference,
        sample_rate
    );


    SinglePrimerSearchEngine engine(
        bifm,
        packed
    );


    std::vector<std::string_view>
        branching;

    branching.reserve(
        primers.size()
    );


    for (
        const auto& primer :
        primers
    ) {
        const auto decision =
            engine.router().decide(
                primer,
                anchor_length,
                3
            );


        if (
            decision.recommended_strategy ==
            SearchStrategy::DirectBranching
        ) {
            branching.push_back(
                primer
            );
        }
    }


    if (
        branching.empty()
    ) {
        throw std::runtime_error(
            "No branching primers."
        );
    }


    /*
     * Correctness checksum outside timing.
     */
    std::uint64_t expected_checksum = 0;

    for (
        const auto primer :
        branching
    ) {
        expected_checksum ^=
            checksum(
                engine.search(
                    primer,
                    anchor_length,
                    3
                )
            );
    }


    std::uint64_t sink = 0;


    const auto begin =
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


    const auto end =
        Clock::now();


    const double ns =
        std::chrono::duration<
            double,
            std::nano
        >(
            end - begin
        ).count();


    /*
     * BiFM contains forward + reverse FM index.
     */
    const std::uint64_t sampled_bytes =
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


    const std::uint64_t sampled_count =
        static_cast<std::uint64_t>(
            bifm
                .forward_index()
                .sampled_sa_count()
        )
        +
        static_cast<std::uint64_t>(
            bifm
                .reverse_index()
                .sampled_sa_count()
        );


    /*
     * Keep benchmark work alive.
     */
    if (
        sink ==
        std::numeric_limits<
            std::uint64_t
        >::max()
    ) {
        std::cerr
            << "sink\t"
            << sink
            << '\n';
    }


    return Result{
        sample_rate,
        branching.size(),
        ns
            /
            static_cast<double>(
                branching.size()
                *
                repetitions
            ),
        sampled_bytes,
        sampled_count,
        expected_checksum
    };
}

}  // namespace


int main(
    const int argc,
    char** argv
) {
    try {
        std::size_t repetitions = 50;

        if (
            argc >= 2
        ) {
            repetitions =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[1]
                    )
                );
        }


        std::size_t requested_rate = 0;

        if (
            argc >= 3
        ) {
            requested_rate =
                static_cast<std::size_t>(
                    std::stoull(
                        argv[2]
                    )
                );

            if (
                requested_rate != 4 &&
                requested_rate != 8 &&
                requested_rate != 16 &&
                requested_rate != 32
            ) {
                throw std::invalid_argument(
                    "Sample rate must be "
                    "4, 8, 16, or 32."
                );
            }
        }


        constexpr std::size_t
            primer_length = 24;


        const std::string reference =
            make_reference(
                200000
            );


        std::vector<std::string>
            primers;

        primers.reserve(
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


            primers.push_back(
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


                primers.push_back(
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
            primers.push_back(
                "TTTTTTTTTTTT"
                "ACGTACGTACGT"
            );

            primers.push_back(
                "AAAAAAAAAAAA"
                "ACACACACACAC"
            );

            primers.push_back(
                "CCCCCCCCCCCC"
                "GTGTGTGTGTGT"
            );
        }


        std::vector<std::size_t> rates;

        if (
            requested_rate == 0
        ) {
            rates = {
                4,
                8,
                16,
                32
            };

        } else {
            rates = {
                requested_rate
            };
        }


        std::uint64_t reference_checksum = 0;

        bool first = true;


        std::cout
            << "sample_rate"
            << '\t'
            << "branching_count"
            << '\t'
            << "ns_per_primer"
            << '\t'
            << "sampled_sa_bytes"
            << '\t'
            << "sampled_sa_count"
            << '\t'
            << "checksum"
            << '\n';


        for (
            const std::size_t rate :
            rates
        ) {
            const Result result =
                run_rate(
                    reference,
                    primers,
                    rate,
                    repetitions
                );


            if (first) {
                reference_checksum =
                    result.checksum;

                first = false;

            } else if (
                result.checksum !=
                reference_checksum
            ) {
                throw std::runtime_error(
                    "Checksum mismatch across "
                    "SA sample rates."
                );
            }


            std::cout
                << result.sample_rate
                << '\t'
                << result.branching_count
                << '\t'
                << result.ns_per_primer
                << '\t'
                << result.sampled_sa_bytes
                << '\t'
                << result.sampled_sa_count
                << '\t'
                << result.checksum
                << '\n';
        }


        std::cout
            << "VERIFY_CHECKSUM\tYES\n";


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
