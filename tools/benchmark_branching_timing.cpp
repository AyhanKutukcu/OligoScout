#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_strategy.hpp"
#include "primerpair/single_primer_search.hpp"

#include <array>
#include <chrono>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using Clock =
    std::chrono::steady_clock;

using primerpair::BidirectionalInterval;


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


std::string
normalize_primer(
    const std::string_view primer
) {
    std::string result;
    result.reserve(primer.size());

    for (const char raw : primer) {
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
                result.push_back(base);
                break;

            default:
                throw std::invalid_argument(
                    "Invalid primer base."
                );
        }
    }

    return result;
}


struct BranchState {
    BidirectionalInterval state{};
    std::size_t mismatches{0};
};


struct Timing {
    double anchor_ns{0.0};

    /*
     * Production branching artık iki extension yolu
     * kullanıyor:
     *
     * 1. mismatch budget dolu:
     *      yalnız expected base -> scalar extend_left()
     *
     * 2. mismatch budget açık:
     *      A/C/G/T -> extend_left_all()
     */
    double scalar_extension_ns{0.0};
    double batched_extension_ns{0.0};

    double locate_ns{0.0};
    double bookkeeping_ns{0.0};

    std::uint64_t
        scalar_expected_extensions{0};

    std::uint64_t
        scalar_empty_extensions{0};

    std::uint64_t
        batched_parent_extensions{0};

    std::uint64_t
        batched_children_generated{0};

    std::uint64_t
        batched_empty_children{0};

    std::uint64_t locate_calls{0};
    std::uint64_t located_positions{0};

    std::uint64_t final_states{0};
};


Timing
run_profile(
    const primerpair::BidirectionalFMIndex& index,
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) {
    Timing timing;

    const auto function_begin =
        Clock::now();


    const std::string normalized =
        normalize_primer(
            primer
        );


    const std::size_t prefix_length =
        normalized.size()
        -
        anchor_length;


    const std::string_view anchor(
        normalized.data()
            +
            prefix_length,
        anchor_length
    );


    const auto anchor_begin =
        Clock::now();


    const BidirectionalInterval
        anchor_state =
            index.search(
                anchor
            );


    const auto anchor_end =
        Clock::now();


    timing.anchor_ns +=
        std::chrono::duration<
            double,
            std::nano
        >(
            anchor_end -
            anchor_begin
        ).count();


    if (
        anchor_state.empty()
    ) {
        timing.bookkeeping_ns =
            std::chrono::duration<
                double,
                std::nano
            >(
                Clock::now() -
                function_begin
            ).count()
            -
            timing.anchor_ns;

        return timing;
    }


    std::vector<BranchState> current;

    current.push_back(
        BranchState{
            anchor_state,
            0
        }
    );


    constexpr std::array<char, 4>
        alphabet{
            'A',
            'C',
            'G',
            'T'
        };


    for (
        std::size_t remaining =
            prefix_length;
        remaining > 0;
        --remaining
    ) {
        const std::size_t primer_position =
            remaining - 1;


        const char expected =
            normalized.at(
                primer_position
            );


        std::vector<BranchState> next;

        next.reserve(
            current.size() * 4
        );


        for (
            const auto& branch :
            current
        ) {
            /*
             * Production SinglePrimerSearchEngine ile
             * aynı karar:
             *
             * mismatch budget tamamen dolmuşsa yalnız
             * primerde beklenen baz legal olabilir.
             */
            if (
                branch.mismatches ==
                max_mismatches
            ) {
                const auto extend_begin =
                    Clock::now();

                const BidirectionalInterval
                    next_state =
                        index.extend_left(
                            branch.state,
                            expected
                        );

                const auto extend_end =
                    Clock::now();

                timing.scalar_extension_ns +=
                    std::chrono::duration<
                        double,
                        std::nano
                    >(
                        extend_end -
                        extend_begin
                    ).count();

                ++timing
                    .scalar_expected_extensions;

                if (
                    next_state.empty()
                ) {
                    ++timing
                        .scalar_empty_extensions;

                    continue;
                }

                next.push_back(
                    BranchState{
                        next_state,
                        branch.mismatches
                    }
                );

                continue;
            }

            /*
             * mismatch budget hâlâ açıksa A/C/G/T
             * çocuklarının tamamı potansiyel olarak
             * legal.
             *
             * Parent başına tek extend_left_all().
             */
            const auto extend_begin =
                Clock::now();

            const auto children =
                index.extend_left_all(
                    branch.state
                );

            const auto extend_end =
                Clock::now();

            timing.batched_extension_ns +=
                std::chrono::duration<
                    double,
                    std::nano
                >(
                    extend_end -
                    extend_begin
                ).count();

            ++timing
                .batched_parent_extensions;

            timing
                .batched_children_generated +=
                    children.size();

            for (
                std::size_t base_index = 0;
                base_index < alphabet.size();
                ++base_index
            ) {
                const char base =
                    alphabet.at(
                        base_index
                    );

                const std::size_t mismatches =
                    branch.mismatches +
                    (
                        base == expected
                            ? 0
                            : 1
                    );

                if (
                    mismatches >
                    max_mismatches
                ) {
                    throw std::logic_error(
                        "Unexpected mismatch-budget "
                        "violation in batched branch."
                    );
                }

                const BidirectionalInterval&
                    next_state =
                        children.at(
                            base_index
                        );

                if (
                    next_state.empty()
                ) {
                    ++timing
                        .batched_empty_children;

                    continue;
                }

                next.push_back(
                    BranchState{
                        next_state,
                        mismatches
                    }
                );
            }
        }


        current =
            std::move(
                next
            );


        if (
            current.empty()
        ) {
            break;
        }
    }


    timing.final_states +=
        current.size();


    for (
        const auto& branch :
        current
    ) {
        const auto locate_begin =
            Clock::now();


        const auto positions =
            index.locate_unsorted(
                branch.state
            );


        const auto locate_end =
            Clock::now();


        timing.locate_ns +=
            std::chrono::duration<
                double,
                std::nano
            >(
                locate_end -
                locate_begin
            ).count();


        ++timing.locate_calls;

        timing.located_positions +=
            positions.size();
    }


    const double total_ns =
        std::chrono::duration<
            double,
            std::nano
        >(
            Clock::now() -
            function_begin
        ).count();


    timing.bookkeeping_ns =
        total_ns
        -
        timing.anchor_ns
        -
        timing.scalar_extension_ns
        -
        timing.batched_extension_ns
        -
        timing.locate_ns;


    return timing;
}


void
add(
    Timing& total,
    const Timing& value
) {
    total.anchor_ns +=
        value.anchor_ns;

    total.scalar_extension_ns +=
        value.scalar_extension_ns;

    total.batched_extension_ns +=
        value.batched_extension_ns;

    total.locate_ns +=
        value.locate_ns;

    total.bookkeeping_ns +=
        value.bookkeeping_ns;

    total.scalar_expected_extensions +=
        value.scalar_expected_extensions;

    total.scalar_empty_extensions +=
        value.scalar_empty_extensions;

    total.batched_parent_extensions +=
        value.batched_parent_extensions;

    total.batched_children_generated +=
        value.batched_children_generated;

    total.batched_empty_children +=
        value.batched_empty_children;

    total.locate_calls +=
        value.locate_calls;

    total.located_positions +=
        value.located_positions;

    total.final_states +=
        value.final_states;
}


}  // namespace


int
main(
    const int argc,
    char** argv
) {
    try {
        using namespace primerpair;


        std::size_t repetitions = 20;

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


        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            primer_length = 24;

        constexpr std::size_t
            max_mismatches = 3;


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
            branching;


        for (
            const auto& primer :
            owned_primers
        ) {
            const auto decision =
                engine.router().decide(
                    primer,
                    anchor_length,
                    max_mismatches
                );


            if (
                decision.recommended_strategy ==
                SearchStrategy::
                    DirectBranching
            ) {
                branching.push_back(
                    primer
                );
            }
        }


        Timing total;


        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            for (
                const auto primer :
                branching
            ) {
                add(
                    total,
                    run_profile(
                        bifm,
                        primer,
                        anchor_length,
                        max_mismatches
                    )
                );
            }
        }


        const double measured_total =
            total.anchor_ns
            +
            total.scalar_extension_ns
            +
            total.batched_extension_ns
            +
            total.locate_ns
            +
            total.bookkeeping_ns;


        const double primer_runs =
            static_cast<double>(
                branching.size()
                *
                repetitions
            );


        std::cout
            << "branching_primers\t"
            << branching.size()
            << '\n';


        std::cout
            << "repetitions\t"
            << repetitions
            << '\n';


        const double
            scalar_extension_total =
                total.scalar_extension_ns;

        const double
            batched_extension_total =
                total.batched_extension_ns;

        const double
            extension_total =
                scalar_extension_total +
                batched_extension_total;


        const double scalar_calls =
            static_cast<double>(
                total.scalar_expected_extensions
            );

        const double batch_calls =
            static_cast<double>(
                total.batched_parent_extensions
            );


        std::cout
            << "anchor_ns_per_primer\t"
            << total.anchor_ns /
               primer_runs
            << '\n';


        std::cout
            << "scalar_extension_ns_per_primer\t"
            << scalar_extension_total /
               primer_runs
            << '\n';


        std::cout
            << "batched_extension_ns_per_primer\t"
            << batched_extension_total /
               primer_runs
            << '\n';


        std::cout
            << "extension_ns_per_primer\t"
            << extension_total /
               primer_runs
            << '\n';


        std::cout
            << "locate_ns_per_primer\t"
            << total.locate_ns /
               primer_runs
            << '\n';


        std::cout
            << "bookkeeping_ns_per_primer\t"
            << total.bookkeeping_ns /
               primer_runs
            << '\n';


        std::cout
            << "anchor_time_fraction\t"
            << total.anchor_ns /
               measured_total
            << '\n';


        std::cout
            << "scalar_extension_time_fraction\t"
            << scalar_extension_total /
               measured_total
            << '\n';


        std::cout
            << "batched_extension_time_fraction\t"
            << batched_extension_total /
               measured_total
            << '\n';


        std::cout
            << "extension_time_fraction\t"
            << extension_total /
               measured_total
            << '\n';


        std::cout
            << "locate_time_fraction\t"
            << total.locate_ns /
               measured_total
            << '\n';


        std::cout
            << "bookkeeping_time_fraction\t"
            << total.bookkeeping_ns /
               measured_total
            << '\n';


        std::cout
            << "scalar_expected_extensions_per_primer\t"
            << scalar_calls /
               primer_runs
            << '\n';


        std::cout
            << "scalar_empty_extensions_per_primer\t"
            << static_cast<double>(
                   total.scalar_empty_extensions
               ) /
               primer_runs
            << '\n';


        std::cout
            << "batched_parent_extensions_per_primer\t"
            << batch_calls /
               primer_runs
            << '\n';


        std::cout
            << "batched_children_generated_per_primer\t"
            << static_cast<double>(
                   total.batched_children_generated
               ) /
               primer_runs
            << '\n';


        std::cout
            << "batched_empty_children_per_primer\t"
            << static_cast<double>(
                   total.batched_empty_children
               ) /
               primer_runs
            << '\n';


        std::cout
            << "index_extension_calls_per_primer\t"
            << (
                   scalar_calls +
                   batch_calls
               ) /
               primer_runs
            << '\n';


        std::cout
            << "logical_children_per_primer\t"
            << (
                   scalar_calls +
                   static_cast<double>(
                       total.batched_children_generated
                   )
               ) /
               primer_runs
            << '\n';


        std::cout
            << "scalar_ns_per_call\t"
            << (
                   scalar_calls > 0.0
                       ? scalar_extension_total /
                         scalar_calls
                       : 0.0
               )
            << '\n';


        std::cout
            << "batched_ns_per_parent_call\t"
            << (
                   batch_calls > 0.0
                       ? batched_extension_total /
                         batch_calls
                       : 0.0
               )
            << '\n';


        std::cout
            << "locate_calls_per_primer\t"
            << static_cast<double>(
                   total.locate_calls
               ) /
               primer_runs
            << '\n';


        std::cout
            << "located_positions_per_primer\t"
            << static_cast<double>(
                   total.located_positions
               ) /
               primer_runs
            << '\n';


        std::cout
            << "final_states_per_primer\t"
            << static_cast<double>(
                   total.final_states
               ) /
               primer_runs
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
