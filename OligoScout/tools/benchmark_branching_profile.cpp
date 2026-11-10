#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_strategy.hpp"
#include "primerpair/single_primer_search.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

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


/*
 * Full state identity:
 *
 * Same BiFM interval
 * + same matched length
 * + same mismatch count.
 */
struct StateKey {
    std::uint64_t f_begin{0};
    std::uint64_t f_end{0};

    std::uint64_t r_begin{0};
    std::uint64_t r_end{0};

    std::size_t length{0};
    std::size_t mismatches{0};

    bool operator==(
        const StateKey&
    ) const = default;
};


/*
 * Interval identity ignoring mismatch count.
 *
 * Used only to measure whether identical genomic
 * search states arise with different mismatch
 * costs.
 */
struct IntervalKey {
    std::uint64_t f_begin{0};
    std::uint64_t f_end{0};

    std::uint64_t r_begin{0};
    std::uint64_t r_end{0};

    std::size_t length{0};

    bool operator==(
        const IntervalKey&
    ) const = default;
};


inline void
hash_combine(
    std::size_t& seed,
    const std::uint64_t value
) noexcept {
    seed ^=
        std::hash<std::uint64_t>{}(value)
        +
        0x9E3779B97F4A7C15ULL
        +
        (seed << 6U)
        +
        (seed >> 2U);
}


struct StateKeyHash {
    std::size_t operator()(
        const StateKey& key
    ) const noexcept {
        std::size_t seed = 0;

        hash_combine(seed, key.f_begin);
        hash_combine(seed, key.f_end);
        hash_combine(seed, key.r_begin);
        hash_combine(seed, key.r_end);
        hash_combine(
            seed,
            static_cast<std::uint64_t>(
                key.length
            )
        );
        hash_combine(
            seed,
            static_cast<std::uint64_t>(
                key.mismatches
            )
        );

        return seed;
    }
};


struct IntervalKeyHash {
    std::size_t operator()(
        const IntervalKey& key
    ) const noexcept {
        std::size_t seed = 0;

        hash_combine(seed, key.f_begin);
        hash_combine(seed, key.f_end);
        hash_combine(seed, key.r_begin);
        hash_combine(seed, key.r_end);
        hash_combine(
            seed,
            static_cast<std::uint64_t>(
                key.length
            )
        );

        return seed;
    }
};


StateKey
make_state_key(
    const BranchState& branch
) {
    return StateKey{
        branch.state.forward.begin,
        branch.state.forward.end,
        branch.state.reverse.begin,
        branch.state.reverse.end,
        branch.state.length,
        branch.mismatches
    };
}


IntervalKey
make_interval_key(
    const BidirectionalInterval& state
) {
    return IntervalKey{
        state.forward.begin,
        state.forward.end,
        state.reverse.begin,
        state.reverse.end,
        state.length
    };
}


struct Profile {
    std::uint64_t alphabet_candidates{0};

    std::uint64_t budget_pruned{0};

    std::uint64_t extension_attempts{0};

    std::uint64_t empty_extensions{0};

    std::uint64_t states_generated{0};

    std::uint64_t duplicate_states{0};

    std::uint64_t interval_collisions{0};

    std::uint64_t final_states{0};

    std::uint64_t final_interval_rows{0};

    std::uint64_t locate_calls{0};

    std::uint64_t located_positions{0};

    std::uint64_t max_live_states{0};

    std::uint64_t max_next_states{0};
};


Profile
profile_branching(
    const primerpair::BidirectionalFMIndex& index,
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) {
    const std::string normalized =
        normalize_primer(
            primer
        );

    if (
        anchor_length == 0 ||
        anchor_length >
            normalized.size()
    ) {
        throw std::invalid_argument(
            "Invalid anchor length."
        );
    }


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


    const BidirectionalInterval
        anchor_state =
            index.search(
                anchor
            );


    Profile profile;


    if (anchor_state.empty()) {
        return profile;
    }


    std::vector<BranchState> current;

    current.push_back(
        BranchState{
            anchor_state,
            0
        }
    );


    profile.max_live_states = 1;


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
            for (
                const char base :
                alphabet
            ) {
                ++profile.alphabet_candidates;


                const std::size_t mismatches =
                    branch.mismatches
                    +
                    (
                        base == expected
                            ? 0
                            : 1
                    );


                if (
                    mismatches >
                    max_mismatches
                ) {
                    ++profile.budget_pruned;
                    continue;
                }


                ++profile.extension_attempts;


                const BidirectionalInterval
                    next_state =
                        index.extend_left(
                            branch.state,
                            base
                        );


                if (
                    next_state.empty()
                ) {
                    ++profile.empty_extensions;
                    continue;
                }


                next.push_back(
                    BranchState{
                        next_state,
                        mismatches
                    }
                );


                ++profile.states_generated;
            }
        }


        /*
         * IMPORTANT:
         *
         * We measure duplicates but DO NOT remove
         * them. Production behavior remains
         * faithfully represented.
         */
        std::unordered_set<
            StateKey,
            StateKeyHash
        > full_states;


        std::unordered_map<
            IntervalKey,
            std::size_t,
            IntervalKeyHash
        > best_mismatch;


        full_states.reserve(
            next.size() * 2 + 1
        );

        best_mismatch.reserve(
            next.size() * 2 + 1
        );


        for (
            const auto& branch :
            next
        ) {
            const StateKey full_key =
                make_state_key(
                    branch
                );


            if (
                !full_states
                    .insert(
                        full_key
                    )
                    .second
            ) {
                ++profile.duplicate_states;
            }


            const IntervalKey interval_key =
                make_interval_key(
                    branch.state
                );


            const auto it =
                best_mismatch.find(
                    interval_key
                );


            if (
                it ==
                best_mismatch.end()
            ) {
                best_mismatch.emplace(
                    interval_key,
                    branch.mismatches
                );

            } else if (
                it->second !=
                branch.mismatches
            ) {
                ++profile.interval_collisions;

                if (
                    branch.mismatches <
                    it->second
                ) {
                    it->second =
                        branch.mismatches;
                }
            }
        }


        profile.max_next_states =
            std::max<std::uint64_t>(
                profile.max_next_states,
                next.size()
            );


        current =
            std::move(
                next
            );


        profile.max_live_states =
            std::max<std::uint64_t>(
                profile.max_live_states,
                current.size()
            );


        if (
            current.empty()
        ) {
            break;
        }
    }


    profile.final_states =
        current.size();


    for (
        const auto& branch :
        current
    ) {
        profile.final_interval_rows +=
            branch.state.size();


        ++profile.locate_calls;


        const auto positions =
            index.locate(
                branch.state
            );


        profile.located_positions +=
            positions.size();
    }


    return profile;
}


void
add_profile(
    Profile& total,
    const Profile& value
) {
    total.alphabet_candidates +=
        value.alphabet_candidates;

    total.budget_pruned +=
        value.budget_pruned;

    total.extension_attempts +=
        value.extension_attempts;

    total.empty_extensions +=
        value.empty_extensions;

    total.states_generated +=
        value.states_generated;

    total.duplicate_states +=
        value.duplicate_states;

    total.interval_collisions +=
        value.interval_collisions;

    total.final_states +=
        value.final_states;

    total.final_interval_rows +=
        value.final_interval_rows;

    total.locate_calls +=
        value.locate_calls;

    total.located_positions +=
        value.located_positions;

    total.max_live_states =
        std::max(
            total.max_live_states,
            value.max_live_states
        );

    total.max_next_states =
        std::max(
            total.max_next_states,
            value.max_next_states
        );
}

}  // namespace


int
main() {
    try {
        using namespace primerpair;


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


        /*
         * Same mixed-panel workload used in the
         * hybrid benchmark.
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


        Profile total;


        std::uint64_t worst_states = 0;

        std::string_view
            worst_primer;


        for (
            const auto primer :
            branching
        ) {
            const Profile profile =
                profile_branching(
                    bifm,
                    primer,
                    anchor_length,
                    max_mismatches
                );


            add_profile(
                total,
                profile
            );


            if (
                profile.max_live_states >
                worst_states
            ) {
                worst_states =
                    profile.max_live_states;

                worst_primer =
                    primer;
            }
        }


        const double empty_rate =
            total.extension_attempts == 0
            ?
            0.0
            :
            static_cast<double>(
                total.empty_extensions
            )
            /
            static_cast<double>(
                total.extension_attempts
            );


        const double duplicate_rate =
            total.states_generated == 0
            ?
            0.0
            :
            static_cast<double>(
                total.duplicate_states
            )
            /
            static_cast<double>(
                total.states_generated
            );


        const double budget_prune_rate =
            total.alphabet_candidates == 0
            ?
            0.0
            :
            static_cast<double>(
                total.budget_pruned
            )
            /
            static_cast<double>(
                total.alphabet_candidates
            );


        const double interval_collision_rate =
            total.states_generated == 0
            ?
            0.0
            :
            static_cast<double>(
                total.interval_collisions
            )
            /
            static_cast<double>(
                total.states_generated
            );


        std::cout
            << "branching_primers\t"
            << branching.size()
            << '\n';


        std::cout
            << "alphabet_candidates\t"
            << total.alphabet_candidates
            << '\n';


        std::cout
            << "budget_pruned\t"
            << total.budget_pruned
            << '\n';


        std::cout
            << "budget_prune_rate\t"
            << budget_prune_rate
            << '\n';


        std::cout
            << "extension_attempts\t"
            << total.extension_attempts
            << '\n';


        std::cout
            << "empty_extensions\t"
            << total.empty_extensions
            << '\n';


        std::cout
            << "empty_extension_rate\t"
            << empty_rate
            << '\n';


        std::cout
            << "states_generated\t"
            << total.states_generated
            << '\n';


        std::cout
            << "duplicate_states\t"
            << total.duplicate_states
            << '\n';


        std::cout
            << "duplicate_state_rate\t"
            << duplicate_rate
            << '\n';


        std::cout
            << "interval_collisions\t"
            << total.interval_collisions
            << '\n';


        std::cout
            << "interval_collision_rate\t"
            << interval_collision_rate
            << '\n';


        std::cout
            << "max_live_states\t"
            << total.max_live_states
            << '\n';


        std::cout
            << "max_next_states\t"
            << total.max_next_states
            << '\n';


        std::cout
            << "final_states\t"
            << total.final_states
            << '\n';


        std::cout
            << "final_interval_rows\t"
            << total.final_interval_rows
            << '\n';


        std::cout
            << "locate_calls\t"
            << total.locate_calls
            << '\n';


        std::cout
            << "located_positions\t"
            << total.located_positions
            << '\n';


        std::cout
            << "worst_primer\t"
            << worst_primer
            << '\n';


        std::cout
            << "worst_max_live_states\t"
            << worst_states
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
