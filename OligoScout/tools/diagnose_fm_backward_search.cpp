#include "primerpair/fm_index.hpp"
#include "primerpair/suffix_array_builder.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

std::string make_reference() {
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
        5000
    );

    for (
        std::size_t i = 0;
        i < 4096;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL
            +
            1442695040888963407ULL;

        const std::size_t index =
            static_cast<std::size_t>(
                (
                    state >>
                    32
                )
                &
                3ULL
            );

        reference.push_back(
            alphabet[index]
        );
    }

    const std::string motif =
        "ACGTACGTACGTACGTACGTACGT"
        "TTTTAAAACCCCGGGG"
        "ACGTACGTACGTACGTACGTACGT";

    for (
        std::size_t repeat = 0;
        repeat < 12;
        ++repeat
    ) {
        reference += motif;
    }

    return reference;
}


bool suffix_starts_with(
    const std::string& text,
    const std::uint32_t start,
    const std::string& pattern
) {
    const std::size_t position =
        static_cast<std::size_t>(
            start
        );

    if (
        position + pattern.size() >
        text.size()
    ) {
        return false;
    }

    for (
        std::size_t i = 0;
        i < pattern.size();
        ++i
    ) {
        if (
            text.at(
                position + i
            )
            !=
            pattern.at(i)
        ) {
            return false;
        }
    }

    return true;
}


primerpair::Interval naive_interval(
    const std::string& text,
    const std::vector<std::uint32_t>& sa,
    const std::string& pattern
) {
    bool found = false;

    std::uint64_t begin = 0;
    std::uint64_t end = 0;

    for (
        std::size_t row = 0;
        row < sa.size();
        ++row
    ) {
        if (
            suffix_starts_with(
                text,
                sa.at(row),
                pattern
            )
        ) {
            if (!found) {
                begin =
                    static_cast<std::uint64_t>(
                        row
                    );

                found = true;
            }

            end =
                static_cast<std::uint64_t>(
                    row + 1
                );
        }
    }

    if (!found) {
        throw std::runtime_error(
            "Naive interval unexpectedly empty "
            "for a suffix of a known-present query."
        );
    }

    return {
        begin,
        end
    };
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string reference =
            make_reference();

        const std::string text =
            reference + "$";

        const std::string query =
            reference.substr(
                0,
                18
            );

        FMIndex fm(
            reference
        );

        const auto sa =
            build_suffix_array_prefix_doubling(
                text
            );


        std::cout
            << "reference_length\t"
            << reference.size()
            << '\n';

        std::cout
            << "query\t"
            << query
            << '\n';

        std::cout
            << "known_position\t0\n";


        Interval fm_interval{
            0,
            static_cast<std::uint64_t>(
                text.size()
            )
        };

        std::string current_pattern;


        for (
            std::size_t offset = 0;
            offset < query.size();
            ++offset
        ) {
            const std::size_t query_position =
                query.size() -
                1 -
                offset;

            const char base =
                query.at(
                    query_position
                );

            current_pattern.insert(
                current_pattern.begin(),
                base
            );


            fm_interval =
                fm.backward_extend(
                    fm_interval,
                    base
                );


            const Interval expected =
                naive_interval(
                    text,
                    sa,
                    current_pattern
                );


            const bool same =
                fm_interval.begin ==
                    expected.begin
                &&
                fm_interval.end ==
                    expected.end;


            std::cout
                << "STEP\t"
                << offset + 1
                << '\t'
                << "base="
                << base
                << '\t'
                << "pattern="
                << current_pattern
                << '\t'
                << "FM="
                << fm_interval.begin
                << ','
                << fm_interval.end
                << '\t'
                << "NAIVE="
                << expected.begin
                << ','
                << expected.end
                << '\t'
                << (
                    same
                    ? "OK"
                    : "MISMATCH"
                )
                << '\n';


            if (!same) {
                std::cout
                    << "FIRST_DIVERGENCE_STEP\t"
                    << offset + 1
                    << '\n';

                std::cout
                    << "FIRST_DIVERGENCE_PATTERN\t"
                    << current_pattern
                    << '\n';

                std::cout
                    << "ALL_CHECKS\tNO\n";

                return 1;
            }
        }


        const auto final_direct =
            fm.backward_search(
                query
            );


        std::cout
            << "direct_backward_search\t"
            << final_direct.begin
            << '\t'
            << final_direct.end
            << '\n';


        std::cout
            << "ALL_CHECKS\tYES\n";

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
