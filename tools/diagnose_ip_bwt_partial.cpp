#include "primerpair/fm_index.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/ip_bwt_rmi.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>

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

void print_interval(
    const char* name,
    const primerpair::Interval interval
) {
    std::cout
        << name
        << '\t'
        << interval.begin
        << '\t'
        << interval.end
        << '\t'
        << "size="
        << interval.size()
        << '\n';
}

}  // namespace


int main() {
    using namespace primerpair;

    const std::string reference =
        make_reference();

    const std::string query =
        reference.substr(
            0,
            18
        );


    FMIndex fm(
        reference
    );

    IPBWTIndex ip(
        reference,
        21
    );

    IPBWTRMI rmi(
        ip,
        64,
        8
    );


    const Interval fm_interval =
        fm.backward_search(
            query
        );

    const Interval ip_interval =
        ip.exact_search(
            query
        );

    const Interval rmi_interval =
        rmi.exact_search(
            query
        );


    std::cout
        << "reference_length\t"
        << reference.size()
        << '\n';

    std::cout
        << "query_length\t"
        << query.size()
        << '\n';

    std::cout
        << "K\t"
        << ip.chunk_length()
        << '\n';

    std::cout
        << "query\t"
        << query
        << '\n';


    print_interval(
        "FM",
        fm_interval
    );

    print_interval(
        "IPBWT_BINARY",
        ip_interval
    );

    print_interval(
        "IPBWT_RMI",
        rmi_interval
    );


    const auto fm_positions =
        fm.locate(
            fm_interval
        );

    const auto ip_positions =
        ip.locate(
            ip_interval
        );

    const auto rmi_positions =
        ip.locate(
            rmi_interval
        );


    std::cout
        << "FM_POSITIONS";

    for (
        const auto position :
        fm_positions
    ) {
        std::cout
            << '\t'
            << position;
    }

    std::cout
        << '\n';


    std::cout
        << "IPBWT_POSITIONS";

    for (
        const auto position :
        ip_positions
    ) {
        std::cout
            << '\t'
            << position;
    }

    std::cout
        << '\n';


    std::cout
        << "RMI_POSITIONS";

    for (
        const auto position :
        rmi_positions
    ) {
        std::cout
            << '\t'
            << position;
    }

    std::cout
        << '\n';


    /*
     * Show the exact artificial K=21 bounds used
     * by the partial-chunk branch.
     */
    std::string low_chunk =
        query;

    low_chunk.push_back(
        '$'
    );

    low_chunk.append(
        21 -
        query.size() -
        1,
        'A'
    );


    std::string high_chunk =
        query;

    high_chunk.append(
        21 -
        query.size(),
        'T'
    );


    std::cout
        << "LOW_CHUNK\t"
        << low_chunk
        << '\n';

    std::cout
        << "HIGH_CHUNK\t"
        << high_chunk
        << '\n';


    const auto binary_low =
        ip.lower_bound(
            low_chunk,
            0
        );

    const auto binary_high =
        ip.lower_bound(
            high_chunk,
            ip.row_count()
        );


    const auto rmi_low =
        rmi.lower_bound_diagnostics(
            low_chunk,
            0
        );

    const auto rmi_high =
        rmi.lower_bound_diagnostics(
            high_chunk,
            ip.row_count()
        );


    std::cout
        << "BINARY_LOW\t"
        << binary_low
        << '\n';

    std::cout
        << "BINARY_HIGH\t"
        << binary_high
        << '\n';

    std::cout
        << "RMI_LOW\t"
        << rmi_low.position
        << '\t'
        << "fallback="
        << (
            rmi_low.used_global_fallback
            ? "YES"
            : "NO"
        )
        << '\n';

    std::cout
        << "RMI_HIGH\t"
        << rmi_high.position
        << '\t'
        << "fallback="
        << (
            rmi_high.used_global_fallback
            ? "YES"
            : "NO"
        )
        << '\n';


    const bool correct =
        fm_interval.begin ==
            ip_interval.begin
        &&
        fm_interval.end ==
            ip_interval.end
        &&
        fm_interval.begin ==
            rmi_interval.begin
        &&
        fm_interval.end ==
            rmi_interval.end;


    std::cout
        << "ALL_CHECKS\t"
        << (
            correct
            ? "YES"
            : "NO"
        )
        << '\n';


    return
        correct
        ? 0
        : 1;
}
