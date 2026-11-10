#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "primerpair/fm_index.hpp"

namespace {

std::string generate_reference(
    const std::size_t length
) {
    static constexpr char alphabet[] = {
        'A', 'C', 'G', 'T'
    };

    std::string sequence;
    sequence.resize(length);

    /*
     * Deterministik xorshift64.
     * Her benchmark çalışmasında aynı referans üretilir.
     */
    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;

    for (std::size_t i = 0;
         i < length;
         ++i) {

        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;

        sequence.at(i) =
            alphabet[state & 3ULL];
    }

    return sequence;
}

}  // namespace

int main(
    int argc,
    char* argv[]
) {
    try {
        if (argc < 2 || argc > 3) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <reference_length_bp>"
                << " [sa_sample_rate]\n";

            return 2;
        }

        const std::size_t reference_length =
            static_cast<std::size_t>(
                std::stoull(argv[1])
            );

        const std::size_t sa_sample_rate =
            argc >= 3
                ? static_cast<std::size_t>(
                      std::stoull(argv[2])
                  )
                : primerpair::FMIndex::kDefaultSuffixArraySampleRate;

        if (reference_length < 100) {
            throw std::invalid_argument(
                "Reference length must be at least 100 bp."
            );
        }

        if (sa_sample_rate == 0) {
            throw std::invalid_argument(
                "SA sample rate must be greater than zero."
            );
        }

        std::string reference =
            generate_reference(
                reference_length
            );

        /*
         * Ortadaki gerçek bir 20-mer'i daha sonra
         * arayarak indeksin çalıştığını doğrularız.
         */
        const std::size_t pattern_start =
            reference_length / 2;

        const std::string pattern =
            reference.substr(
                pattern_start,
                20
            );

        /*
         * Timer yalnızca FM-index oluşturmayı ölçer.
         * Referans üretimi bu sürenin dışında.
         */
        const auto build_start =
            std::chrono::steady_clock::now();

        primerpair::FMIndex index(
            std::move(reference),
            sa_sample_rate
        );

        const auto build_stop =
            std::chrono::steady_clock::now();

        const auto build_ns =
            std::chrono::duration_cast<
                std::chrono::nanoseconds
            >(
                build_stop - build_start
            ).count();

        /*
         * İndeks doğruluk sanity-check.
         */
        const auto interval =
            index.backward_search(
                pattern
            );

        const auto positions =
            index.locate(
                interval
            );

        bool expected_position_found = false;

        for (const auto position : positions) {
            if (position == pattern_start) {
                expected_position_found = true;
                break;
            }
        }

        const std::size_t text_bytes =
            index.text_memory_bytes();

        const std::size_t bwt_bytes =
            index.bwt_memory_bytes();

        const std::size_t packed_bwt_bytes =
            index.packed_bwt_memory_bytes();

        const std::size_t rank_bytes =
            index.rank_memory_bytes();

        const std::size_t sampled_sa_bytes =
            index.sampled_sa_memory_bytes();

        const std::size_t sampled_sa_marker_bytes =
            index.sampled_sa_marker_memory_bytes();

        const std::size_t sampled_sa_prefix_bytes =
            index.sampled_sa_prefix_memory_bytes();

        const std::size_t sampled_sa_values_bytes =
            index.sampled_sa_values_memory_bytes();

        const std::size_t approximate_persistent_bytes =
            text_bytes +
            bwt_bytes +
            packed_bwt_bytes +
            rank_bytes +
            sampled_sa_bytes;

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "reference_length_bp\t"
            << reference_length
            << '\n';

        std::cout
            << "build_time_ms\t"
            << static_cast<double>(
                   build_ns
               ) /
               1'000'000.0
            << '\n';

        std::cout
            << "sa_sample_rate\t"
            << index.suffix_array_sample_rate()
            << '\n';

        std::cout
            << "sampled_sa_count\t"
            << index.sampled_sa_count()
            << '\n';

        std::cout
            << "text_bytes\t"
            << text_bytes
            << '\n';

        std::cout
            << "bwt_bytes\t"
            << bwt_bytes
            << '\n';

        std::cout
            << "packed_bwt_bytes\t"
            << packed_bwt_bytes
            << '\n';

        std::cout
            << "rank_bytes\t"
            << rank_bytes
            << '\n';

        std::cout
            << "sampled_sa_bytes\t"
            << sampled_sa_bytes
            << '\n';

        std::cout
            << "sampled_sa_marker_bytes\t"
            << sampled_sa_marker_bytes
            << '\n';

        std::cout
            << "sampled_sa_prefix_bytes\t"
            << sampled_sa_prefix_bytes
            << '\n';

        std::cout
            << "sampled_sa_values_bytes\t"
            << sampled_sa_values_bytes
            << '\n';

        std::cout
            << "approx_persistent_bytes\t"
            << approximate_persistent_bytes
            << '\n';

        std::cout
            << "test_pattern\t"
            << pattern
            << '\n';

        std::cout
            << "match_count\t"
            << interval.size()
            << '\n';

        std::cout
            << "expected_position\t"
            << pattern_start
            << '\n';

        std::cout
            << "expected_position_found\t"
            << (
                expected_position_found
                    ? "YES"
                    : "NO"
            )
            << '\n';

        if (!expected_position_found) {
            std::cerr
                << "ERROR: expected genomic position "
                << "was not recovered.\n";

            return 1;
        }

        return 0;

    } catch (const std::exception& exception) {
        std::cerr
            << "Fatal error: "
            << exception.what()
            << '\n';

        return 1;
    }
}
