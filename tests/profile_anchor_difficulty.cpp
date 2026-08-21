#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/bidirectional_fm_index.hpp"

namespace {

std::string load_fasta(
    const std::string& path
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open FASTA."
        );
    }

    std::string sequence;
    std::string line;

    while (std::getline(input, line)) {

        if (
            line.empty() ||
            line.front() == '>'
        ) {
            continue;
        }

        for (const char raw : line) {

            if (
                std::isspace(
                    static_cast<unsigned char>(raw)
                )
            ) {
                continue;
            }

            const char base =
                static_cast<char>(
                    std::toupper(
                        static_cast<unsigned char>(raw)
                    )
                );

            switch (base) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':
                    sequence.push_back(base);
                    break;

                case 'R':
                case 'Y':
                case 'S':
                case 'W':
                case 'K':
                case 'M':
                case 'B':
                case 'D':
                case 'H':
                case 'V':
                    sequence.push_back('N');
                    break;

                default:
                    throw std::runtime_error(
                        "Unsupported FASTA symbol."
                    );
            }
        }
    }

    return sequence;
}

bool all_acgt(
    const std::string_view sequence
) {
    return std::all_of(
        sequence.begin(),
        sequence.end(),
        [](const char base) {
            return
                base == 'A' ||
                base == 'C' ||
                base == 'G' ||
                base == 'T';
        }
    );
}

double percentile(
    const std::vector<std::uint64_t>& sorted_values,
    const double fraction
) {
    if (sorted_values.empty()) {
        return 0.0;
    }

    const double position =
        fraction *
        static_cast<double>(
            sorted_values.size() - 1
        );

    const std::size_t lower =
        static_cast<std::size_t>(
            position
        );

    const std::size_t upper =
        std::min(
            lower + 1,
            sorted_values.size() - 1
        );

    const double remainder =
        position -
        static_cast<double>(lower);

    return
        static_cast<double>(
            sorted_values.at(lower)
        ) *
        (1.0 - remainder) +
        static_cast<double>(
            sorted_values.at(upper)
        ) *
        remainder;
}

}  // namespace

int main(
    int argc,
    char* argv[]
) {
    try {
        if (argc != 2) {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <fasta>\n";

            return 2;
        }

        constexpr std::size_t
            primer_length = 20;

        constexpr std::size_t
            anchor_length = 12;

        constexpr std::size_t
            target_pool = 4096;

        constexpr std::uint64_t
            stride = 104729ULL;

        std::string reference =
            load_fasta(argv[1]);

        const std::size_t available =
            reference.size() -
            primer_length +
            1;

        std::vector<std::string> primers;

        primers.reserve(
            target_pool
        );

        std::uint64_t candidate = 0;

        for (
            std::size_t attempts = 0;
            attempts < available &&
            primers.size() < target_pool;
            ++attempts
        ) {
            candidate =
                (
                    candidate +
                    stride
                ) %
                available;

            const std::string_view window(
                reference.data() + candidate,
                primer_length
            );

            if (!all_acgt(window)) {
                continue;
            }

            primers.emplace_back(
                window
            );
        }

        if (primers.size() != target_pool) {
            throw std::runtime_error(
                "Could not collect requested primer pool."
            );
        }

        primerpair::BidirectionalFMIndex index(
            std::move(reference)
        );

        std::vector<std::uint64_t>
            anchor_counts;

        std::vector<std::uint64_t>
            full_counts;

        anchor_counts.reserve(target_pool);
        full_counts.reserve(target_pool);

        for (const auto& primer : primers) {

            const std::string_view anchor(
                primer.data() +
                    primer.size() -
                    anchor_length,
                anchor_length
            );

            anchor_counts.push_back(
                index.search(anchor).size()
            );

            full_counts.push_back(
                index.search(primer).size()
            );
        }

        std::sort(
            anchor_counts.begin(),
            anchor_counts.end()
        );

        std::sort(
            full_counts.begin(),
            full_counts.end()
        );

        const auto count_above =
            [](
                const std::vector<std::uint64_t>& values,
                const std::uint64_t threshold
            ) {
                return static_cast<std::size_t>(
                    std::count_if(
                        values.begin(),
                        values.end(),
                        [threshold](
                            const std::uint64_t value
                        ) {
                            return value > threshold;
                        }
                    )
                );
            };

        const auto mean =
            [](
                const std::vector<std::uint64_t>& values
            ) {
                long double total = 0.0;

                for (const auto value : values) {
                    total += value;
                }

                return
                    static_cast<double>(
                        total /
                        values.size()
                    );
            };

        std::cout
            << std::fixed
            << std::setprecision(3);

        std::cout
            << "primer_count\t"
            << primers.size()
            << '\n';

        std::cout
            << "anchor_length\t"
            << anchor_length
            << '\n';

        std::cout
            << "anchor_mean\t"
            << mean(anchor_counts)
            << '\n';

        std::cout
            << "anchor_p50\t"
            << percentile(anchor_counts, 0.50)
            << '\n';

        std::cout
            << "anchor_p75\t"
            << percentile(anchor_counts, 0.75)
            << '\n';

        std::cout
            << "anchor_p90\t"
            << percentile(anchor_counts, 0.90)
            << '\n';

        std::cout
            << "anchor_p95\t"
            << percentile(anchor_counts, 0.95)
            << '\n';

        std::cout
            << "anchor_p99\t"
            << percentile(anchor_counts, 0.99)
            << '\n';

        std::cout
            << "anchor_max\t"
            << anchor_counts.back()
            << '\n';

        for (const auto threshold :
             {1ULL, 10ULL, 100ULL,
              1000ULL, 10000ULL}) {

            const std::size_t number =
                count_above(
                    anchor_counts,
                    threshold
                );

            std::cout
                << "anchor_gt_"
                << threshold
                << '\t'
                << number
                << '\t'
                << (
                    100.0 *
                    static_cast<double>(number) /
                    primers.size()
                )
                << "%\n";
        }

        std::cout
            << "full20_mean\t"
            << mean(full_counts)
            << '\n';

        std::cout
            << "full20_p50\t"
            << percentile(full_counts, 0.50)
            << '\n';

        std::cout
            << "full20_p90\t"
            << percentile(full_counts, 0.90)
            << '\n';

        std::cout
            << "full20_p95\t"
            << percentile(full_counts, 0.95)
            << '\n';

        std::cout
            << "full20_p99\t"
            << percentile(full_counts, 0.99)
            << '\n';

        std::cout
            << "full20_max\t"
            << full_counts.back()
            << '\n';

        return 0;

    } catch (const std::exception& error) {

        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
