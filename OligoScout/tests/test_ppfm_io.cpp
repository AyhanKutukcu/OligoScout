#include "primerpair/ppfm_io.hpp"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect(
    const bool condition,
    const char* message
) {
    if (!condition) {
        throw std::runtime_error(
            message
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string sequence =
            "ACGTACGT"
            "NNNN"
            "GACCTTAGGCTA"
            "ACGTACGT";

        const std::string pattern =
            "ACGTACGT";

        PackedReference reference(
            sequence
        );

        BidirectionalFMIndex index(
            sequence,
            4
        );

        const Interval before_interval =
            index
                .forward_index()
                .backward_search(
                    pattern
                );

        const auto before_positions =
            index
                .forward_index()
                .locate(
                    before_interval
                );

        expect(
            before_positions ==
                std::vector<std::uint64_t>{
                    0,
                    24
                },
            "Baseline exact locate is correct"
        );

        const auto before_bi =
            index.search(
                pattern
            );

        expect(
            before_bi.forward.size() == 2,
            "Baseline bidirectional occurrence count is two"
        );

        const auto unique_suffix =
            std::chrono::
                high_resolution_clock::
                now().
                time_since_epoch().
                count();

        const auto path =
            std::filesystem::
                temp_directory_path()
            /
            (
                "primerpair_ppfm_" +
                std::to_string(
                    unique_suffix
                ) +
                ".ppfm"
            );

        PpfmIO::save_shard(
            path,
            "chrTest",
            reference,
            index
        );

        expect(
            std::filesystem::exists(
                path
            ),
            "PPFM file created"
        );

        expect(
            std::filesystem::file_size(
                path
            ) > 0,
            "PPFM file is non-empty"
        );

        auto loaded =
            PpfmIO::load_shard(
                path
            );

        expect(
            loaded.chromosome ==
                "chrTest",
            "Chromosome metadata restored"
        );

        expect(
            loaded.reference.size() ==
                reference.size(),
            "Reference length restored"
        );

        for (
            std::uint64_t i = 0;
            i < reference.size();
            ++i
        ) {
            expect(
                loaded.reference
                    .base_at(i)
                ==
                reference.base_at(i),
                "Packed reference base restored"
            );
        }

        expect(
            loaded.index
                .forward_index()
                .suffix_array_sample_rate()
            ==
            4,
            "Forward SA sample rate restored"
        );

        expect(
            loaded.index
                .reverse_index()
                .suffix_array_sample_rate()
            ==
            4,
            "Reverse SA sample rate restored"
        );

        const Interval after_interval =
            loaded.index
                .forward_index()
                .backward_search(
                    pattern
                );

        expect(
            after_interval.begin ==
                before_interval.begin
            &&
            after_interval.end ==
                before_interval.end,
            "Backward-search interval restored"
        );

        const auto after_positions =
            loaded.index
                .forward_index()
                .locate(
                    after_interval
                );

        expect(
            after_positions ==
                before_positions,
            "Locate results identical after reload"
        );

        const auto after_bi =
            loaded.index.search(
                pattern
            );

        expect(
            after_bi.forward.begin ==
                before_bi.forward.begin
            &&
            after_bi.forward.end ==
                before_bi.forward.end
            &&
            after_bi.reverse.begin ==
                before_bi.reverse.begin
            &&
            after_bi.reverse.end ==
                before_bi.reverse.end
            &&
            after_bi.length ==
                before_bi.length,
            "Bidirectional search state identical after reload"
        );

        expect(
            loaded.reference
                .bounded_hamming_distance(
                    0,
                    pattern,
                    0
                )
            ==
            0,
            "Packed reference verification works after reload"
        );

        std::error_code error;

        std::filesystem::remove(
            path,
            error
        );

        std::cout
            << "PPFM round-trip tests passed.\n";

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
