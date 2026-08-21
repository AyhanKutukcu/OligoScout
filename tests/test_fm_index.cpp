#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "primerpair/fm_index.hpp"

namespace {

bool check(
    const bool condition,
    const std::string& name
) {
    if (condition) {
        std::cout
            << "[PASS] "
            << name
            << '\n';

        return true;
    }

    std::cerr
        << "[FAIL] "
        << name
        << '\n';

    return false;
}

}  // namespace

int main() {
    int failures = 0;

    /*
     * Küçük sample rate, LF tabanlı locate'ın
     * birden fazla sample arasında çalışmasını test eder.
     */
    const primerpair::FMIndex gattaca_index(
        "GATTACA",
        4
    );

    if (!check(
            gattaca_index.bwt_string() == "ACTGA$TA",
            "GATTACA BWT"
        )) {
        ++failures;
    }

    const std::vector<std::uint64_t>
        expected_lf{
            1, 4, 6, 5,
            2, 0, 7, 3
        };

    bool lf_mapping_correct = true;

    for (std::uint64_t row = 0;
         row < expected_lf.size();
         ++row) {

        if (
            gattaca_index.lf(row) !=
            expected_lf.at(
                static_cast<std::size_t>(row)
            )
        ) {
            lf_mapping_correct = false;
            break;
        }
    }

    if (!check(
            lf_mapping_correct,
            "GATTACA LF mapping"
        )) {
        ++failures;
    }

    if (!check(
            gattaca_index.suffix_array_sample_rate() == 4,
            "Suffix-array sample rate"
        )) {
        ++failures;
    }

    if (!check(
            gattaca_index.sampled_sa_count() == 2,
            "GATTACA sampled SA count"
        )) {
        ++failures;
    }

    /*
     * Bütün FM satırlarını locate ederek tam suffix
     * array'in artık gerekli olmadığını doğrular.
     *
     * Sentinel koordinatı locate() tarafından çıkarılır.
     */
    const primerpair::Interval all_rows{
        .begin = 0,
        .end = static_cast<std::uint64_t>(
            gattaca_index.bwt_size()
        )
    };

    const std::vector<std::uint64_t>
        all_positions =
            gattaca_index.locate(all_rows);

    std::vector<std::uint64_t>
        expected_all_positions(7);

    std::iota(
        expected_all_positions.begin(),
        expected_all_positions.end(),
        std::uint64_t{0}
    );

    if (!check(
            all_positions ==
                expected_all_positions,
            "LF-based locate for all rows"
        )) {
        ++failures;
    }

    bool lf_out_of_range_rejected = false;

    try {
        static_cast<void>(
            gattaca_index.lf(8)
        );
    } catch (const std::out_of_range&) {
        lf_out_of_range_rejected = true;
    }

    if (!check(
            lf_out_of_range_rejected,
            "LF out-of-range rejection"
        )) {
        ++failures;
    }

    bool zero_sample_rate_rejected = false;

    try {
        const primerpair::FMIndex invalid_index(
            "ACGT",
            0
        );

        static_cast<void>(invalid_index);
    } catch (const std::invalid_argument&) {
        zero_sample_rate_rejected = true;
    }

    if (!check(
            zero_sample_rate_rejected,
            "Zero sample-rate rejection"
        )) {
        ++failures;
    }

    const std::string reference =
        "ACGTACGTACGTGATTACAGATTACACCCCC"
        "GGGGGTTTTTAAAAAACGTACGTACGT";

    const std::string pattern =
        "ACGTACGTACGT";

    const primerpair::FMIndex index(
        reference
    );

    if (!check(
            index.suffix_array_sample_rate() ==
                primerpair::FMIndex::
                    kDefaultSuffixArraySampleRate,
            "Default sample rate"
        )) {
        ++failures;
    }

    if (!check(
            index.sampled_sa_count() ==
                (
                    reference.size() + 1
                    +
                    primerpair::FMIndex::
                        kDefaultSuffixArraySampleRate
                    - 1
                )
                /
                primerpair::FMIndex::
                    kDefaultSuffixArraySampleRate,
            "58-bp sampled SA count"
        )) {
        ++failures;
    }

    const primerpair::Interval interval =
        index.backward_search(pattern);

    if (!check(
            interval.begin == 11 &&
            interval.end == 13 &&
            interval.size() == 2,
            "Exact interval"
        )) {
        ++failures;
    }

    const std::vector<std::uint64_t> positions =
        index.locate(interval);

    const std::vector<std::uint64_t>
        expected_positions{0, 46};

    if (!check(
            positions == expected_positions,
            "Sampled-SA genomic positions"
        )) {
        ++failures;
    }

    const primerpair::Interval absent_interval =
        index.backward_search(
            "TTTTCCCCAAAA"
        );

    if (!check(
            absent_interval.empty(),
            "Absent pattern interval"
        )) {
        ++failures;
    }

    if (!check(
            index.locate(absent_interval).empty(),
            "Absent pattern coordinates"
        )) {
        ++failures;
    }

    const primerpair::Interval lowercase_interval =
        index.backward_search(
            "acgtacgtacgt"
        );

    if (!check(
            index.locate(lowercase_interval) ==
                expected_positions,
            "Lowercase query normalization"
        )) {
        ++failures;
    }

    bool invalid_nucleotide_rejected = false;

    try {
        static_cast<void>(
            index.backward_search(
                "ACGTX"
            )
        );
    } catch (const std::invalid_argument&) {
        invalid_nucleotide_rejected = true;
    }

    if (!check(
            invalid_nucleotide_rejected,
            "Invalid nucleotide rejection"
        )) {
        ++failures;
    }

    if (failures != 0) {
        std::cerr
            << failures
            << " FM-index test(s) failed.\n";

        return 1;
    }

    std::cout
        << "All sampled-SA FM-index tests passed.\n";

    return 0;
}
