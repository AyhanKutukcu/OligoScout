#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

namespace primerpair {

/*
 * Build a suffix array using prefix doubling with
 * stable radix/counting sort.
 *
 * Input text must:
 *
 *   - be non-empty
 *   - contain symbols from:
 *       $, A, C, G, N, T
 *   - contain exactly one '$'
 *   - have '$' as the final character
 *
 * Output:
 *
 *   suffix_array[row] = suffix start coordinate
 *
 * Chromosome-shard coordinates use uint32_t.
 */
[[nodiscard]]
std::vector<std::uint32_t>
build_suffix_array_prefix_doubling(
    std::string_view text
);

}  // namespace primerpair
