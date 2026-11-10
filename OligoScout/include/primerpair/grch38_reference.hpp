#pragma once

#include <optional>
#include <string>
#include <string_view>

namespace primerpair {


/*
 * Resolve a GRCh38.p14 NCBI primary chromosome
 * record to:
 *
 *   chr1 .. chr22
 *   chrX
 *   chrY
 *
 * Unlocalized, alternate and patch records return
 * std::nullopt.
 *
 * Accession versions are intentionally not
 * hard-coded.
 */
[[nodiscard]]
std::optional<std::string>
grch38_primary_chromosome_alias(
    std::string_view record_name,
    std::string_view description
);


}  // namespace primerpair
