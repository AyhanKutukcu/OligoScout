#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/primer_pair_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {

/*
 * Precomputed strand-aware primer hitlerinden
 * PCR-uyumlu primer-pair ürünlerini oluşturur.
 *
 * Bu fonksiyon genom araması YAPMAZ.
 *
 * Join:
 *
 *   FORWARD --->       <--- REVERSE
 *
 * Her iki orientation konfigürasyonu için
 * monoton sweep-line kullanılır.
 *
 * Complexity:
 *
 *   O(F + R + K)
 *
 * K = üretilen PCR-compatible amplicon sayısı.
 */
[[nodiscard]]
PrimerPairSearchResult
assemble_primer_pair_hits(
    std::string_view primer1,
    const std::vector<
        OrientedPrimerSearchHit
    >& primer1_hits,

    std::string_view primer2,
    const std::vector<
        OrientedPrimerSearchHit
    >& primer2_hits,

    std::uint64_t min_amplicon_length,
    std::uint64_t max_amplicon_length
);

}  // namespace primerpair
