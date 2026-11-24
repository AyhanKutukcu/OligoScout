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

/*
 * Tek primer ürünleri (single-primer products).
 *
 * Aynı primerin bir sitesi ileri, başka bir sitesi geri iplikte bağlanıp
 * içe bakıyorsa, o primer tek başına amplikon üretir:
 *
 *   P --->       <--- P
 *
 * Primer-BLAST bu ürünleri de raporlar. assemble_primer_pair_hits()
 * yalnızca P1->/<-P2 ve P2->/<-P1 birleştirir; bu fonksiyon P->/<-P
 * yönelimini aynı geometri kurallarıyla ekler. left_primer ve
 * right_primer ikisi de `identity` olur.
 *
 * Complexity: O(H + K)
 */
[[nodiscard]]
std::vector<PrimerPairHit>
assemble_single_primer_hits(
    std::string_view primer,
    const std::vector<
        OrientedPrimerSearchHit
    >& primer_hits,
    PrimerIdentity identity,

    std::uint64_t min_amplicon_length,
    std::uint64_t max_amplicon_length
);

}  // namespace primerpair
