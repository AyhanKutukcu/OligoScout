#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "primerpair/strand_aware_primer_search.hpp"

namespace primerpair {


/*
 * Bir unique primer için daha önce hesaplanmış
 * both-strand genomic hit listesi.
 *
 * Bu katman genom araması yapmaz.
 */
struct GlobalMultiplexPrimerHits {
    std::size_t primer_id{0};
    std::size_t primer_length{0};

    std::span<
        const OrientedPrimerSearchHit
    > hits{};
};


/*
 * Global sweep sonucunda bulunan fiziksel olarak
 * PCR-compatible cross-primer ürünü.
 *
 * forward_primer_id:
 *     genome üzerinde ---> bağlanan primer.
 *
 * reverse_primer_id:
 *     genome üzerinde <--- bağlanan primer.
 */
struct GlobalMultiplexCrossProduct {
    std::size_t forward_primer_id{0};
    std::size_t reverse_primer_id{0};

    std::uint64_t forward_position{0};
    std::uint64_t reverse_position{0};

    std::uint64_t amplicon_start{0};
    std::uint64_t amplicon_end_exclusive{0};
    std::uint64_t amplicon_length{0};

    std::size_t forward_mismatches{0};
    std::size_t reverse_mismatches{0};

    std::uint64_t forward_mismatch_mask{0};
    std::uint64_t reverse_mismatch_mask{0};

    bool operator==(
        const GlobalMultiplexCrossProduct&
    ) const = default;
};


struct GlobalMultiplexCrossJoinStats {
    std::size_t unique_primers{0};

    std::size_t forward_hits{0};
    std::size_t reverse_hits{0};

    /*
     * Amplicon endpoint window'a giren hit-hit
     * kombinasyonları.
     */
    std::size_t window_candidates{0};

    /*
     * Non-overlap + PCR geometry filtresinden
     * geçen ürün sayısı; normalization öncesi.
     */
    std::size_t emitted_products{0};

    std::size_t unique_products{0};
};


struct GlobalMultiplexCrossJoinResult {
    std::vector<
        GlobalMultiplexCrossProduct
    > products{};

    GlobalMultiplexCrossJoinStats stats{};
};


/*
 * Tüm primer-pair kombinasyonlarını P^2 dolaşmak
 * yerine bütün forward ve reverse genomic hitleri
 * global koordinat düzeninde eşleştirir.
 *
 * Endpoint constraint:
 *
 *   forward_position + min_amplicon
 *       <= reverse_end
 *       <=
 *   forward_position + max_amplicon
 *
 * Ek PCR geometry:
 *
 *   reverse_position
 *       >=
 *   forward_position + forward_primer_length
 *
 * Complexity:
 *
 *   O(H log H + C + K)
 *
 * H = toplam hit sayısı
 * C = endpoint window candidate sayısı
 * K = üretilen ürün sayısı
 */
[[nodiscard]]
GlobalMultiplexCrossJoinResult
global_multiplex_cross_join(
    const std::vector<
        GlobalMultiplexPrimerHits
    >& primers,

    std::uint64_t min_amplicon_length,
    std::uint64_t max_amplicon_length
);


}  // namespace primerpair
