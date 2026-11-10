#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/primer_pair_search.hpp"

namespace primerpair {


enum class MultiplexPrimerSlot :
    std::uint8_t {
    Primer1 = 0,
    Primer2 = 1
};


struct MultiplexPrimerPairRequest {
    std::string_view primer1{};
    std::string_view primer2{};

    std::size_t max_mismatches{3};

    std::uint64_t min_amplicon_length{50};
    std::uint64_t max_amplicon_length{3000};
};


/*
 * Cross-pair PCR product.
 *
 * forward_*:
 *   genome üzerinde ---> yönündeki primer.
 *
 * reverse_*:
 *   genome üzerinde <--- yönündeki primer.
 */
struct MultiplexCrossAmplicon {
    std::size_t forward_pair_index{0};

    MultiplexPrimerSlot forward_slot{
        MultiplexPrimerSlot::Primer1
    };

    std::size_t reverse_pair_index{0};

    MultiplexPrimerSlot reverse_slot{
        MultiplexPrimerSlot::Primer1
    };

    std::uint64_t amplicon_start{0};
    std::uint64_t amplicon_end_exclusive{0};
    std::uint64_t amplicon_length{0};

    std::size_t forward_mismatches{0};
    std::size_t reverse_mismatches{0};

    std::uint64_t forward_mismatch_mask{0};
    std::uint64_t reverse_mismatch_mask{0};

    bool operator==(
        const MultiplexCrossAmplicon&
    ) const = default;
};


struct MultiplexSearchStats {
    std::size_t pair_requests{0};

    std::size_t total_primer_slots{0};

    /*
     * Gerçekte HybridStrandAware motoruna gönderilen
     * distinct (sequence,k) query sayısı.
     */
    std::size_t unique_primer_queries{0};

    std::size_t reused_primer_slots{0};

    std::size_t intended_join_computations{0};

    /*
     * Pair_i x Pair_j için:
     *
     *   P1-P1
     *   P1-P2
     *   P2-P1
     *   P2-P2
     */
    std::size_t cross_slot_pair_requests{0};

    /*
     * Aynı unique primer kombinasyonu birden fazla
     * panel slotunda görülürse sweep join tekrar
     * hesaplanmaz.
     */
    std::size_t unique_cross_join_computations{0};

    std::size_t reused_cross_join_requests{0};

    std::size_t cross_amplicon_records{0};
};


struct MultiplexPrimerSearchResult {
    /*
     * Girdi pair sırasıyla birebir eşleşir.
     */
    std::vector<PrimerPairSearchResult>
        intended_pairs{};

    std::vector<MultiplexCrossAmplicon>
        cross_amplicons{};

    MultiplexSearchStats stats{};
};


class MultiplexPrimerSearchEngine {
public:
    explicit
    MultiplexPrimerSearchEngine(
        const HybridStrandAwarePrimerSearchEngine&
            primer_engine
    ) noexcept;

    /*
     * Pipeline:
     *
     * panel
     *   ↓
     * distinct (primer sequence, mismatch budget)
     *   ↓
     * ONE HybridStrandAware batch
     *   ↓
     * unique primer hit cache
     *   ↓
     * intended pair sweep joins
     *   ↓
     * optional cross-pair sweep joins
     *
     * Cross pair range panel-genelidir.
     */
    [[nodiscard]]
    MultiplexPrimerSearchResult
    search(
        const std::vector<
            MultiplexPrimerPairRequest
        >& requests,

        std::size_t anchor_length = 12,

        bool include_cross_pairs = true,

        std::uint64_t
            cross_min_amplicon_length = 50,

        std::uint64_t
            cross_max_amplicon_length = 3000
    ) const;

private:
    const HybridStrandAwarePrimerSearchEngine&
        primer_engine_;
};


}  // namespace primerpair
