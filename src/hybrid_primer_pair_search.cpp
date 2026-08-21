#include "primerpair/hybrid_primer_pair_search.hpp"

#include "primerpair/primer_pair_assembler.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace primerpair {


HybridPrimerPairSearchEngine::
HybridPrimerPairSearchEngine(
    const HybridStrandAwarePrimerSearchEngine&
        primer_engine
) noexcept
    : primer_engine_(
          primer_engine
      ) {
}


std::vector<HybridPrimerPairResult>
HybridPrimerPairSearchEngine::search(
    const std::vector<
        HybridPrimerPairRequest
    >& requests,
    const std::size_t anchor_length
) const {
    if (requests.empty()) {
        return {};
    }


    std::vector<
        HybridStrandAwarePrimerRequest
    > primer_requests;

    primer_requests.reserve(
        requests.size() * 2
    );


    /*
     * Pair panelini tek both-strand primer batch'ine
     * dönüştür.
     */
    for (
        const auto& request :
        requests
    ) {
        primer_requests.push_back(
            HybridStrandAwarePrimerRequest{
                request.primer1,
                request.max_mismatches
            }
        );

        primer_requests.push_back(
            HybridStrandAwarePrimerRequest{
                request.primer2,
                request.max_mismatches
            }
        );
    }


    const auto primer_results =
        primer_engine_.search(
            primer_requests,
            anchor_length
        );


    if (
        primer_results.size() !=
        primer_requests.size()
    ) {
        throw std::logic_error(
            "Hybrid primer-pair primer "
            "result-count mismatch."
        );
    }


    std::vector<
        HybridPrimerPairResult
    > results;

    results.reserve(
        requests.size()
    );


    for (
        std::size_t pair_index = 0;
        pair_index < requests.size();
        ++pair_index
    ) {
        const std::size_t
            primer1_index =
                pair_index * 2;

        const std::size_t
            primer2_index =
                primer1_index + 1;


        const auto& request =
            requests.at(
                pair_index
            );

        const auto& primer1_result =
            primer_results.at(
                primer1_index
            );

        const auto& primer2_result =
            primer_results.at(
                primer2_index
            );


        auto assembled =
            assemble_primer_pair_hits(
                request.primer1,
                primer1_result.hits,

                request.primer2,
                primer2_result.hits,

                request.min_amplicon_length,
                request.max_amplicon_length
            );


        results.push_back(
            HybridPrimerPairResult{
                primer1_result,
                primer2_result,
                std::move(
                    assembled
                )
            }
        );
    }


    return results;
}


}  // namespace primerpair
