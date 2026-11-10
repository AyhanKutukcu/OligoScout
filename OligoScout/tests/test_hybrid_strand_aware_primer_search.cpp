#include "primerpair/anchor_candidate_search.hpp"
#include "primerpair/batched_anchor_lookup.hpp"
#include "primerpair/batched_candidate_search.hpp"
#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/hybrid_batched_primer_search.hpp"
#include "primerpair/hybrid_strand_aware_primer_search.hpp"
#include "primerpair/ip_bwt_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/search_difficulty.hpp"
#include "primerpair/single_primer_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect(
    const bool condition,
    const std::string& name
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " +
            name
        );
    }

    std::cout
        << "[PASS] "
        << name
        << '\n';
}


void normalize(
    std::vector<
        primerpair::OrientedPrimerSearchHit
    >& hits
) {
    std::sort(
        hits.begin(),
        hits.end(),
        [](
            const auto& lhs,
            const auto& rhs
        ) {
            if (
                lhs.position !=
                rhs.position
            ) {
                return
                    lhs.position <
                    rhs.position;
            }

            if (
                lhs.orientation !=
                rhs.orientation
            ) {
                return
                    static_cast<int>(
                        lhs.orientation
                    ) <
                    static_cast<int>(
                        rhs.orientation
                    );
            }

            if (
                lhs.mismatches !=
                rhs.mismatches
            ) {
                return
                    lhs.mismatches <
                    rhs.mismatches;
            }

            return
                lhs.mismatch_mask <
                rhs.mismatch_mask;
        }
    );

    hits.erase(
        std::unique(
            hits.begin(),
            hits.end()
        ),
        hits.end()
    );
}

}  // namespace


int main() {
    try {
        using namespace primerpair;

        constexpr std::size_t
            anchor_length = 12;

        const std::string primer =
            "ACGTTGCAACGTACGT"
            "GATCTGCA";

        const std::string reverse =
            reverse_complement(
                primer
            );

        std::string reverse_5prime_mismatch =
            reverse;

        /*
         * End of reverse_query =
         * original biological 5-prime side.
         */
        reverse_5prime_mismatch.back() =
            reverse_5prime_mismatch.back() == 'A'
                ? 'C'
                : 'A';

        std::string reverse_3prime_mismatch =
            reverse;

        /*
         * Beginning of reverse_query =
         * original biological 3-prime exact anchor.
         */
        reverse_3prime_mismatch.front() =
            reverse_3prime_mismatch.front() == 'A'
                ? 'C'
                : 'A';

        const std::string reference =
            "GCGCGTATATGC" +
            primer +
            "CGATCGATCGAT" +
            reverse +
            "TGCATGCATGCA" +
            reverse_5prime_mismatch +
            "GATCGATCGATC" +
            reverse_3prime_mismatch +
            "CGCGATATCGCG";


        PackedReference packed(
            reference
        );

        BidirectionalFMIndex bifm(
            reference
        );

        IPBWTIndex ipbwt(
            reference,
            anchor_length
        );

        SearchDifficultyEstimator estimator(
            bifm
        );

        BatchedAnchorLookup anchor_lookup(
            ipbwt,
            estimator
        );

        AnchorCandidateSearcher verifier(
            bifm,
            packed
        );

        BatchedCandidateSearchEngine
            candidate_engine(
                ipbwt,
                anchor_lookup,
                verifier
            );

        SinglePrimerSearchEngine
            single_engine(
                bifm,
                packed
            );

        HybridBatchedPrimerSearchEngine
            forward_hybrid(
                candidate_engine,
                single_engine
            );

        StrandAwarePrimerSearchEngine
            legacy_strand(
                bifm,
                packed
            );

        HybridStrandAwarePrimerSearchEngine
            hybrid(
                forward_hybrid,
                anchor_lookup,
                ipbwt,
                packed,
                legacy_strand
            );


        /*
         * k=1 guarantees DirectBranching.
         */
        const std::vector<
            HybridStrandAwarePrimerRequest
        > branching_requests{
            {
                primer,
                1
            }
        };

        const auto branching =
            hybrid.search(
                branching_requests,
                anchor_length
            );

        expect(
            branching.size() == 1,
            "Branching result count"
        );

        auto expected_k1 =
            legacy_strand.search(
                primer,
                anchor_length,
                1
            );

        auto observed_k1 =
            branching.at(0).hits;

        normalize(
            expected_k1.hits
        );

        normalize(
            observed_k1
        );

        expect(
            expected_k1.hits ==
                observed_k1,
            "Hybrid branching equals legacy strand search"
        );

        expect(
            !branching.at(0)
                .reverse_candidate_backend,
            "Reverse k1 uses DirectBranching"
        );


        /*
         * k=3 on a low-frequency anchor should
         * exercise the IP-BWT candidate route.
         */
        const std::vector<
            HybridStrandAwarePrimerRequest
        > candidate_requests{
            {
                primer,
                3
            }
        };

        const auto candidate =
            hybrid.search(
                candidate_requests,
                anchor_length
            );

        expect(
            candidate.size() == 1,
            "Candidate result count"
        );

        auto expected_k3 =
            legacy_strand.search(
                primer,
                anchor_length,
                3
            );

        auto observed_k3 =
            candidate.at(0).hits;

        normalize(
            expected_k3.hits
        );

        normalize(
            observed_k3
        );

        expect(
            expected_k3.hits ==
                observed_k3,
            "Hybrid candidate equals legacy strand search"
        );

        const auto bad =
            std::find_if(
                observed_k3.begin(),
                observed_k3.end(),
                [&](
                    const auto& hit
                ) {
                    const std::size_t bad_position =
                        reference.find(
                            reverse_3prime_mismatch
                        );

                    return
                        bad_position !=
                            std::string::npos &&
                        hit.position ==
                            bad_position &&
                        hit.orientation ==
                            PrimerOrientation::Reverse;
                }
            );

        expect(
            bad ==
                observed_k3.end(),
            "Reverse biological 3-prime mismatch rejected"
        );

        const auto allowed_position =
            reference.find(
                reverse_5prime_mismatch
            );

        const bool allowed =
            std::any_of(
                observed_k3.begin(),
                observed_k3.end(),
                [&](
                    const auto& hit
                ) {
                    return
                        allowed_position !=
                            std::string::npos &&
                        hit.position ==
                            allowed_position &&
                        hit.orientation ==
                            PrimerOrientation::Reverse &&
                        hit.mismatches == 1;
                }
            );

        expect(
            allowed,
            "Reverse biological 5-prime mismatch accepted"
        );

        std::cout
            << "reverse_candidate_backend\t"
            << (
                candidate.at(0)
                    .reverse_candidate_backend
                    ? "YES"
                    : "NO"
            )
            << '\n';

        std::cout
            << "ALL_CHECKS\tYES\n";

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
