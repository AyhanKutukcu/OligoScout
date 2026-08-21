#include "primerpair/persistent_multiplex_primer_search_v2.hpp"

#include "primerpair/bidirectional_fm_index.hpp"
#include "primerpair/packed_reference.hpp"
#include "primerpair/ppfm_io.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/primer_pair_assembler.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <utility>
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


struct TemporaryDirectory {
    std::filesystem::path path;

    ~TemporaryDirectory() {
        std::error_code error;

        std::filesystem::remove_all(
            path,
            error
        );
    }
};


using CrossKey =
    std::tuple<
        std::size_t,
        int,
        std::size_t,
        int,
        std::uint64_t,
        std::uint64_t,
        std::uint64_t,
        std::size_t,
        std::size_t,
        std::uint64_t,
        std::uint64_t
    >;


CrossKey cross_key(
    const primerpair::MultiplexCrossAmplicon&
        product
) {
    return CrossKey{
        product.forward_pair_index,
        static_cast<int>(
            product.forward_slot
        ),

        product.reverse_pair_index,
        static_cast<int>(
            product.reverse_slot
        ),

        product.amplicon_start,
        product.amplicon_end_exclusive,
        product.amplicon_length,

        product.forward_mismatches,
        product.reverse_mismatches,

        product.forward_mismatch_mask,
        product.reverse_mismatch_mask
    };
}


std::vector<CrossKey>
normalized_keys(
    const std::vector<
        primerpair::MultiplexCrossAmplicon
    >& products
) {
    std::vector<CrossKey> keys;

    keys.reserve(
        products.size()
    );

    for (
        const auto& product :
        products
    ) {
        keys.push_back(
            cross_key(
                product
            )
        );
    }

    std::sort(
        keys.begin(),
        keys.end()
    );

    keys.erase(
        std::unique(
            keys.begin(),
            keys.end()
        ),
        keys.end()
    );

    return keys;
}


/*
 * Independent chromosome-local production oracle.
 *
 * This deliberately uses:
 *
 *   StrandAwarePrimerSearchEngine
 *       +
 *   assemble_primer_pair_hits()
 *
 * and performs the old logical pair/slot enumeration.
 *
 * It therefore validates the persistent V2 global sweep
 * against the existing production V1 pair assembler.
 */
std::vector<
    primerpair::MultiplexCrossAmplicon
>
build_local_oracle(
    const std::string& sequence,

    const std::vector<
        primerpair::MultiplexPrimerPairRequest
    >& requests,

    const std::size_t anchor_length,

    const std::uint64_t min_amplicon_length,

    const std::uint64_t max_amplicon_length
) {
    using namespace primerpair;


    PackedReference reference(
        sequence
    );

    BidirectionalFMIndex index(
        sequence,
        8
    );

    StrandAwarePrimerSearchEngine engine(
        index,
        reference
    );


    /*
     * Synthetic test uses unique primer sequences,
     * therefore one hit list per panel slot is enough.
     */
    std::vector<
        StrandAwarePrimerSearchResult
    > slot_results;

    slot_results.reserve(
        requests.size() * 2
    );


    for (
        const auto& request :
        requests
    ) {
        slot_results.push_back(
            engine.search(
                request.primer1,
                anchor_length,
                request.max_mismatches
            )
        );

        slot_results.push_back(
            engine.search(
                request.primer2,
                anchor_length,
                request.max_mismatches
            )
        );
    }


    std::vector<
        MultiplexCrossAmplicon
    > output;


    for (
        std::size_t pair_a = 0;
        pair_a < requests.size();
        ++pair_a
    ) {
        for (
            std::size_t pair_b =
                pair_a + 1;

            pair_b < requests.size();
            ++pair_b
        ) {
            for (
                std::size_t slot_a = 0;
                slot_a < 2;
                ++slot_a
            ) {
                for (
                    std::size_t slot_b = 0;
                    slot_b < 2;
                    ++slot_b
                ) {
                    const std::size_t index_a =
                        pair_a * 2 +
                        slot_a;

                    const std::size_t index_b =
                        pair_b * 2 +
                        slot_b;


                    const std::string_view
                        sequence_a =
                            slot_a == 0
                            ?
                            requests
                                .at(pair_a)
                                .primer1
                            :
                            requests
                                .at(pair_a)
                                .primer2;


                    const std::string_view
                        sequence_b =
                            slot_b == 0
                            ?
                            requests
                                .at(pair_b)
                                .primer1
                            :
                            requests
                                .at(pair_b)
                                .primer2;


                    const auto assembled =
                        assemble_primer_pair_hits(
                            sequence_a,
                            slot_results
                                .at(index_a)
                                .hits,

                            sequence_b,
                            slot_results
                                .at(index_b)
                                .hits,

                            min_amplicon_length,
                            max_amplicon_length
                        );


                    for (
                        const auto& product :
                        assembled.amplicons
                    ) {
                        std::size_t
                            forward_pair = 0;

                        std::size_t
                            reverse_pair = 0;

                        MultiplexPrimerSlot
                            forward_slot =
                                MultiplexPrimerSlot::
                                    Primer1;

                        MultiplexPrimerSlot
                            reverse_slot =
                                MultiplexPrimerSlot::
                                    Primer1;


                        if (
                            product.left_primer ==
                                PrimerIdentity::
                                    Primer1
                            &&
                            product.right_primer ==
                                PrimerIdentity::
                                    Primer2
                        ) {
                            forward_pair =
                                pair_a;

                            reverse_pair =
                                pair_b;

                            forward_slot =
                                slot_a == 0
                                ?
                                MultiplexPrimerSlot::
                                    Primer1
                                :
                                MultiplexPrimerSlot::
                                    Primer2;

                            reverse_slot =
                                slot_b == 0
                                ?
                                MultiplexPrimerSlot::
                                    Primer1
                                :
                                MultiplexPrimerSlot::
                                    Primer2;

                        } else if (
                            product.left_primer ==
                                PrimerIdentity::
                                    Primer2
                            &&
                            product.right_primer ==
                                PrimerIdentity::
                                    Primer1
                        ) {
                            forward_pair =
                                pair_b;

                            reverse_pair =
                                pair_a;

                            forward_slot =
                                slot_b == 0
                                ?
                                MultiplexPrimerSlot::
                                    Primer1
                                :
                                MultiplexPrimerSlot::
                                    Primer2;

                            reverse_slot =
                                slot_a == 0
                                ?
                                MultiplexPrimerSlot::
                                    Primer1
                                :
                                MultiplexPrimerSlot::
                                    Primer2;

                        } else {
                            throw std::runtime_error(
                                "Unexpected oracle "
                                "primer identities."
                            );
                        }


                        output.push_back(
                            MultiplexCrossAmplicon{
                                forward_pair,
                                forward_slot,

                                reverse_pair,
                                reverse_slot,

                                product
                                    .amplicon_start,

                                product
                                    .amplicon_end_exclusive,

                                product
                                    .amplicon_length,

                                product
                                    .left_mismatches,

                                product
                                    .right_mismatches,

                                product
                                    .left_mismatch_mask,

                                product
                                    .right_mismatch_mask
                            }
                        );
                    }
                }
            }
        }
    }


    return output;
}


const primerpair::
PersistentMultiplexShardResultV2&
find_shard(
    const primerpair::
        PersistentMultiplexSearchResultV2& result,

    const std::string_view chromosome
) {
    for (
        const auto& shard :
        result.shards
    ) {
        if (
            shard.chromosome ==
            chromosome
        ) {
            return shard;
        }
    }

    throw std::runtime_error(
        "Requested shard not found."
    );
}


}  // namespace


int main() {
    using namespace primerpair;

    try {
        /*
         * Four deliberately distinct 20mers.
         *
         * Exact-match mode is used in this regression
         * test so accidental approximate matches cannot
         * obscure chromosome-boundary behavior.
         */
        const std::string primer_a =
            "ACGTTGCAACGATCGTACGA";

        const std::string primer_b =
            "TGCACCTAGGCTAACGTGCA";

        const std::string primer_c =
            "GATCCGATGCTAGTCAGTCA";

        const std::string primer_d =
            "CTAGGTCAGCATCGTACGTT";


        /*
         * Pair 0:
         *   primer_a + primer_b
         *
         * Pair 1:
         *   primer_c + primer_d
         *
         * Important:
         *
         * pair0's compatible sites are split between
         * chrA and chrB.
         *
         * pair1's compatible sites are also split
         * between chrB and chrA.
         *
         * Therefore a buggy whole-genome hit pool could
         * fabricate cross-chromosome intended products.
         */
        const std::vector<
            MultiplexPrimerPairRequest
        > requests{
            MultiplexPrimerPairRequest{
                primer_a,
                primer_b,
                0,
                80,
                150
            },

            MultiplexPrimerPairRequest{
                primer_c,
                primer_d,
                0,
                80,
                150
            }
        };


        const std::string sequence_a =
            std::string(
                50,
                'N'
            )
            +
            primer_a
            +
            std::string(
                60,
                'N'
            )
            +
            reverse_complement(
                primer_d
            )
            +
            std::string(
                50,
                'N'
            );


        const std::string sequence_b =
            std::string(
                50,
                'N'
            )
            +
            primer_c
            +
            std::string(
                70,
                'N'
            )
            +
            reverse_complement(
                primer_b
            )
            +
            std::string(
                50,
                'N'
            );


        expect(
            sequence_a.size() == 200,
            "chrA synthetic length correct"
        );

        expect(
            sequence_b.size() == 210,
            "chrB synthetic length correct"
        );


        const auto suffix =
            std::chrono::
                high_resolution_clock::
                now().
                time_since_epoch().
                count();


        TemporaryDirectory temporary{
            std::filesystem::
                temp_directory_path()
            /
            (
                "primerpair_persistent_"
                "multiplex_v2_test_" +
                std::to_string(
                    suffix
                )
            )
        };


        std::filesystem::
            create_directories(
                temporary.path
            );


        const std::array<
            std::pair<
                std::string,
                std::string
            >,
            2
        > chromosomes{
            std::pair{
                std::string("chrA"),
                sequence_a
            },

            std::pair{
                std::string("chrB"),
                sequence_b
            }
        };


        /*
         * Create two real persistent PPFM shards.
         */
        for (
            const auto& [
                chromosome,
                sequence
            ] :
            chromosomes
        ) {
            PackedReference reference(
                sequence
            );

            BidirectionalFMIndex index(
                sequence,
                8
            );

            PpfmIO::save_shard(
                temporary.path
                    /
                    (
                        chromosome +
                        ".sa8.ppfm"
                    ),

                chromosome,
                reference,
                index
            );
        }


        const auto manifest_path =
            temporary.path
            /
            "manifest.sa8.tsv";


        {
            std::ofstream output(
                manifest_path
            );

            if (!output) {
                throw std::runtime_error(
                    "Cannot create synthetic manifest."
                );
            }


            output
                << "chromosome\t"
                << "file_bytes\t"
                << "sha256\t"
                << "status\n";


            for (
                const auto& [
                    chromosome,
                    sequence
                ] :
                chromosomes
            ) {
                (void)sequence;

                const auto path =
                    temporary.path
                    /
                    (
                        chromosome +
                        ".sa8.ppfm"
                    );


                output
                    << chromosome
                    << '\t'
                    << std::filesystem::
                        file_size(
                            path
                        )
                    << '\t'
                    << std::string(
                        64,
                        '0'
                    )
                    << '\t'
                    << "BUILT\n";
            }
        }


        auto manifest =
            PpfmManifest::load(
                manifest_path,
                temporary.path
            );


        expect(
            manifest.size() == 2,
            "Synthetic manifest has two shards"
        );


        /*
         * Capacity one is intentional:
         *
         * chrA load
         * chrB load -> chrA eviction
         */
        PersistentMultiplexPrimerSearchEngineV2
            engine(
                std::move(
                    manifest
                ),
                1,
                8
            );


        const auto result =
            engine.search(
                requests,
                12,
                true,
                80,
                150
            );


        expect(
            result.shards.size() == 2,
            "Persistent V2 returns two shards"
        );


        const auto& shard_a =
            find_shard(
                result,
                "chrA"
            );


        const auto& shard_b =
            find_shard(
                result,
                "chrB"
            );


        /*
         * Exact single-primer hit expectations.
         *
         * chrA:
         *   primer_a Forward @ 50
         *   primer_d Reverse @ 130
         *
         * chrB:
         *   primer_c Forward @ 50
         *   primer_b Reverse @ 140
         */
        expect(
            shard_a.stats
                .total_primer_hits == 2,
            "chrA has exactly two primer hits"
        );

        expect(
            shard_a.stats
                .forward_primer_hits == 1,
            "chrA has one forward primer hit"
        );

        expect(
            shard_a.stats
                .reverse_primer_hits == 1,
            "chrA has one reverse primer hit"
        );


        expect(
            shard_b.stats
                .total_primer_hits == 2,
            "chrB has exactly two primer hits"
        );

        expect(
            shard_b.stats
                .forward_primer_hits == 1,
            "chrB has one forward primer hit"
        );

        expect(
            shard_b.stats
                .reverse_primer_hits == 1,
            "chrB has one reverse primer hit"
        );


        /*
         * No intended pair is complete inside either
         * chromosome.
         *
         * This explicitly protects the biological
         * no-cross-chromosome-pairing invariant.
         */
        expect(
            shard_a.intended_pairs.size() == 2,
            "chrA contains two intended result slots"
        );

        expect(
            shard_b.intended_pairs.size() == 2,
            "chrB contains two intended result slots"
        );


        for (
            const auto& intended :
            shard_a.intended_pairs
        ) {
            expect(
                intended.amplicons.empty(),
                "chrA has no cross-chromosome "
                "intended product"
            );
        }


        for (
            const auto& intended :
            shard_b.intended_pairs
        ) {
            expect(
                intended.amplicons.empty(),
                "chrB has no cross-chromosome "
                "intended product"
            );
        }


        /*
         * Exact V2 physical products.
         */
        expect(
            shard_a.cross_amplicons.size() == 1,
            "chrA has exactly one cross product"
        );

        expect(
            shard_b.cross_amplicons.size() == 1,
            "chrB has exactly one cross product"
        );


        const auto& product_a =
            shard_a.cross_amplicons.front();


        expect(
            product_a.forward_pair_index == 0,
            "chrA forward pair is Pair0"
        );

        expect(
            product_a.forward_slot ==
                MultiplexPrimerSlot::Primer1,
            "chrA forward slot is Pair0 Primer1"
        );

        expect(
            product_a.reverse_pair_index == 1,
            "chrA reverse pair is Pair1"
        );

        expect(
            product_a.reverse_slot ==
                MultiplexPrimerSlot::Primer2,
            "chrA reverse slot is Pair1 Primer2"
        );

        expect(
            product_a.amplicon_start == 50,
            "chrA amplicon start correct"
        );

        expect(
            product_a.amplicon_end_exclusive == 150,
            "chrA amplicon end correct"
        );

        expect(
            product_a.amplicon_length == 100,
            "chrA amplicon length correct"
        );

        expect(
            product_a.forward_mismatches == 0 &&
            product_a.reverse_mismatches == 0,
            "chrA product is exact-match"
        );


        const auto& product_b =
            shard_b.cross_amplicons.front();


        expect(
            product_b.forward_pair_index == 1,
            "chrB forward pair is Pair1"
        );

        expect(
            product_b.forward_slot ==
                MultiplexPrimerSlot::Primer1,
            "chrB forward slot is Pair1 Primer1"
        );

        expect(
            product_b.reverse_pair_index == 0,
            "chrB reverse pair is Pair0"
        );

        expect(
            product_b.reverse_slot ==
                MultiplexPrimerSlot::Primer2,
            "chrB reverse slot is Pair0 Primer2"
        );

        expect(
            product_b.amplicon_start == 50,
            "chrB amplicon start correct"
        );

        expect(
            product_b.amplicon_end_exclusive == 160,
            "chrB amplicon end correct"
        );

        expect(
            product_b.amplicon_length == 110,
            "chrB amplicon length correct"
        );

        expect(
            product_b.forward_mismatches == 0 &&
            product_b.reverse_mismatches == 0,
            "chrB product is exact-match"
        );


        /*
         * Independent chromosome-local production V1
         * oracle versus persistent Global Sweep V2.
         */
        const auto oracle_a =
            build_local_oracle(
                sequence_a,
                requests,
                12,
                80,
                150
            );


        const auto oracle_b =
            build_local_oracle(
                sequence_b,
                requests,
                12,
                80,
                150
            );


        expect(
            normalized_keys(
                shard_a.cross_amplicons
            )
            ==
            normalized_keys(
                oracle_a
            ),
            "chrA V2 equals chromosome-local "
            "production assembler oracle"
        );


        expect(
            normalized_keys(
                shard_b.cross_amplicons
            )
            ==
            normalized_keys(
                oracle_b
            ),
            "chrB V2 equals chromosome-local "
            "production assembler oracle"
        );


        /*
         * Whole-result aggregate invariants.
         */
        expect(
            result.total_primer_hits() == 4,
            "Whole persistent search has four hits"
        );

        expect(
            result.total_cross_amplicons() == 2,
            "Whole persistent search has two "
            "cross products"
        );

        expect(
            result.total_window_candidates() == 2,
            "Global sweep sees exactly two "
            "physical candidates"
        );


        /*
         * LRU/cache behavior.
         */
        expect(
            engine.cache().capacity() == 1,
            "Persistent test cache capacity is one"
        );

        expect(
            engine.cache().load_count() == 2,
            "Both PPFM shards loaded exactly once"
        );

        expect(
            engine.cache().eviction_count() == 1,
            "Capacity-one cache performs one eviction"
        );

        expect(
            engine.cache().size() == 1,
            "Only one shard remains resident"
        );

        expect(
            engine.cache().contains(
                "chrB"
            ),
            "Final chromosome remains resident"
        );

        expect(
            !engine.cache().contains(
                "chrA"
            ),
            "Earlier chromosome was evicted"
        );


        std::cout
            << "shards\t"
            << result.shards.size()
            << '\n';

        std::cout
            << "total_hits\t"
            << result.total_primer_hits()
            << '\n';

        std::cout
            << "window_candidates\t"
            << result.total_window_candidates()
            << '\n';

        std::cout
            << "cross_amplicons\t"
            << result.total_cross_amplicons()
            << '\n';

        std::cout
            << "cache_loads\t"
            << engine.cache().load_count()
            << '\n';

        std::cout
            << "cache_evictions\t"
            << engine.cache().eviction_count()
            << '\n';

        std::cout
            << "NO_CROSS_CHROMOSOME_PAIRING\tYES\n";

        std::cout
            << "LOCAL_ORACLE_EQUIVALENCE\tYES\n";

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
