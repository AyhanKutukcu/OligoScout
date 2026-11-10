#include "primerpair/global_multiplex_cross_join.hpp"
#include "primerpair/primer_pair_assembler.hpp"
#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct PanelPair {
    std::string primer1;
    std::string primer2;
};

struct CrossRecord {
    std::size_t forward_uid{0};
    std::size_t reverse_uid{0};

    std::uint64_t forward_position{0};
    std::uint64_t reverse_position{0};

    std::uint64_t amplicon_end_exclusive{0};
    std::uint64_t amplicon_length{0};

    std::size_t forward_mismatches{0};
    std::size_t reverse_mismatches{0};

    std::uint64_t forward_mismatch_mask{0};
    std::uint64_t reverse_mismatch_mask{0};

    bool operator==(
        const CrossRecord&
    ) const = default;
};


bool record_less(
    const CrossRecord& lhs,
    const CrossRecord& rhs
) noexcept {
    return std::tie(
        lhs.forward_uid,
        lhs.reverse_uid,
        lhs.forward_position,
        lhs.reverse_position,
        lhs.amplicon_end_exclusive,
        lhs.amplicon_length,
        lhs.forward_mismatches,
        lhs.reverse_mismatches,
        lhs.forward_mismatch_mask,
        lhs.reverse_mismatch_mask
    )
    <
    std::tie(
        rhs.forward_uid,
        rhs.reverse_uid,
        rhs.forward_position,
        rhs.reverse_position,
        rhs.amplicon_end_exclusive,
        rhs.amplicon_length,
        rhs.forward_mismatches,
        rhs.reverse_mismatches,
        rhs.forward_mismatch_mask,
        rhs.reverse_mismatch_mask
    );
}


void normalize(
    std::vector<CrossRecord>& records
) {
    std::sort(
        records.begin(),
        records.end(),
        record_less
    );

    records.erase(
        std::unique(
            records.begin(),
            records.end()
        ),
        records.end()
    );
}


double elapsed_ms(
    const Clock::time_point begin,
    const Clock::time_point end
) {
    return std::chrono::duration<
        double,
        std::milli
    >(
        end - begin
    ).count();
}


std::vector<PanelPair>
load_panel(
    const std::string& path,
    const std::size_t pair_limit
) {
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Cannot open panel: " + path
        );
    }

    std::string line;

    if (!std::getline(input, line)) {
        throw std::runtime_error(
            "Panel is empty."
        );
    }

    std::vector<PanelPair> panel;

    while (std::getline(input, line)) {

        if (line.empty()) {
            continue;
        }

        std::istringstream row(line);

        std::size_t pair_id = 0;
        std::string primer1;
        std::string primer2;

        if (!(row >> pair_id >> primer1 >> primer2)) {
            throw std::runtime_error(
                "Malformed panel row."
            );
        }

        panel.push_back(
            PanelPair{
                primer1,
                primer2
            }
        );

        if (
            pair_limit > 0 &&
            panel.size() >= pair_limit
        ) {
            break;
        }
    }

    if (panel.empty()) {
        throw std::runtime_error(
            "No primer pairs loaded."
        );
    }

    return panel;
}


void append_directed_baseline(
    const std::size_t forward_uid,
    const std::size_t reverse_uid,

    const std::vector<
        primerpair::OrientedPrimerSearchHit
    >& forward_hits,

    const std::vector<
        primerpair::OrientedPrimerSearchHit
    >& reverse_hits,

    const std::size_t forward_length,
    const std::size_t reverse_length,

    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon,

    std::vector<CrossRecord>& output
) {
    using namespace primerpair;

    for (const auto& forward : forward_hits) {

        if (
            forward.orientation !=
            PrimerOrientation::Forward
        ) {
            continue;
        }

        if (
            forward.position >
            std::numeric_limits<
                std::uint64_t
            >::max() - forward_length
        ) {
            continue;
        }

        const std::uint64_t
            minimum_reverse_start =
                forward.position +
                static_cast<std::uint64_t>(
                    forward_length
                );

        for (const auto& reverse : reverse_hits) {

            if (
                reverse.orientation !=
                PrimerOrientation::Reverse
            ) {
                continue;
            }

            if (
                reverse.position <
                minimum_reverse_start
            ) {
                continue;
            }

            if (
                reverse.position >
                std::numeric_limits<
                    std::uint64_t
                >::max() - reverse_length
            ) {
                continue;
            }

            const std::uint64_t reverse_end =
                reverse.position +
                static_cast<std::uint64_t>(
                    reverse_length
                );

            if (
                reverse_end <=
                forward.position
            ) {
                continue;
            }

            const std::uint64_t amplicon_length =
                reverse_end -
                forward.position;

            if (
                amplicon_length <
                    min_amplicon ||
                amplicon_length >
                    max_amplicon
            ) {
                continue;
            }

            output.push_back(
                CrossRecord{
                    forward_uid,
                    reverse_uid,

                    forward.position,
                    reverse.position,

                    reverse_end,
                    amplicon_length,

                    forward.mismatches,
                    reverse.mismatches,

                    forward.mismatch_mask,
                    reverse.mismatch_mask
                }
            );
        }
    }
}


std::vector<CrossRecord>
legacy_logical_join(
    const std::vector<std::string>& primers,

    const std::vector<
        std::vector<
            primerpair::OrientedPrimerSearchHit
        >
    >& hits,

    const std::vector<std::size_t>& lengths,

    const std::size_t pair_count,

    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon
) {
    using namespace primerpair;

    std::vector<CrossRecord> output;


    /*
     * Production V1 logical multiplex baseline:
     *
     *   pair_a < pair_b
     *       four slot combinations
     *
     * Each logical primer combination uses the same
     * monotonic assembler as production Multiplex V1:
     *
     *   assemble_primer_pair_hits()
     *
     * The GRCh38 benchmark panel contains unique
     * primer sequences, so V1's duplicate join cache
     * does not collapse these logical requests.
     */
    for (
        std::size_t pair_a = 0;
        pair_a < pair_count;
        ++pair_a
    ) {
        for (
            std::size_t pair_b =
                pair_a + 1;

            pair_b < pair_count;
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
                    const std::size_t uid_a =
                        pair_a * 2 +
                        slot_a;

                    const std::size_t uid_b =
                        pair_b * 2 +
                        slot_b;


                    const auto assembled =
                        assemble_primer_pair_hits(
                            primers.at(uid_a),
                            hits.at(uid_a),

                            primers.at(uid_b),
                            hits.at(uid_b),

                            min_amplicon,
                            max_amplicon
                        );


                    for (
                        const auto& product :
                        assembled.amplicons
                    ) {
                        std::size_t forward_uid = 0;
                        std::size_t reverse_uid = 0;


                        if (
                            product.left_primer ==
                                PrimerIdentity::Primer1
                            &&
                            product.right_primer ==
                                PrimerIdentity::Primer2
                        ) {
                            forward_uid = uid_a;
                            reverse_uid = uid_b;

                        } else if (
                            product.left_primer ==
                                PrimerIdentity::Primer2
                            &&
                            product.right_primer ==
                                PrimerIdentity::Primer1
                        ) {
                            forward_uid = uid_b;
                            reverse_uid = uid_a;

                        } else {
                            throw std::runtime_error(
                                "Unexpected primer identity "
                                "from production assembler."
                            );
                        }


                        const std::uint64_t reverse_length =
                            static_cast<std::uint64_t>(
                                lengths.at(
                                    reverse_uid
                                )
                            );


                        if (
                            product.amplicon_end_exclusive
                            <
                            reverse_length
                        ) {
                            throw std::runtime_error(
                                "Invalid production V1 "
                                "amplicon coordinates."
                            );
                        }


                        const std::uint64_t reverse_position =
                            product.amplicon_end_exclusive
                            -
                            reverse_length;


                        output.push_back(
                            CrossRecord{
                                forward_uid,
                                reverse_uid,

                                product.amplicon_start,
                                reverse_position,

                                product.amplicon_end_exclusive,
                                product.amplicon_length,

                                product.left_mismatches,
                                product.right_mismatches,

                                product.left_mismatch_mask,
                                product.right_mismatch_mask
                            }
                        );
                    }
                }
            }
        }
    }


    normalize(
        output
    );


    return output;
}

std::vector<CrossRecord>
global_join(
    const std::vector<
        primerpair::GlobalMultiplexPrimerHits
    >& inputs,

    const std::vector<std::size_t>&
        uid_to_pair,

    const std::uint64_t min_amplicon,
    const std::uint64_t max_amplicon,

    primerpair::
        GlobalMultiplexCrossJoinStats*
            output_stats
) {
    using namespace primerpair;

    const auto result =
        global_multiplex_cross_join(
            inputs,
            min_amplicon,
            max_amplicon
        );

    if (output_stats != nullptr) {
        *output_stats =
            result.stats;
    }

    std::vector<CrossRecord> output;

    output.reserve(
        result.products.size()
    );

    for (const auto& product : result.products) {

        /*
         * Cross-pair benchmark:
         * intended pair içindeki ürünleri çıkar.
         */
        if (
            uid_to_pair.at(
                product.forward_primer_id
            )
            ==
            uid_to_pair.at(
                product.reverse_primer_id
            )
        ) {
            continue;
        }

        output.push_back(
            CrossRecord{
                product.forward_primer_id,
                product.reverse_primer_id,

                product.forward_position,
                product.reverse_position,

                product.amplicon_end_exclusive,
                product.amplicon_length,

                product.forward_mismatches,
                product.reverse_mismatches,

                product.forward_mismatch_mask,
                product.reverse_mismatch_mask
            }
        );
    }

    normalize(output);

    return output;
}

}  // namespace


int main(
    const int argc,
    char** argv
) {
    try {
        using namespace primerpair;

        if (
            argc < 5 ||
            argc > 7
        ) {
            std::cerr
                << "Usage:\n"
                << "benchmark_grch38_shard_multiplex_v2 "
                << "<manifest.tsv> "
                << "<index_dir> "
                << "<chromosome> "
                << "<panel.tsv> "
                << "[pair_limit=8] "
                << "[repetitions=5]\n";

            return 2;
        }

        const std::string manifest_path =
            argv[1];

        const std::string index_dir =
            argv[2];

        const std::string chromosome =
            argv[3];

        const std::string panel_path =
            argv[4];

        const std::size_t pair_limit =
            argc >= 6
                ? static_cast<std::size_t>(
                      std::stoull(argv[5])
                  )
                : 8;

        const std::size_t repetitions =
            argc >= 7
                ? static_cast<std::size_t>(
                      std::stoull(argv[6])
                  )
                : 5;

        if (repetitions == 0) {
            throw std::invalid_argument(
                "Repetitions must be > 0."
            );
        }

        constexpr std::size_t anchor_length = 12;
        constexpr std::size_t max_mismatches = 3;

        constexpr std::uint64_t min_amplicon = 50;
        constexpr std::uint64_t max_amplicon = 3000;


        const auto panel =
            load_panel(
                panel_path,
                pair_limit
            );


        std::set<std::string>
            unique_primers;

        for (const auto& pair : panel) {
            unique_primers.insert(
                pair.primer1
            );

            unique_primers.insert(
                pair.primer2
            );
        }

        if (
            unique_primers.size() !=
            panel.size() * 2
        ) {
            throw std::runtime_error(
                "Pilot panel is not all-unique."
            );
        }


        const auto manifest =
            PpfmManifest::load(
                manifest_path,
                index_dir
            );


        PpfmShardCache cache(
            manifest,
            1,
            8
        );


        const auto load_begin =
            Clock::now();

        const GenomeShard& shard =
            cache.get(chromosome);

        const auto load_end =
            Clock::now();


        StrandAwarePrimerSearchEngine engine(
            shard.index(),
            shard.reference()
        );


        std::vector<std::string> primers;
        std::vector<std::size_t> uid_to_pair;

        primers.reserve(
            panel.size() * 2
        );

        uid_to_pair.reserve(
            panel.size() * 2
        );


        for (
            std::size_t pair_index = 0;
            pair_index < panel.size();
            ++pair_index
        ) {
            primers.push_back(
                panel.at(pair_index).primer1
            );

            uid_to_pair.push_back(
                pair_index
            );

            primers.push_back(
                panel.at(pair_index).primer2
            );

            uid_to_pair.push_back(
                pair_index
            );
        }


        std::vector<
            std::vector<
                OrientedPrimerSearchHit
            >
        > hit_lists;

        hit_lists.reserve(
            primers.size()
        );


        std::size_t total_hits = 0;
        std::size_t forward_hits = 0;
        std::size_t reverse_hits = 0;


        const auto search_begin =
            Clock::now();


        for (const auto& primer : primers) {

            const auto search =
                engine.search(
                    primer,
                    anchor_length,
                    max_mismatches
                );

            for (const auto& hit : search.hits) {

                ++total_hits;

                if (
                    hit.orientation ==
                    PrimerOrientation::Forward
                ) {
                    ++forward_hits;
                } else {
                    ++reverse_hits;
                }
            }

            hit_lists.push_back(
                search.hits
            );
        }


        const auto search_end =
            Clock::now();


        std::vector<std::size_t>
            primer_lengths;

        primer_lengths.reserve(
            primers.size()
        );


        std::vector<
            GlobalMultiplexPrimerHits
        > global_inputs;

        global_inputs.reserve(
            primers.size()
        );


        for (
            std::size_t uid = 0;
            uid < primers.size();
            ++uid
        ) {
            primer_lengths.push_back(
                primers.at(uid).size()
            );

            global_inputs.push_back(
                GlobalMultiplexPrimerHits{
                    uid,
                    primers.at(uid).size(),
                    hit_lists.at(uid)
                }
            );
        }


        /*
         * Correctness check.
         */
        const auto baseline =
            legacy_logical_join(
                primers,
                hit_lists,
                primer_lengths,
                panel.size(),
                min_amplicon,
                max_amplicon
            );


        GlobalMultiplexCrossJoinStats
            global_stats;


        const auto observed =
            global_join(
                global_inputs,
                uid_to_pair,
                min_amplicon,
                max_amplicon,
                &global_stats
            );


        if (
            baseline !=
            observed
        ) {
            std::cerr
                << "SHARD_EQUIVALENCE_MISMATCH"
                << '\t'
                << "baseline="
                << baseline.size()
                << '\t'
                << "global="
                << observed.size()
                << '\n';

            return 1;
        }


        std::uint64_t sink = 0;


        const auto baseline_begin =
            Clock::now();


        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            const auto result =
                legacy_logical_join(
                primers,
                hit_lists,
                    primer_lengths,
                    panel.size(),
                    min_amplicon,
                    max_amplicon
                );

            sink += result.size();
        }


        const auto baseline_end =
            Clock::now();


        const auto global_begin =
            Clock::now();


        for (
            std::size_t rep = 0;
            rep < repetitions;
            ++rep
        ) {
            const auto result =
                global_join(
                    global_inputs,
                    uid_to_pair,
                    min_amplicon,
                    max_amplicon,
                    nullptr
                );

            sink += result.size();
        }


        const auto global_end =
            Clock::now();


        const double baseline_ms =
            elapsed_ms(
                baseline_begin,
                baseline_end
            )
            /
            repetitions;


        const double global_ms =
            elapsed_ms(
                global_begin,
                global_end
            )
            /
            repetitions;


        const std::size_t
            logical_cross_requests =
                (
                    panel.size() *
                    (
                        panel.size() - 1
                    )
                    /
                    2
                )
                *
                4;


        std::cout
            << "chromosome\t"
            << shard.chromosome()
            << '\n';

        std::cout
            << "sequence_length\t"
            << shard.sequence_length()
            << '\n';

        std::cout
            << "sa_rate\t"
            << shard.suffix_array_sample_rate()
            << '\n';

        std::cout
            << "pair_count\t"
            << panel.size()
            << '\n';

        std::cout
            << "primer_count\t"
            << primers.size()
            << '\n';

        std::cout
            << "unique_primers\t"
            << unique_primers.size()
            << '\n';

        std::cout
            << "total_hits\t"
            << total_hits
            << '\n';

        std::cout
            << "forward_hits\t"
            << forward_hits
            << '\n';

        std::cout
            << "reverse_hits\t"
            << reverse_hits
            << '\n';

        std::cout
            << "logical_cross_requests\t"
            << logical_cross_requests
            << '\n';

        std::cout
            << "global_window_candidates\t"
            << global_stats.window_candidates
            << '\n';

        std::cout
            << "global_unique_products\t"
            << global_stats.unique_products
            << '\n';

        std::cout
            << "cross_products\t"
            << observed.size()
            << '\n';

        std::cout
            << "ppfm_load_ms\t"
            << elapsed_ms(
                load_begin,
                load_end
            )
            << '\n';

        std::cout
            << "primer_search_ms\t"
            << elapsed_ms(
                search_begin,
                search_end
            )
            << '\n';

        std::cout
            << "baseline_join_ms\t"
            << baseline_ms
            << '\n';

        std::cout
            << "global_join_ms\t"
            << global_ms
            << '\n';

        std::cout
            << "join_speedup\t"
            << (
                baseline_ms /
                global_ms
            )
            << '\n';

        std::cout
            << "cache_loads\t"
            << cache.load_count()
            << '\n';

        std::cout
            << "cache_evictions\t"
            << cache.eviction_count()
            << '\n';

        std::cout
            << "VERIFY_SHARD_EQUIVALENT\tYES\n";


        std::cerr
            << "sink\t"
            << sink
            << '\n';


        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}
