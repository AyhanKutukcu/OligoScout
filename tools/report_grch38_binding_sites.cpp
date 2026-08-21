#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"
#include "primerpair/primer_binding_site.hpp"
#include "primerpair/primer_pair_search.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

int main(
    const int argc,
    char** argv
) {
    using namespace primerpair;

    try {
        if (argc != 5) {
            std::cerr
                << "Usage:\n"
                << "  report_grch38_binding_sites "
                << "<manifest> "
                << "<index_dir> "
                << "<primer1> "
                << "<primer2>\n";

            return 2;
        }

        const std::string manifest_path =
            argv[1];

        const std::string index_dir =
            argv[2];

        const std::string primer1 =
            argv[3];

        const std::string primer2 =
            argv[4];


        auto manifest =
            PpfmManifest::load(
                manifest_path,
                index_dir
            );

        PpfmShardCache cache(
            manifest,
            2,
            8
        );

        std::uint64_t total_candidates = 0;


        for (
            const auto& entry :
            manifest.entries()
        ) {
            const GenomeShard& shard =
                cache.get(
                    entry.chromosome
                );

            const auto search =
                shard.search_pair(
                    primer1,
                    primer2,
                    3,
                    50,
                    3000
                );

            const auto& hits =
                search
                    .pair_result
                    .amplicons;


            for (
                std::size_t i = 0;
                i < hits.size();
                ++i
            ) {
                const auto& hit =
                    hits.at(i);

                const auto sites =
                    extract_pair_binding_sites(
                        shard.reference(),
                        hit,
                        primer1,
                        primer2
                    );


                std::cout
                    << "CANDIDATE\t"
                    << shard.id()
                    << '\t'
                    << shard.chromosome()
                    << '\t'
                    << i
                    << '\t'
                    << hit.amplicon_start
                    << '\t'
                    << hit.amplicon_end_exclusive
                    << '\t'
                    << hit.amplicon_length
                    << '\t'
                    << hit.total_mismatches()
                    << '\n';


                std::cout
                    << "BINDING\t"
                    << shard.id()
                    << '\t'
                    << shard.chromosome()
                    << '\t'
                    << i
                    << '\t'
                    << "LEFT"
                    << '\t'
                    << to_string(
                        sites.left.primer
                    )
                    << '\t'
                    << (
                        sites.left.reverse_strand
                        ? "REVERSE"
                        : "FORWARD"
                    )
                    << '\t'
                    << sites.left.genomic_start
                    << '\t'
                    << sites.left.genomic_end_exclusive
                    << '\t'
                    << sites.left.primer_sequence
                    << '\t'
                    << sites.left.genomic_sequence
                    << '\t'
                    << sites.left.binding_sequence
                    << '\t'
                    << sites.left.mismatch_count()
                    << '\t'
                    << sites.left.mismatch_mask
                    << '\n';


                std::cout
                    << "BINDING\t"
                    << shard.id()
                    << '\t'
                    << shard.chromosome()
                    << '\t'
                    << i
                    << '\t'
                    << "RIGHT"
                    << '\t'
                    << to_string(
                        sites.right.primer
                    )
                    << '\t'
                    << (
                        sites.right.reverse_strand
                        ? "REVERSE"
                        : "FORWARD"
                    )
                    << '\t'
                    << sites.right.genomic_start
                    << '\t'
                    << sites.right.genomic_end_exclusive
                    << '\t'
                    << sites.right.primer_sequence
                    << '\t'
                    << sites.right.genomic_sequence
                    << '\t'
                    << sites.right.binding_sequence
                    << '\t'
                    << sites.right.mismatch_count()
                    << '\t'
                    << sites.right.mismatch_mask
                    << '\n';


                for (
                    const auto& mismatch :
                    sites.left.mismatches
                ) {
                    std::cout
                        << "MISMATCH\t"
                        << shard.id()
                        << '\t'
                        << shard.chromosome()
                        << '\t'
                        << i
                        << '\t'
                        << "LEFT"
                        << '\t'
                        << to_string(
                            sites.left.primer
                        )
                        << '\t'
                        << mismatch.primer_position
                        << '\t'
                        << (
                            mismatch.primer_position +
                            1
                        )
                        << '\t'
                        << mismatch.distance_from_3prime
                        << '\t'
                        << mismatch.primer_base
                        << '\t'
                        << mismatch.target_base
                        << '\n';
                }


                for (
                    const auto& mismatch :
                    sites.right.mismatches
                ) {
                    std::cout
                        << "MISMATCH\t"
                        << shard.id()
                        << '\t'
                        << shard.chromosome()
                        << '\t'
                        << i
                        << '\t'
                        << "RIGHT"
                        << '\t'
                        << to_string(
                            sites.right.primer
                        )
                        << '\t'
                        << mismatch.primer_position
                        << '\t'
                        << (
                            mismatch.primer_position +
                            1
                        )
                        << '\t'
                        << mismatch.distance_from_3prime
                        << '\t'
                        << mismatch.primer_base
                        << '\t'
                        << mismatch.target_base
                        << '\n';
                }


                ++total_candidates;
            }
        }


        std::cout
            << "total_candidates\t"
            << total_candidates
            << '\n';

        std::cout
            << "cache_loads\t"
            << cache.load_count()
            << '\n';

        std::cout
            << "cache_evictions\t"
            << cache.eviction_count()
            << '\n';

        if (total_candidates == 0) {
            throw std::runtime_error(
                "No genome-wide candidates found."
            );
        }

        std::cout
            << "ALL_CHECKS\tYES\n";

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
