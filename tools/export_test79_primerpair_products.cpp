#include "primerpair/ppfm_manifest.hpp"
#include "primerpair/ppfm_shard_cache.hpp"
#include "primerpair/strand_aware_primer_search.hpp"
#include "primerpair/primer_pair_assembler.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>


namespace {


struct PanelPair {
    std::string pair_id;
    std::string primer1;
    std::string primer2;
};


std::vector<std::string>
split_tab(
    const std::string& line
) {
    std::vector<std::string> fields;

    std::size_t begin = 0;

    while (true) {

        const std::size_t end =
            line.find(
                '\t',
                begin
            );

        if (
            end ==
            std::string::npos
        ) {

            fields.emplace_back(
                line.substr(
                    begin
                )
            );

            break;
        }

        fields.emplace_back(
            line.substr(
                begin,
                end - begin
            )
        );

        begin =
            end + 1;
    }

    return fields;
}


std::vector<PanelPair>
read_panel(
    const std::string& path
) {
    std::ifstream input{
        path
    };

    if (!input) {
        throw std::runtime_error(
            "Cannot open panel TSV."
        );
    }

    std::string line;

    if (!std::getline(
        input,
        line
    )) {
        throw std::runtime_error(
            "Panel TSV is empty."
        );
    }

    const auto header =
        split_tab(
            line
        );

    if (
        header.size() < 3
        ||
        header.at(0) != "pair_id"
        ||
        header.at(1) != "primer1"
        ||
        header.at(2) != "primer2"
    ) {
        throw std::runtime_error(
            "Unexpected panel TSV header."
        );
    }

    std::vector<PanelPair> rows;

    while (
        std::getline(
            input,
            line
        )
    ) {

        if (line.empty()) {
            continue;
        }

        const auto fields =
            split_tab(
                line
            );

        if (fields.size() < 3) {
            throw std::runtime_error(
                "Malformed panel TSV row."
            );
        }

        rows.push_back(
            PanelPair{
                fields.at(0),
                fields.at(1),
                fields.at(2)
            }
        );
    }

    if (
        rows.size() !=
        32
    ) {
        throw std::runtime_error(
            "Expected exactly 32 panel pairs."
        );
    }

    return rows;
}


std::string
mask_hex(
    const std::uint64_t value
) {
    std::ostringstream stream;

    stream
        << "0x"
        << std::hex
        << std::nouppercase
        << std::setw(16)
        << std::setfill('0')
        << value;

    return stream.str();
}


const char*
primer_identity_name(
    const primerpair::PrimerIdentity identity
) {
    using primerpair::PrimerIdentity;

    switch (identity) {

        case PrimerIdentity::Primer1:
            return "P1";

        case PrimerIdentity::Primer2:
            return "P2";
    }

    return "UNKNOWN";
}


}  // namespace


int
main(
    const int argc,
    char** argv
) {
    try {

        if (argc != 5) {

            std::cerr
                << "Usage: "
                << argv[0]
                << " <manifest.tsv>"
                << " <ppfm_directory>"
                << " <panel.tsv>"
                << " <output.tsv>\n";

            return 2;
        }

        const std::string
            manifest_path{
                argv[1]
            };

        const std::string
            ppfm_directory{
                argv[2]
            };

        const std::string
            panel_path{
                argv[3]
            };

        const std::string
            output_path{
                argv[4]
            };


        auto panel =
            read_panel(
                panel_path
            );


        auto manifest =
            primerpair::PpfmManifest::load(
                manifest_path,
                ppfm_directory
            );


        if (
            manifest.size() !=
            24
        ) {
            throw std::runtime_error(
                "Expected 24 PPFM shards."
            );
        }


        std::vector<std::string>
            chromosomes;

        chromosomes.reserve(
            manifest.size()
        );


        for (
            const auto& entry :
            manifest.entries()
        ) {

            chromosomes.push_back(
                entry.chromosome
            );
        }


        primerpair::PpfmShardCache
            cache{
                std::move(
                    manifest
                ),
                1,
                8
            };


        std::ofstream output{
            output_path
        };


        if (!output) {
            throw std::runtime_error(
                "Cannot open output TSV."
            );
        }


        output
            << "pair_id\t"
            << "chromosome\t"
            << "left_primer\t"
            << "right_primer\t"
            << "left_position\t"
            << "right_position\t"
            << "amplicon_start\t"
            << "amplicon_end_exclusive\t"
            << "amplicon_length\t"
            << "left_mismatches\t"
            << "right_mismatches\t"
            << "total_mismatches\t"
            << "left_mismatch_mask\t"
            << "right_mismatch_mask\n";


        std::uint64_t
            search_instances = 0;

        std::uint64_t
            product_count = 0;


        for (
            const auto& chromosome :
            chromosomes
        ) {

            const auto& shard =
                cache.get(
                    chromosome
                );


            primerpair::StrandAwarePrimerSearchEngine
                searcher{
                    shard.index(),
                    shard.reference()
                };


            for (
                const auto& pair :
                panel
            ) {

                const auto p1 =
                    searcher.search(
                        pair.primer1,
                        12,
                        3
                    );

                const auto p2 =
                    searcher.search(
                        pair.primer2,
                        12,
                        3
                    );


                search_instances += 2;


                const auto assembled =
                    primerpair::
                    assemble_primer_pair_hits(
                        pair.primer1,
                        p1.hits,
                        pair.primer2,
                        p2.hits,
                        50,
                        3000
                    );


                for (
                    const auto& hit :
                    assembled.amplicons
                ) {

                    output
                        << pair.pair_id
                        << '\t'
                        << chromosome
                        << '\t'
                        << primer_identity_name(
                            hit.left_primer
                        )
                        << '\t'
                        << primer_identity_name(
                            hit.right_primer
                        )
                        << '\t'
                        << hit.left_position
                        << '\t'
                        << hit.right_position
                        << '\t'
                        << hit.amplicon_start
                        << '\t'
                        << hit.amplicon_end_exclusive
                        << '\t'
                        << hit.amplicon_length
                        << '\t'
                        << hit.left_mismatches
                        << '\t'
                        << hit.right_mismatches
                        << '\t'
                        << hit.total_mismatches()
                        << '\t'
                        << mask_hex(
                            hit.left_mismatch_mask
                        )
                        << '\t'
                        << mask_hex(
                            hit.right_mismatch_mask
                        )
                        << '\n';


                    ++product_count;
                }
            }
        }


        output.flush();


        if (!output) {
            throw std::runtime_error(
                "Failed while writing output TSV."
            );
        }


        std::cout
            << "panel_pairs\t"
            << panel.size()
            << '\n';

        std::cout
            << "ppfm_shards\t"
            << chromosomes.size()
            << '\n';

        std::cout
            << "primer_search_instances\t"
            << search_instances
            << '\n';

        std::cout
            << "product_count\t"
            << product_count
            << '\n';

        std::cout
            << "ppfm_loads\t"
            << cache.load_count()
            << '\n';

        std::cout
            << "ppfm_evictions\t"
            << cache.eviction_count()
            << '\n';

        std::cout
            << "search_anchor_length\t12\n";

        std::cout
            << "max_mismatches_per_primer\t3\n";

        std::cout
            << "amplicon_range\t50..3000\n";

        std::cout
            << "TEST79_PPS_PRODUCT_EXPORT\tYES\n";

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
