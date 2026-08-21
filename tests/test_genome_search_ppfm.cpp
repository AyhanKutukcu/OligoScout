#include "primerpair/genome_search.hpp"
#include "primerpair/ppfm_io.hpp"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

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


std::string reverse_complement(
    const std::string& sequence
) {
    std::string result;

    result.reserve(
        sequence.size()
    );

    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {
        switch (*it) {

            case 'A':
                result.push_back('T');
                break;

            case 'C':
                result.push_back('G');
                break;

            case 'G':
                result.push_back('C');
                break;

            case 'T':
                result.push_back('A');
                break;

            default:
                throw std::runtime_error(
                    "Unexpected primer base."
                );
        }
    }

    return result;
}

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const std::string primer1 =
            "ACGTACGTTGCAACGTACGT";

        const std::string primer2 =
            "TTGGAACCTTGGAACCTTGG";

        const std::string right_site =
            reverse_complement(
                primer2
            );

        /*
         * Exact 100-bp product:
         *
         * start          = 50
         * primer1        = [50,70)
         * primer2 site   = [130,150)
         * amplicon       = [50,150)
         */
        const std::string sequence =
            std::string(
                50,
                'N'
            )
            +
            primer1
            +
            std::string(
                60,
                'N'
            )
            +
            right_site
            +
            std::string(
                50,
                'N'
            );

        const auto suffix =
            std::chrono::
                high_resolution_clock::
                now().
                time_since_epoch().
                count();

        const auto temp =
            std::filesystem::
                temp_directory_path();

        const auto path_a =
            temp /
            (
                "primerpair_ppfm_chrA_" +
                std::to_string(
                    suffix
                ) +
                ".ppfm"
            );

        const auto path_b =
            temp /
            (
                "primerpair_ppfm_chrB_" +
                std::to_string(
                    suffix
                ) +
                ".ppfm"
            );

        const auto path_rate4 =
            temp /
            (
                "primerpair_ppfm_chrC_" +
                std::to_string(
                    suffix
                ) +
                ".ppfm"
            );


        {
            PackedReference reference(
                sequence
            );

            BidirectionalFMIndex index(
                sequence,
                8
            );

            PpfmIO::save_shard(
                path_a,
                "chrA",
                reference,
                index
            );

            PpfmIO::save_shard(
                path_b,
                "chrB",
                reference,
                index
            );
        }


        GenomeSearchEngine genome(
            8
        );

        const std::size_t id_a =
            genome.load_ppfm_shard(
                path_a
            );

        const std::size_t id_b =
            genome.load_ppfm_shard(
                path_b
            );

        expect(
            id_a == 0,
            "First persistent shard gets id 0"
        );

        expect(
            id_b == 1,
            "Second persistent shard gets id 1"
        );

        expect(
            genome.shard(0).chromosome() ==
                "chrA",
            "First chromosome metadata restored"
        );

        expect(
            genome.shard(1).chromosome() ==
                "chrB",
            "Second chromosome metadata restored"
        );


        const auto result =
            genome.search_pair(
                primer1,
                primer2,
                3,
                100,
                100
            );

        expect(
            result.shards.size() == 2,
            "Two persistent shards searched"
        );

        expect(
            result.total_amplicon_count() ==
                2,
            "Persistent genome search finds two products"
        );


        for (
            const auto& shard_result :
            result.shards
        ) {
            const auto& hits =
                shard_result
                    .search_result
                    .pair_result
                    .amplicons;

            expect(
                hits.size() == 1,
                "Each persistent shard has one product"
            );

            const auto& hit =
                hits.front();

            expect(
                hit.left_position == 50,
                "Persistent left position correct"
            );

            expect(
                hit.right_position == 130,
                "Persistent right position correct"
            );

            expect(
                hit.amplicon_start == 50,
                "Persistent amplicon start correct"
            );

            expect(
                hit.amplicon_end_exclusive ==
                    150,
                "Persistent amplicon end correct"
            );

            expect(
                hit.amplicon_length == 100,
                "Persistent amplicon length correct"
            );

            expect(
                hit.total_mismatches() == 0,
                "Persistent amplicon is exact"
            );
        }


        bool duplicate_rejected = false;

        try {
            genome.load_ppfm_shard(
                path_a
            );

        } catch (
            const std::invalid_argument&
        ) {
            duplicate_rejected = true;
        }

        expect(
            duplicate_rejected,
            "Duplicate persistent chromosome rejected"
        );


        /*
         * SA-rate mismatch must not be accepted
         * silently.
         */
        {
            PackedReference reference(
                sequence
            );

            BidirectionalFMIndex index(
                sequence,
                4
            );

            PpfmIO::save_shard(
                path_rate4,
                "chrC",
                reference,
                index
            );
        }

        bool rate_rejected = false;

        try {
            GenomeSearchEngine genome8(
                8
            );

            genome8.load_ppfm_shard(
                path_rate4
            );

        } catch (
            const std::invalid_argument&
        ) {
            rate_rejected = true;
        }

        expect(
            rate_rejected,
            "Persistent SA-rate mismatch rejected"
        );


        std::error_code error;

        std::filesystem::remove(
            path_a,
            error
        );

        std::filesystem::remove(
            path_b,
            error
        );

        std::filesystem::remove(
            path_rate4,
            error
        );


        std::cout
            << "Genome PPFM loader tests passed.\n";

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
