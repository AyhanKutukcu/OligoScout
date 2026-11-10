#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "primerpair/genome_search.hpp"

namespace {


void expect(
    const bool condition,
    const std::string& message
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: " +
            message
        );
    }

    std::cout
        << "[PASS] "
        << message
        << '\n';
}


char complement(
    const char base
) {
    switch (base) {
        case 'A':
            return 'T';

        case 'C':
            return 'G';

        case 'G':
            return 'C';

        case 'T':
            return 'A';

        default:
            throw std::invalid_argument(
                "Invalid DNA nucleotide."
            );
    }
}


std::string reverse_complement(
    const std::string_view sequence
) {
    std::string output;

    output.reserve(
        sequence.size()
    );

    for (
        auto it = sequence.rbegin();
        it != sequence.rend();
        ++it
    ) {
        output.push_back(
            complement(
                *it
            )
        );
    }

    return output;
}


/*
 * Reproducible high-complexity primer using only
 * C/G/T so it cannot approximately match the
 * A-only synthetic chromosome background.
 */
std::string make_primer(
    const std::size_t length,
    std::uint32_t state
) {
    static constexpr
        char bases[3]{
            'C',
            'G',
            'T'
        };

    std::string output;

    output.reserve(
        length
    );

    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;

        output.push_back(
            bases[
                static_cast<std::size_t>(
                    state % 3u
                )
            ]
        );
    }

    return output;
}


void insert_pair(
    std::string& chromosome,
    const std::string_view primer1,
    const std::string_view primer2,
    const std::size_t forward_start,
    const std::size_t reverse_start
) {
    chromosome.replace(
        forward_start,
        primer1.size(),
        primer1
    );

    const std::string reverse_target =
        reverse_complement(
            primer2
        );

    chromosome.replace(
        reverse_start,
        reverse_target.size(),
        reverse_target
    );
}


const primerpair::
GenomeShardPairSearchResult&
find_result(
    const primerpair::
        GenomePairSearchResult& result,
    const std::string_view chromosome
) {
    for (const auto& shard : result.shards) {

        if (
            shard.chromosome ==
            chromosome
        ) {
            return shard;
        }
    }

    throw std::runtime_error(
        "Chromosome result not found."
    );
}


std::size_t amplicon_count(
    const primerpair::
        GenomeShardPairSearchResult& result
) {
    return result
        .search_result
        .pair_result
        .amplicons
        .size();
}


}  // namespace


int main() {
    try {
        constexpr std::size_t
            primer_length = 30;

        const std::string
            primer_a1 =
                make_primer(
                    primer_length,
                    0x13579BDFu
                );

        const std::string
            primer_a2 =
                make_primer(
                    primer_length,
                    0x2468ACE1u
                );

        const std::string
            primer_b1 =
                make_primer(
                    primer_length,
                    0x10293847u
                );

        const std::string
            primer_b2 =
                make_primer(
                    primer_length,
                    0x56473829u
                );


        std::string chr_a(
            500,
            'A'
        );

        std::string chr_b(
            500,
            'A'
        );


        insert_pair(
            chr_a,
            primer_a1,
            primer_a2,
            60,
            180
        );

        insert_pair(
            chr_b,
            primer_b1,
            primer_b2,
            60,
            180
        );


        /*
         * SA8 is the current balanced whole-genome
         * candidate profile.
         */
        primerpair::GenomeSearchEngine
            genome(
                8
            );

        const std::size_t chr_a_id =
            genome.add_shard(
                "chrA",
                chr_a
            );

        const std::size_t chr_b_id =
            genome.add_shard(
                "chrB",
                chr_b
            );


        expect(
            chr_a_id == 0,
            "chrA receives shard id 0"
        );

        expect(
            chr_b_id == 1,
            "chrB receives shard id 1"
        );

        expect(
            genome.shard_count() == 2,
            "Genome contains two shards"
        );

        expect(
            genome.shard(0).chromosome() ==
                "chrA",
            "Shard 0 chromosome metadata"
        );

        expect(
            genome.shard(1).chromosome() ==
                "chrB",
            "Shard 1 chromosome metadata"
        );

        expect(
            genome.shard(0)
                .suffix_array_sample_rate()
                ==
                8,
            "Shard preserves SA8"
        );


        /*
         * --------------------------------------------------
         * chrA-specific pair
         * --------------------------------------------------
         */

        const auto result_a =
            genome.search_pair(
                primer_a1,
                primer_a2,
                3,
                50,
                300
            );

        expect(
            result_a.shard_count() == 2,
            "Pair search returns one result per shard"
        );

        const auto& a_on_a =
            find_result(
                result_a,
                "chrA"
            );

        const auto& a_on_b =
            find_result(
                result_a,
                "chrB"
            );

        expect(
            amplicon_count(
                a_on_a
            ) > 0,
            "chrA pair recovered on chrA"
        );

        expect(
            amplicon_count(
                a_on_b
            ) == 0,
            "chrA pair absent from chrB"
        );


        /*
         * --------------------------------------------------
         * chrB-specific pair
         * --------------------------------------------------
         */

        const auto result_b =
            genome.search_pair(
                primer_b1,
                primer_b2,
                3,
                50,
                300
            );

        expect(
            amplicon_count(
                find_result(
                    result_b,
                    "chrA"
                )
            ) == 0,
            "chrB pair absent from chrA"
        );

        expect(
            amplicon_count(
                find_result(
                    result_b,
                    "chrB"
                )
            ) > 0,
            "chrB pair recovered on chrB"
        );


        /*
         * --------------------------------------------------
         * Cross-chromosome false-pair test
         *
         * primer_a1 exists on chrA.
         * primer_b2 exists on chrB.
         *
         * A global hit-pool implementation could
         * accidentally pair these.
         *
         * Shard-local pairing must return zero.
         * --------------------------------------------------
         */

        const auto cross =
            genome.search_pair(
                primer_a1,
                primer_b2,
                3,
                50,
                300
            );

        expect(
            cross.total_amplicon_count() ==
                0,
            "Cross-chromosome hits never form "
            "an amplicon"
        );


        /*
         * --------------------------------------------------
         * Same pair, same coordinate, two chromosomes.
         *
         * Must remain two distinct chromosome-local
         * results rather than collapsing coordinates.
         * --------------------------------------------------
         */

        std::string chr_x(
            500,
            'A'
        );

        std::string chr_y(
            500,
            'A'
        );

        insert_pair(
            chr_x,
            primer_a1,
            primer_a2,
            80,
            200
        );

        insert_pair(
            chr_y,
            primer_a1,
            primer_a2,
            80,
            200
        );

        primerpair::GenomeSearchEngine
            duplicated_genome(
                8
            );

        duplicated_genome.add_shard(
            "chrX",
            chr_x
        );

        duplicated_genome.add_shard(
            "chrY",
            chr_y
        );

        const auto duplicate_result =
            duplicated_genome.search_pair(
                primer_a1,
                primer_a2,
                3,
                50,
                300
            );

        expect(
            amplicon_count(
                find_result(
                    duplicate_result,
                    "chrX"
                )
            ) > 0,
            "Same-coordinate product retained on chrX"
        );

        expect(
            amplicon_count(
                find_result(
                    duplicate_result,
                    "chrY"
                )
            ) > 0,
            "Same-coordinate product retained on chrY"
        );

        expect(
            duplicate_result
                .total_amplicon_count()
            >=
            2,
            "Identical coordinates on different "
            "chromosomes remain distinct"
        );


        /*
         * --------------------------------------------------
         * Metadata/input guards
         * --------------------------------------------------
         */

        bool duplicate_rejected =
            false;

        try {
            genome.add_shard(
                "chrA",
                std::string(
                    100,
                    'A'
                )
            );

        } catch (
            const std::invalid_argument&
        ) {
            duplicate_rejected =
                true;
        }

        expect(
            duplicate_rejected,
            "Duplicate chromosome name rejected"
        );


        bool zero_rate_rejected =
            false;

        try {
            primerpair::GenomeSearchEngine
                invalid(
                    0
                );

            static_cast<void>(
                invalid
            );

        } catch (
            const std::invalid_argument&
        ) {
            zero_rate_rejected =
                true;
        }

        expect(
            zero_rate_rejected,
            "Zero SA sampling rate rejected"
        );


        std::cout
            << "genome_shards\t"
            << genome.shard_count()
            << '\n';

        std::cout
            << "Genome shard isolation tests passed.\n";

        return 0;

    } catch (
        const std::exception& error
    ) {
        std::cerr
            << error.what()
            << '\n';

        return 1;
    }
}
