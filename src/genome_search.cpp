#include "primerpair/genome_search.hpp"

#include "primerpair/ppfm_io.hpp"
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

namespace primerpair {

namespace {


std::string normalize_sequence(
    std::string sequence
) {
    if (sequence.empty()) {
        throw std::invalid_argument(
            "Genome shard sequence cannot be empty."
        );
    }

    for (char& raw : sequence) {

        const char base =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        raw
                    )
                )
            );

        switch (base) {

            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
                raw = base;
                break;

            /*
             * Keep the same ambiguity convention
             * already used by the real-FASTA
             * benchmarks: IUPAC ambiguity -> N.
             */
            case 'R':
            case 'Y':
            case 'S':
            case 'W':
            case 'K':
            case 'M':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
                raw = 'N';
                break;

            default:
                throw std::invalid_argument(
                    "Genome shard contains an "
                    "unsupported nucleotide."
                );
        }
    }

    return sequence;
}


}  // namespace


GenomeShard::GenomeShard(
    const std::size_t id,
    std::string chromosome,
    std::string sequence,
    const std::size_t suffix_array_sample_rate
)
    : id_(id),
      chromosome_(
          std::move(
              chromosome
          )
      ),
      sequence_(
          normalize_sequence(
              std::move(
                  sequence
              )
          )
      ),
      reference_(
          sequence_
      ),
      index_(
          sequence_,
          suffix_array_sample_rate
      ),
      pair_searcher_(
          index_,
          reference_
      ) {

    if (chromosome_.empty()) {
        throw std::invalid_argument(
            "Genome shard chromosome name "
            "cannot be empty."
        );
    }

    /*
     * PackedReference and both FM indexes have now
     * been fully constructed from sequence_.
     *
     * Retaining the original chromosome string
     * would cost roughly one additional byte/base
     * across the genome for no search-time benefit.
     */
    sequence_length_ =
        static_cast<std::uint64_t>(
            sequence_.size()
        );

    std::string empty;
    sequence_.swap(
        empty
    );
}



GenomeShard::GenomeShard(
    const std::size_t id,
    std::string chromosome,
    PackedReference reference,
    BidirectionalFMIndex index
)
    : id_(id),
      chromosome_(
          std::move(
              chromosome
          )
      ),
      sequence_(),
      reference_(
          std::move(
              reference
          )
      ),
      index_(
          std::move(
              index
          )
      ),
      pair_searcher_(
          index_,
          reference_
      ) {

    if (chromosome_.empty()) {
        throw std::invalid_argument(
            "Genome shard chromosome name "
            "cannot be empty."
        );
    }

    sequence_length_ =
        reference_.size();

    if (sequence_length_ == 0) {
        throw std::invalid_argument(
            "Persistent genome shard reference "
            "cannot be empty."
        );
    }
}


SensitivePairConstrainedSearchResult
GenomeShard::search_pair(
    const std::string_view primer1,
    const std::string_view primer2,
    const std::size_t max_mismatches,
    const std::uint64_t min_amplicon_length,
    const std::uint64_t max_amplicon_length,
    const SensitivePairAnchorPolicy anchor_policy
) const {
    return pair_searcher_.search(
        primer1,
        primer2,
        max_mismatches,
        min_amplicon_length,
        max_amplicon_length,
        anchor_policy
    );
}


std::uint64_t
GenomePairSearchResult::
total_amplicon_count() const noexcept {
    std::uint64_t total = 0;

    for (const auto& shard : shards) {
        total +=
            static_cast<std::uint64_t>(
                shard
                    .search_result
                    .pair_result
                    .amplicons
                    .size()
            );
    }

    return total;
}


GenomeSearchEngine::GenomeSearchEngine(
    const std::size_t suffix_array_sample_rate
)
    : suffix_array_sample_rate_(
          suffix_array_sample_rate
      ) {

    if (
        suffix_array_sample_rate_ == 0
    ) {
        throw std::invalid_argument(
            "Genome shard suffix-array "
            "sample rate must be > 0."
        );
    }
}


std::size_t
GenomeSearchEngine::add_shard(
    std::string chromosome,
    std::string sequence
) {
    if (chromosome.empty()) {
        throw std::invalid_argument(
            "Genome shard chromosome name "
            "cannot be empty."
        );
    }

    if (sequence.empty()) {
        throw std::invalid_argument(
            "Genome shard sequence cannot be empty."
        );
    }

    for (const auto& existing : shards_) {

        if (
            existing->chromosome() ==
            chromosome
        ) {
            throw std::invalid_argument(
                "Duplicate genome shard "
                "chromosome name: " +
                chromosome
            );
        }
    }


    const std::size_t id =
        shards_.size();

    shards_.push_back(
        std::make_unique<GenomeShard>(
            id,
            std::move(
                chromosome
            ),
            std::move(
                sequence
            ),
            suffix_array_sample_rate_
        )
    );

    return id;
}



std::size_t
GenomeSearchEngine::load_ppfm_shard(
    const std::filesystem::path& path
) {
    auto loaded =
        PpfmIO::load_shard(
            path
        );

    const std::size_t forward_rate =
        loaded.index
            .forward_index()
            .suffix_array_sample_rate();

    const std::size_t reverse_rate =
        loaded.index
            .reverse_index()
            .suffix_array_sample_rate();

    if (
        forward_rate !=
            suffix_array_sample_rate_
        ||
        reverse_rate !=
            suffix_array_sample_rate_
    ) {
        throw std::invalid_argument(
            "Persistent shard SA sample rate "
            "does not match GenomeSearchEngine."
        );
    }

    for (const auto& existing : shards_) {

        if (
            existing->chromosome() ==
            loaded.chromosome
        ) {
            throw std::invalid_argument(
                "Duplicate genome shard "
                "chromosome name: " +
                loaded.chromosome
            );
        }
    }

    const std::size_t id =
        shards_.size();

    shards_.push_back(
        std::make_unique<GenomeShard>(
            id,
            std::move(
                loaded.chromosome
            ),
            std::move(
                loaded.reference
            ),
            std::move(
                loaded.index
            )
        )
    );

    return id;
}


const GenomeShard&
GenomeSearchEngine::shard(
    const std::size_t shard_id
) const {
    if (
        shard_id >=
        shards_.size()
    ) {
        throw std::out_of_range(
            "Genome shard id out of range."
        );
    }

    return *shards_.at(
        shard_id
    );
}


GenomePairSearchResult
GenomeSearchEngine::search_pair(
    const std::string_view primer1,
    const std::string_view primer2,
    const std::size_t max_mismatches,
    const std::uint64_t min_amplicon_length,
    const std::uint64_t max_amplicon_length,
    const SensitivePairAnchorPolicy anchor_policy
) const {
    GenomePairSearchResult result;

    result.shards.reserve(
        shards_.size()
    );

    /*
     * Deliberately search+pair each chromosome
     * independently.
     *
     * Never collect primer hits globally and pair
     * afterwards: that could create biologically
     * impossible cross-chromosome amplicons.
     */
    for (const auto& shard_ptr : shards_) {

        auto shard_result =
            shard_ptr->search_pair(
                primer1,
                primer2,
                max_mismatches,
                min_amplicon_length,
                max_amplicon_length,
                anchor_policy
            );

        result.shards.push_back(
            GenomeShardPairSearchResult{
                shard_ptr->id(),
                shard_ptr->chromosome(),
                std::move(
                    shard_result
                )
            }
        );
    }

    return result;
}

}  // namespace primerpair
