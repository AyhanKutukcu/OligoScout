#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace primerpair {


struct FastaRecord {
    /*
     * First whitespace-delimited token after '>'.
     */
    std::string name;

    /*
     * Complete header text after '>'.
     */
    std::string description;

    /*
     * Uppercase normalized sequence.
     *
     * A/C/G/T/N preserved.
     * IUPAC ambiguity symbols converted to N.
     */
    std::string sequence;
};


using FastaRecordSelector =
    std::function<
        bool(
            std::string_view name,
            std::string_view description
        )
    >;


using FastaRecordConsumer =
    std::function<
        void(
            FastaRecord&& record
        )
    >;


/*
 * Stream a FASTA file and materialize sequence
 * only for records accepted by selector.
 *
 * Important for multi-gigabyte genome FASTA files:
 * unselected records are never accumulated in RAM.
 *
 * Returns the number of selected records delivered
 * to consumer.
 */
[[nodiscard]]
std::size_t
stream_selected_fasta_records(
    const std::string& path,
    const FastaRecordSelector& selector,
    const FastaRecordConsumer& consumer
);


/*
 * Convenience API for small FASTA files.
 *
 * This intentionally materializes every record.
 * Do not use it to load the complete GRCh38 FASTA
 * for whole-genome indexing.
 */
[[nodiscard]]
std::vector<FastaRecord>
load_fasta_records(
    const std::string& path
);


}  // namespace primerpair
