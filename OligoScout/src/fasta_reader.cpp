#include "primerpair/fasta_reader.hpp"

#include <cctype>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace primerpair {

namespace {


char normalize_base(
    const char raw
) {
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
            return base;

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
            return 'N';

        default:
            throw std::runtime_error(
                "Unsupported FASTA nucleotide."
            );
    }
}


std::string header_name(
    const std::string& description
) {
    std::size_t begin = 0;

    while (
        begin < description.size()
        &&
        std::isspace(
            static_cast<unsigned char>(
                description[begin]
            )
        )
    ) {
        ++begin;
    }

    std::size_t end = begin;

    while (
        end < description.size()
        &&
        !std::isspace(
            static_cast<unsigned char>(
                description[end]
            )
        )
    ) {
        ++end;
    }

    if (begin == end) {
        throw std::runtime_error(
            "FASTA record has an empty name."
        );
    }

    return description.substr(
        begin,
        end - begin
    );
}


}  // namespace


std::size_t
stream_selected_fasta_records(
    const std::string& path,
    const FastaRecordSelector& selector,
    const FastaRecordConsumer& consumer
) {
    if (!selector) {
        throw std::invalid_argument(
            "FASTA selector cannot be empty."
        );
    }

    if (!consumer) {
        throw std::invalid_argument(
            "FASTA consumer cannot be empty."
        );
    }


    std::ifstream input(
        path
    );

    if (!input) {
        throw std::runtime_error(
            "Cannot open FASTA: "
            +
            path
        );
    }


    FastaRecord current;

    bool active =
        false;

    bool selected =
        false;

    bool saw_header =
        false;

    std::size_t selected_records =
        0;

    std::string line;


    auto finish_current =
        [&]() {
            if (
                active
                &&
                selected
            ) {
                if (
                    current.sequence.empty()
                ) {
                    throw std::runtime_error(
                        "FASTA record contains no sequence: "
                        +
                        current.name
                    );
                }

                consumer(
                    std::move(
                        current
                    )
                );

                ++selected_records;
            }

            current =
                FastaRecord{};

            active =
                false;

            selected =
                false;
        };


    while (
        std::getline(
            input,
            line
        )
    ) {
        if (
            !line.empty()
            &&
            line.back() == '\r'
        ) {
            line.pop_back();
        }


        if (line.empty()) {
            continue;
        }


        if (
            line.front() == '>'
        ) {
            finish_current();

            current.description =
                line.substr(
                    1
                );

            current.name =
                header_name(
                    current.description
                );

            selected =
                selector(
                    current.name,
                    current.description
                );

            active =
                true;

            saw_header =
                true;

            continue;
        }


        if (!active) {
            throw std::runtime_error(
                "FASTA sequence encountered "
                "before first header."
            );
        }


        /*
         * Critical memory property:
         *
         * sequence from rejected records is not
         * accumulated or normalized.
         */
        if (!selected) {
            continue;
        }


        for (const char raw : line) {

            if (
                std::isspace(
                    static_cast<unsigned char>(
                        raw
                    )
                )
            ) {
                continue;
            }

            current.sequence.push_back(
                normalize_base(
                    raw
                )
            );
        }
    }


    finish_current();


    if (!saw_header) {
        throw std::runtime_error(
            "FASTA contains no records."
        );
    }


    return selected_records;
}


std::vector<FastaRecord>
load_fasta_records(
    const std::string& path
) {
    std::vector<FastaRecord>
        records;


    static_cast<void>(
        stream_selected_fasta_records(
            path,

            [](
                std::string_view,
                std::string_view
            ) {
                return true;
            },

            [&records](
                FastaRecord&& record
            ) {
                records.push_back(
                    std::move(
                        record
                    )
                );
            }
        )
    );


    return records;
}


}  // namespace primerpair
