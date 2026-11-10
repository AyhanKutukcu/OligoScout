#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "primerpair/fasta_reader.hpp"

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


}  // namespace


int main() {
    try {
        const auto path =
            std::filesystem::
                temp_directory_path()
            /
            "primerpair_test_multifasta.fa";


        {
            std::ofstream output(
                path
            );

            output
                << ">chrA synthetic chromosome A\n"
                << "acgtACGTnnrysw\n"
                << ">chrB second chromosome\n"
                << "TTTTCCCCAAAAGGGG\n";
        }


        const auto records =
            primerpair::
                load_fasta_records(
                    path.string()
                );


        expect(
            records.size() == 2,
            "Two FASTA records loaded"
        );

        expect(
            records.at(0).name ==
                "chrA",
            "First FASTA record name"
        );

        expect(
            records.at(1).name ==
                "chrB",
            "Second FASTA record name"
        );

        expect(
            records.at(0).description ==
                "chrA synthetic chromosome A",
            "Complete FASTA description retained"
        );

        expect(
            records.at(0).sequence ==
                "ACGTACGTNNNNNN",
            "Lowercase and IUPAC normalization"
        );

        expect(
            records.at(1).sequence ==
                "TTTTCCCCAAAAGGGG",
            "Second sequence preserved"
        );


        std::size_t streamed_records = 0;

        std::string streamed_sequence;

        const std::size_t selected_count =
            primerpair::
                stream_selected_fasta_records(
                    path.string(),

                    [](
                        std::string_view name,
                        std::string_view
                    ) {
                        return name == "chrB";
                    },

                    [&](
                        primerpair::FastaRecord&& record
                    ) {
                        ++streamed_records;

                        streamed_sequence =
                            std::move(
                                record.sequence
                            );
                    }
                );


        expect(
            selected_count == 1,
            "Selective FASTA streaming selects one record"
        );

        expect(
            streamed_records == 1,
            "Selective FASTA consumer called once"
        );

        expect(
            streamed_sequence ==
                "TTTTCCCCAAAAGGGG",
            "Selective FASTA sequence correct"
        );


        std::filesystem::remove(
            path
        );


        std::cout
            << "Multi-FASTA reader tests passed.\n";

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
