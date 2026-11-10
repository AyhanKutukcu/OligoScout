#include "primerpair/fasta_reader.hpp"
#include "primerpair/grch38_reference.hpp"
#include "primerpair/ppfm_io.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

double elapsed_seconds(
    const std::chrono::steady_clock::time_point begin,
    const std::chrono::steady_clock::time_point end
) {
    return
        std::chrono::duration<double>(
            end - begin
        ).count();
}

}  // namespace


int main(
    const int argc,
    char** argv
) {
    using namespace primerpair;

    try {
        if (
            argc != 4
            &&
            argc != 5
        ) {
            std::cerr
                << "Usage:\n"
                << "  build_grch38_ppfm "
                << "<GRCh38.fa> "
                << "<chromosome> "
                << "<output.ppfm> "
                << "[sa_rate=8]\n";

            return 2;
        }

        const std::string fasta_path =
            argv[1];

        const std::string chromosome =
            argv[2];

        const std::filesystem::path
            output_path =
                argv[3];

        const std::size_t sa_rate =
            argc == 5
                ? static_cast<std::size_t>(
                      std::stoull(
                          argv[4]
                      )
                  )
                : std::size_t{8};

        if (sa_rate == 0) {
            throw std::invalid_argument(
                "SA sample rate must be greater than zero."
            );
        }

        if (
            output_path.has_parent_path()
        ) {
            std::filesystem::
                create_directories(
                    output_path.parent_path()
                );
        }

        std::optional<FastaRecord>
            selected_record;

        const auto read_begin =
            std::chrono::
                steady_clock::now();

        const std::size_t matched =
            stream_selected_fasta_records(
                fasta_path,

                [&chromosome](
                    const std::string_view name,
                    const std::string_view description
                ) {
                    const auto alias =
                        grch38_primary_chromosome_alias(
                            name,
                            description
                        );

                    return
                        alias.has_value()
                        &&
                        *alias ==
                            chromosome;
                },

                [&selected_record](
                    FastaRecord&& record
                ) {
                    if (selected_record) {
                        throw std::runtime_error(
                            "Requested chromosome "
                            "matched more than once."
                        );
                    }

                    selected_record =
                        std::move(
                            record
                        );
                }
            );

        const auto read_end =
            std::chrono::
                steady_clock::now();

        if (
            matched != 1
            ||
            !selected_record.has_value()
        ) {
            throw std::runtime_error(
                "Requested GRCh38 primary chromosome "
                "was not found exactly once."
            );
        }

        std::string sequence =
            std::move(
                selected_record->sequence
            );

        const std::uint64_t
            sequence_length =
                static_cast<std::uint64_t>(
                    sequence.size()
                );

        std::cout
            << "chromosome\t"
            << chromosome
            << '\n';

        std::cout
            << "sequence_length\t"
            << sequence_length
            << '\n';

        std::cout
            << "sa_rate\t"
            << sa_rate
            << '\n';

        std::cout
            << "fasta_read_seconds\t"
            << elapsed_seconds(
                read_begin,
                read_end
            )
            << '\n';

        /*
         * PackedReference must be created before
         * sequence is moved to BidirectionalFMIndex.
         */
        const auto build_begin =
            std::chrono::
                steady_clock::now();

        PackedReference reference{
            sequence
        };

        BidirectionalFMIndex index{
            std::move(
                sequence
            ),
            sa_rate
        };

        const auto build_end =
            std::chrono::
                steady_clock::now();

        std::cout
            << "build_seconds\t"
            << elapsed_seconds(
                build_begin,
                build_end
            )
            << '\n';

        const auto save_begin =
            std::chrono::
                steady_clock::now();

        PpfmIO::save_shard(
            output_path,
            chromosome,
            reference,
            index
        );

        const auto save_end =
            std::chrono::
                steady_clock::now();

        std::cout
            << "save_seconds\t"
            << elapsed_seconds(
                save_begin,
                save_end
            )
            << '\n';

        std::cout
            << "ppfm_bytes\t"
            << std::filesystem::
                file_size(
                    output_path
                )
            << '\n';

        std::cout
            << "forward_sampled_sa_bytes\t"
            << index
                .forward_index()
                .sampled_sa_memory_bytes()
            << '\n';

        std::cout
            << "reverse_sampled_sa_bytes\t"
            << index
                .reverse_index()
                .sampled_sa_memory_bytes()
            << '\n';

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
