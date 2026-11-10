#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "primerpair/fm_index.hpp"

namespace {

void print_usage(const char* executable) {
    std::cout
        << "OligoScout\n\n"
        << "Usage:\n"
        << "  " << executable << " demo\n"
        << "  " << executable
        << " build --reference <reference.fa>"
        << " --output <index.ppfm>\n"
        << "  " << executable
        << " search --index <index.ppfm>"
        << " --pattern <DNA> [--count-only]\n"
        << "  " << executable
        << " search-genome --index-dir <directory>"
        << " --pattern <DNA> [--count-only]\n";
}

int run_demo() {
    const std::string reference =
        "ACGTACGTACGTGATTACAGATTACACCCCC"
        "GGGGGTTTTTAAAAAACGTACGTACGT";

    const std::string pattern =
        "ACGTACGTACGT";

    const primerpair::FMIndex index(reference);

    const auto interval =
        index.backward_search(pattern);

    const std::vector<std::uint64_t> positions =
        index.locate(interval);

    std::cout
        << "OligoScout exact FM-index test\n";

#ifdef PRIMERPAIR_HAVE_AVX2_BACKEND
    std::cout << "AVX2 backend: compiled\n";
#else
    std::cout << "AVX2 backend: not compiled\n";
#endif

    /*
     * GATTACA, BWT doğrulaması için ayrı küçük indeks.
     */
    const primerpair::FMIndex gattaca_index("GATTACA");

    std::cout
        << "reference\tGATTACA\n"
        << "bwt\t"
        << gattaca_index.bwt_string()
        << '\n';

    std::cout
        << "pattern\t"
        << pattern
        << '\n';

    std::cout
        << "interval_begin\t"
        << interval.begin
        << '\n';

    std::cout
        << "interval_end\t"
        << interval.end
        << '\n';

    std::cout
        << "match_count\t"
        << interval.size()
        << '\n';

    std::cout
        << "start_1based\tend_1based\n";

    for (const std::uint64_t position : positions) {
        const std::uint64_t start_1based =
            position + 1;

        const std::uint64_t end_1based =
            position +
            static_cast<std::uint64_t>(
                pattern.size()
            );

        std::cout
            << start_1based
            << '\t'
            << end_1based
            << '\n';
    }

    const auto absent_interval =
        index.backward_search("TTTTCCCCAAAA");

    std::cout
        << "absent_match_count\t"
        << absent_interval.size()
        << '\n';

    return 0;
}

int not_implemented(
    std::string_view command
) {
    std::cerr
        << "Error: the '"
        << command
        << "' command is not implemented yet.\n"
        << "No FASTA file was read and no index was created.\n"
        << "Current implementation supports only the demo command.\n";

    return 2;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc == 1) {
            print_usage(argv[0]);
            return 0;
        }

        const std::string_view command =
            argv[1];

        if (command == "demo") {
            return run_demo();
        }

        if (command == "build") {
            return not_implemented(command);
        }

        if (command == "search") {
            return not_implemented(command);
        }

        if (command == "search-genome") {
            return not_implemented(command);
        }

        if (command == "--help" ||
            command == "-h" ||
            command == "help") {

            print_usage(argv[0]);
            return 0;
        }

        std::cerr
            << "Error: unknown command: "
            << command
            << "\n\n";

        print_usage(argv[0]);
        return 2;

    } catch (const std::exception& exception) {
        std::cerr
            << "Fatal error: "
            << exception.what()
            << '\n';

        return 1;
    }
}
