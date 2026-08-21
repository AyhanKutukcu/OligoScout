#include "primerpair/ppfm_manifest.hpp"

#include <iostream>
#include <stdexcept>

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

}  // namespace


int main() {
    using namespace primerpair;

    try {
        const auto manifest =
            PpfmManifest::load(
                "results/benchmarks/grch38_ppfm_batch/manifest.sa8.tsv",
                "data/indexes/ppfm"
            );

        expect(
            manifest.size() == 24,
            "Manifest contains 24 canonical shards"
        );

        expect(
            manifest.find("chr1") != nullptr,
            "chr1 present"
        );

        expect(
            manifest.find("chr21") != nullptr,
            "chr21 present"
        );

        expect(
            manifest.find("chr22") != nullptr,
            "chr22 present"
        );

        expect(
            manifest.find("chrX") != nullptr,
            "chrX present"
        );

        expect(
            manifest.find("chrY") != nullptr,
            "chrY present"
        );

        expect(
            manifest.find("chrM") == nullptr,
            "chrM correctly absent from canonical 24"
        );

        for (
            const auto& entry :
            manifest.entries()
        ) {
            expect(
                entry.file_bytes > 0,
                "Shard has non-zero byte size"
            );

            expect(
                entry.sha256.size() == 64,
                "Shard has SHA256 metadata"
            );

            expect(
                !entry.path.empty(),
                "Shard path resolved"
            );
        }

        std::cout
            << "manifest_entries\t"
            << manifest.size()
            << '\n';

        std::cout
            << "ALL_CHECKS\tYES\n";

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
