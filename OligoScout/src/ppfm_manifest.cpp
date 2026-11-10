#include "primerpair/ppfm_manifest.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace primerpair {

namespace {

std::vector<std::string>
split_tsv(
    const std::string& line
) {
    std::vector<std::string> fields;

    std::size_t begin = 0;

    while (true) {
        const std::size_t tab =
            line.find(
                '\t',
                begin
            );

        if (
            tab ==
            std::string::npos
        ) {
            fields.push_back(
                line.substr(begin)
            );

            break;
        }

        fields.push_back(
            line.substr(
                begin,
                tab - begin
            )
        );

        begin =
            tab + 1;
    }

    return fields;
}


bool valid_sha256(
    const std::string& value
) {
    if (value.size() != 64) {
        return false;
    }

    for (const char c : value) {
        const bool digit =
            c >= '0' &&
            c <= '9';

        const bool lower_hex =
            c >= 'a' &&
            c <= 'f';

        const bool upper_hex =
            c >= 'A' &&
            c <= 'F';

        if (
            !digit &&
            !lower_hex &&
            !upper_hex
        ) {
            return false;
        }
    }

    return true;
}

}  // namespace


PpfmManifest PpfmManifest::load(
    const std::filesystem::path& manifest_path,
    const std::filesystem::path& index_directory
) {
    std::ifstream input(
        manifest_path
    );

    if (!input) {
        throw std::runtime_error(
            "Unable to open PPFM manifest: " +
            manifest_path.string()
        );
    }

    std::string header;

    if (
        !std::getline(
            input,
            header
        )
    ) {
        throw std::runtime_error(
            "PPFM manifest is empty."
        );
    }

    if (
        header !=
        "chromosome\tfile_bytes\tsha256\tstatus"
    ) {
        throw std::runtime_error(
            "Unexpected PPFM manifest header."
        );
    }

    PpfmManifest manifest;

    std::unordered_set<std::string>
        seen_chromosomes;

    std::string line;

    std::size_t line_number = 1;

    while (
        std::getline(
            input,
            line
        )
    ) {
        ++line_number;

        if (line.empty()) {
            continue;
        }

        const auto fields =
            split_tsv(
                line
            );

        if (fields.size() != 4) {
            throw std::runtime_error(
                "Invalid PPFM manifest row at line " +
                std::to_string(
                    line_number
                )
            );
        }

        PpfmManifestEntry entry;

        entry.chromosome =
            fields[0];

        if (entry.chromosome.empty()) {
            throw std::runtime_error(
                "Empty chromosome in PPFM manifest."
            );
        }

        if (
            !seen_chromosomes
                .insert(
                    entry.chromosome
                )
                .second
        ) {
            throw std::runtime_error(
                "Duplicate chromosome in PPFM manifest: " +
                entry.chromosome
            );
        }

        try {
            std::size_t consumed = 0;

            entry.file_bytes =
                std::stoull(
                    fields[1],
                    &consumed
                );

            if (
                consumed !=
                fields[1].size()
            ) {
                throw std::invalid_argument(
                    "trailing characters"
                );
            }

        } catch (...) {
            throw std::runtime_error(
                "Invalid file_bytes for chromosome " +
                entry.chromosome
            );
        }

        if (entry.file_bytes == 0) {
            throw std::runtime_error(
                "Zero file_bytes for chromosome " +
                entry.chromosome
            );
        }

        entry.sha256 =
            fields[2];

        if (
            !valid_sha256(
                entry.sha256
            )
        ) {
            throw std::runtime_error(
                "Invalid SHA256 for chromosome " +
                entry.chromosome
            );
        }

        entry.status =
            fields[3];

        if (
            entry.status != "BUILT" &&
            entry.status != "EXISTING"
        ) {
            throw std::runtime_error(
                "Unexpected manifest status for chromosome " +
                entry.chromosome
            );
        }

        entry.path =
            index_directory /
            (
                entry.chromosome +
                ".sa8.ppfm"
            );

        if (
            !std::filesystem::is_regular_file(
                entry.path
            )
        ) {
            throw std::runtime_error(
                "PPFM shard file missing: " +
                entry.path.string()
            );
        }

        const std::uint64_t actual_bytes =
            std::filesystem::file_size(
                entry.path
            );

        if (
            actual_bytes !=
            entry.file_bytes
        ) {
            throw std::runtime_error(
                "PPFM shard size mismatch for " +
                entry.chromosome
            );
        }

        manifest.entries_.push_back(
            std::move(
                entry
            )
        );
    }

    if (manifest.entries_.empty()) {
        throw std::runtime_error(
            "PPFM manifest contains no shards."
        );
    }

    return manifest;
}


const PpfmManifestEntry*
PpfmManifest::find(
    const std::string_view chromosome
) const noexcept {
    for (const auto& entry : entries_) {
        if (
            entry.chromosome ==
            chromosome
        ) {
            return &entry;
        }
    }

    return nullptr;
}

}  // namespace primerpair
