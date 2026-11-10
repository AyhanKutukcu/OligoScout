#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace primerpair {

struct PpfmManifestEntry {
    std::string chromosome;
    std::uint64_t file_bytes{0};
    std::string sha256;
    std::string status;
    std::filesystem::path path;
};


class PpfmManifest {
public:
    [[nodiscard]]
    static PpfmManifest load(
        const std::filesystem::path& manifest_path,
        const std::filesystem::path& index_directory
    );

    [[nodiscard]]
    std::size_t size() const noexcept {
        return entries_.size();
    }

    [[nodiscard]]
    bool empty() const noexcept {
        return entries_.empty();
    }

    [[nodiscard]]
    const std::vector<PpfmManifestEntry>&
    entries() const noexcept {
        return entries_;
    }

    [[nodiscard]]
    const PpfmManifestEntry& at(
        std::size_t index
    ) const {
        return entries_.at(index);
    }

    [[nodiscard]]
    const PpfmManifestEntry* find(
        std::string_view chromosome
    ) const noexcept;

private:
    std::vector<PpfmManifestEntry> entries_;
};

}  // namespace primerpair
