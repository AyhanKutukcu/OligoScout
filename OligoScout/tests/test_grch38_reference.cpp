#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include "primerpair/grch38_reference.hpp"

namespace {


void expect(
    const bool condition,
    const std::string& message
) {
    if (!condition) {
        throw std::runtime_error(
            "FAILED: "
            +
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
        using primerpair::
            grch38_primary_chromosome_alias;


        const auto chr1 =
            grch38_primary_chromosome_alias(
                "NC_000001.11",
                "NC_000001.11 Homo sapiens chromosome 1, "
                "GRCh38.p14 Primary Assembly"
            );

        expect(
            chr1 &&
            *chr1 == "chr1",
            "NCBI chromosome 1 resolves to chr1"
        );


        const auto chr22 =
            grch38_primary_chromosome_alias(
                "NC_000022.11",
                "NC_000022.11 Homo sapiens chromosome 22, "
                "GRCh38.p14 Primary Assembly"
            );

        expect(
            chr22 &&
            *chr22 == "chr22",
            "NCBI chromosome 22 resolves to chr22"
        );


        const auto chr_x =
            grch38_primary_chromosome_alias(
                "NC_000023.11",
                "NC_000023.11 Homo sapiens chromosome X, "
                "GRCh38.p14 Primary Assembly"
            );

        expect(
            chr_x &&
            *chr_x == "chrX",
            "NCBI chromosome X resolves to chrX"
        );


        const auto chr_y =
            grch38_primary_chromosome_alias(
                "NC_000024.10",
                "NC_000024.10 Homo sapiens chromosome Y, "
                "GRCh38.p14 Primary Assembly"
            );

        expect(
            chr_y &&
            *chr_y == "chrY",
            "NCBI chromosome Y resolves to chrY"
        );


        const auto unlocalized =
            grch38_primary_chromosome_alias(
                "NT_187361.1",
                "NT_187361.1 Homo sapiens chromosome 1 "
                "unlocalized genomic scaffold, "
                "GRCh38.p14 Primary Assembly"
            );

        expect(
            !unlocalized,
            "Unlocalized chromosome record rejected"
        );


        const auto wrong_assembly =
            grch38_primary_chromosome_alias(
                "NC_000001.11",
                "NC_000001.11 Homo sapiens chromosome 1, "
                "GRCh38.p13 Primary Assembly"
            );

        expect(
            !wrong_assembly,
            "Wrong assembly context rejected"
        );


        const auto wrong_assembly_x =
            grch38_primary_chromosome_alias(
                "NC_000023.11",
                "NC_000023.11 Homo sapiens chromosome X, "
                "GRCh38.p13 Primary Assembly"
            );

        expect(
            !wrong_assembly_x,
            "Wrong assembly chrX rejected"
        );


        const auto wrong_assembly_y =
            grch38_primary_chromosome_alias(
                "NC_000024.10",
                "NC_000024.10 Homo sapiens chromosome Y, "
                "GRCh38.p13 Primary Assembly"
            );

        expect(
            !wrong_assembly_y,
            "Wrong assembly chrY rejected"
        );


        std::cout
            << "GRCh38 chromosome resolver tests passed.\n";

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
