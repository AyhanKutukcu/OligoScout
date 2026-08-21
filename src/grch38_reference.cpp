#include "primerpair/grch38_reference.hpp"

#include <cctype>
#include <optional>
#include <string>
#include <string_view>

namespace primerpair {


std::optional<std::string>
grch38_primary_chromosome_alias(
    const std::string_view record_name,
    const std::string_view description
) {
    /*
     * Canonical NCBI chromosome records use NC_.
     */
    if (
        !record_name.starts_with(
            "NC_"
        )
    ) {
        return std::nullopt;
    }


    static constexpr
        std::string_view marker =
            "Homo sapiens chromosome ";


    const std::size_t begin =
        description.find(
            marker
        );

    if (
        begin ==
        std::string_view::npos
    ) {
        return std::nullopt;
    }


    const std::size_t chromosome_begin =
        begin +
        marker.size();


    const std::size_t comma =
        description.find(
            ',',
            chromosome_begin
        );

    if (
        comma ==
        std::string_view::npos
    ) {
        return std::nullopt;
    }


    const std::string_view token =
        description.substr(
            chromosome_begin,
            comma -
            chromosome_begin
        );


    /*
     * Require the GRCh38.p14 primary-assembly
     * context before accepting ANY canonical
     * chromosome, including X and Y.
     */
    if (
        description.find(
            "GRCh38.p14 Primary Assembly"
        )
        ==
        std::string_view::npos
    ) {
        return std::nullopt;
    }


    /*
     * Reject unlocalized/alternate descriptions
     * naturally because token must be exactly one
     * canonical chromosome label.
     */
    if (
        token == "X"
        ||
        token == "Y"
    ) {
        return
            std::string(
                "chr"
            )
            +
            std::string(
                token
            );
    }


    if (
        token.empty()
        ||
        token.size() > 2
    ) {
        return std::nullopt;
    }


    unsigned int chromosome =
        0;


    for (const char c : token) {

        if (
            !std::isdigit(
                static_cast<unsigned char>(
                    c
                )
            )
        ) {
            return std::nullopt;
        }

        chromosome =
            chromosome * 10u
            +
            static_cast<unsigned int>(
                c - '0'
            );
    }


    if (
        chromosome < 1
        ||
        chromosome > 22
    ) {
        return std::nullopt;
    }


    return
        std::string(
            "chr"
        )
        +
        std::to_string(
            chromosome
        );
}


}  // namespace primerpair
