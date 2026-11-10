#include "primerpair/anchor_search.hpp"

#include <stdexcept>

namespace primerpair {

AnchorSearcher::AnchorSearcher(
    const BidirectionalFMIndex& index
) noexcept
    : index_(index) {
}

AnchorSearchResult AnchorSearcher::search_exact(
    const std::string_view primer,
    const std::size_t anchor_length
) const {
    if (primer.empty()) {
        throw std::invalid_argument(
            "Primer cannot be empty."
        );
    }

    if (anchor_length == 0) {
        throw std::invalid_argument(
            "Anchor length must be greater than zero."
        );
    }

    if (anchor_length > primer.size()) {
        throw std::invalid_argument(
            "Anchor length cannot exceed primer length."
        );
    }

    const std::size_t anchor_begin =
        primer.size() -
        anchor_length;

    const std::string_view anchor =
        primer.substr(
            anchor_begin,
            anchor_length
        );

    BidirectionalInterval state =
        index_.search(
            anchor
        );

    AnchorSearchResult result{
        state,
        primer.size(),
        anchor_length,
        0
    };

    /*
     * Anchor genomda hiç yoksa primerin geri kalanını
     * extend etmeye gerek yok.
     */
    if (state.empty()) {
        return result;
    }

    /*
     * 3' anchor'ın solundaki bazları 5' yönüne doğru
     * ters sırayla BiFM state'e ekleriz.
     *
     * primer:
     *
     * [0 ... anchor_begin-1][anchor]
     *
     * İlk extension:
     * primer[anchor_begin - 1]
     */
    for (std::size_t i = anchor_begin;
         i > 0;
         --i) {

        state =
            index_.extend_left(
                state,
                primer.at(i - 1)
            );

        ++result.extension_steps;

        if (state.empty()) {
            break;
        }
    }

    result.state =
        state;

    return result;
}

std::vector<std::uint64_t>
AnchorSearcher::locate(
    const AnchorSearchResult& result
) const {
    if (result.empty()) {
        return {};
    }

    return index_.locate(
        result.state
    );
}

}  // namespace primerpair
