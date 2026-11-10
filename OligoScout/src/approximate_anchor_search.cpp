#include "primerpair/approximate_anchor_search.hpp"

#include <array>
#include <cctype>
#include <stdexcept>
#include <string>
#include <utility>

namespace primerpair {

namespace {

constexpr std::array<char, 4>
    kDnaAlphabet{
        'A',
        'C',
        'G',
        'T'
    };

std::string normalize_primer(
    const std::string_view primer
) {
    if (primer.empty()) {
        throw std::invalid_argument(
            "Primer cannot be empty."
        );
    }

    std::string normalized;
    normalized.reserve(
        primer.size()
    );

    for (const char base : primer) {
        const char upper =
            static_cast<char>(
                std::toupper(
                    static_cast<unsigned char>(
                        base
                    )
                )
            );

        switch (upper) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                normalized.push_back(
                    upper
                );
                break;

            default:
                throw std::invalid_argument(
                    "Approximate primer search currently "
                    "accepts only A/C/G/T."
                );
        }
    }

    return normalized;
}

}  // namespace

ApproximateAnchorSearcher::
ApproximateAnchorSearcher(
    const BidirectionalFMIndex& index
) noexcept
    : index_(index) {
}

ApproximateAnchorSearchResult
ApproximateAnchorSearcher::
search_5prime_mismatches(
    const std::string_view primer,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    if (anchor_length == 0) {
        throw std::invalid_argument(
            "Anchor length must be greater than zero."
        );
    }

    if (max_mismatches > 3) {
        throw std::invalid_argument(
            "Current PrimerPair MVP supports "
            "at most 3 mismatches per primer."
        );
    }

    const std::string normalized =
        normalize_primer(
            primer
        );

    if (anchor_length > normalized.size()) {
        throw std::invalid_argument(
            "Anchor length cannot exceed primer length."
        );
    }

    const std::size_t anchor_begin =
        normalized.size() -
        anchor_length;

    const std::string_view anchor(
        normalized.data() +
            anchor_begin,
        anchor_length
    );

    const BidirectionalInterval
        anchor_state =
            index_.search(
                anchor
            );

    ApproximateAnchorSearchResult result;

    result.primer_length =
        normalized.size();

    result.anchor_length =
        anchor_length;

    result.max_mismatches =
        max_mismatches;

    /*
     * Exact 3' anchor genomda yoksa bütün
     * approximate arama burada biter.
     */
    if (anchor_state.empty()) {
        return result;
    }

    std::vector<ApproximateAnchorHit>
        current;

    current.push_back(
        ApproximateAnchorHit{
            anchor_state,
            0,
            {}
        }
    );

    /*
     * Anchor'ın hemen solundan primerin
     * 5' ucuna doğru ilerliyoruz.
     */
    for (std::size_t i = anchor_begin;
         i > 0;
         --i) {

        const std::size_t primer_position =
            i - 1;

        const char expected_base =
            normalized.at(
                primer_position
            );

        std::vector<ApproximateAnchorHit>
            next;

        /*
         * En kötü durumda her branch dört dala
         * ayrılabilir.
         */
        next.reserve(
            current.size() * 4
        );

        for (const auto& branch : current) {

            for (const char genomic_base :
                 kDnaAlphabet) {

                const bool is_mismatch =
                    genomic_base !=
                    expected_base;

                const std::size_t
                    new_mismatch_count =
                        branch.mismatches +
                        (
                            is_mismatch
                                ? 1
                                : 0
                        );

                if (
                    new_mismatch_count >
                    max_mismatches
                ) {
                    continue;
                }

                const BidirectionalInterval
                    extended =
                        index_.extend_left(
                            branch.state,
                            genomic_base
                        );

                /*
                 * Genomda bu prefix yoksa branch
                 * hemen öldürülür.
                 */
                if (extended.empty()) {
                    continue;
                }

                ApproximateAnchorHit
                    new_hit;

                new_hit.state =
                    extended;

                new_hit.mismatches =
                    new_mismatch_count;

                new_hit.mismatch_positions =
                    branch.mismatch_positions;

                if (is_mismatch) {
                    new_hit
                        .mismatch_positions
                        .push_back(
                            primer_position
                        );
                }

                next.push_back(
                    std::move(
                        new_hit
                    )
                );
            }
        }

        current =
            std::move(
                next
            );

        /*
         * Hiç canlı branch kalmadıysa erken bitir.
         */
        if (current.empty()) {
            break;
        }
    }

    result.hits =
        std::move(
            current
        );

    return result;
}

std::vector<std::uint64_t>
ApproximateAnchorSearcher::locate(
    const ApproximateAnchorHit& hit
) const {
    if (hit.empty()) {
        return {};
    }

    return index_.locate(
        hit.state
    );
}

}  // namespace primerpair
