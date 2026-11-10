#include "primerpair/search_profile.hpp"

#include <stdexcept>

namespace primerpair {

const char* to_string(
    const SearchProfile profile
) noexcept {
    switch (profile) {

        case SearchProfile::Strict:
            return "STRICT";

        case SearchProfile::Sensitive:
            return "SENSITIVE";
    }

    return "UNKNOWN";
}


ProfiledPrimerSearchEngine::
ProfiledPrimerSearchEngine(
    const BidirectionalFMIndex& index,
    const PackedReference& reference
)
    : strict_searcher_(
          index,
          reference
      ),
      sensitive_searcher_(
          index,
          reference
      ) {
}


ProfiledPrimerSearchResult
ProfiledPrimerSearchEngine::search(
    const std::string_view primer,
    const SearchProfile profile,
    const std::size_t anchor_length,
    const std::size_t max_mismatches
) const {
    switch (profile) {

        case SearchProfile::Strict: {
            ProfiledPrimerSearchResult
                result;

            result.profile =
                SearchProfile::Strict;

            /*
             * Preserve the validated legacy path
             * exactly.
             */
            result.search_result =
                strict_searcher_.search(
                    primer,
                    anchor_length,
                    max_mismatches
                );

            result.sensitive_result_available =
                false;

            return result;
        }


        case SearchProfile::Sensitive: {
            ProfiledPrimerSearchResult
                result;

            result.profile =
                SearchProfile::Sensitive;

            /*
             * IMPORTANT:
             *
             * SENSITIVE intentionally ignores
             * anchor_length because it allows
             * mismatches anywhere in the full
             * primer, including the biological
             * 3-prime region.
             */
            result.sensitive_result =
                sensitive_searcher_.search(
                    primer,
                    max_mismatches
                );

            result.sensitive_result_available =
                true;

            return result;
        }
    }

    throw std::logic_error(
        "Unknown search profile."
    );
}

}  // namespace primerpair
