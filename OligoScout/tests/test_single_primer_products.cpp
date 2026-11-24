// Single-primer products (P-> <-P): geometry, length window, identity and
// equivalence with a brute-force pairing over the same hits.
#include "primerpair/primer_pair_assembler.hpp"
#include "primerpair/primer_pair_search.hpp"
#include "primerpair/strand_aware_primer_search.hpp"

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(const bool ok, const std::string& what) {
    if (!ok) {
        throw std::runtime_error(what);
    }
}

}  // namespace

int main() {
    try {
        using namespace primerpair;

        const std::string primer(20, 'A');
        const auto F = PrimerOrientation::Forward;
        const auto R = PrimerOrientation::Reverse;

        // forward at 100 and 5000; reverse at 100 (overlap), 130 (too short:
        // 130+20-100 = 50 -> exactly min, kept), 200, 2900 (2900+20-100 = 2820),
        // 3200 (too long from 100), 5100.
        const std::vector<OrientedPrimerSearchHit> hits = {
            {100, 0, F, 0}, {5000, 1, F, 1},
            {100, 0, R, 0}, {130, 0, R, 0}, {200, 2, R, 6},
            {2900, 0, R, 0}, {3200, 0, R, 0}, {5100, 0, R, 0},
        };

        const auto products = assemble_single_primer_hits(
            primer, hits, PrimerIdentity::Primer2, 50, 3000);

        // brute force with the same geometry rules
        std::size_t expected = 0;
        for (const auto& f : hits) {
            if (f.orientation != F) continue;
            for (const auto& r : hits) {
                if (r.orientation != R) continue;
                if (r.position < f.position + primer.size()) continue;
                const auto len = r.position + primer.size() - f.position;
                if (len >= 50 && len <= 3000) ++expected;
            }
        }

        require(products.size() == expected, "count differs from brute force");
        require(products.size() == 4, "expected 4 single-primer products");

        for (const auto& p : products) {
            require(p.left_primer == PrimerIdentity::Primer2 &&
                    p.right_primer == PrimerIdentity::Primer2,
                    "identity must be the single primer on both sides");
            require(p.amplicon_length >= 50 && p.amplicon_length <= 3000,
                    "length window violated");
            require(p.right_position >= p.left_position + primer.size(),
                    "overlapping geometry accepted");
        }

        // the 200 bp product carries the reverse site's mismatches
        bool found = false;
        for (const auto& p : products) {
            if (p.left_position == 100 && p.right_position == 200) {
                found = p.right_mismatches == 2 && p.right_mismatch_mask == 6 &&
                        p.amplicon_length == 120;
            }
        }
        require(found, "100->200 product missing or wrong");

        // no forward or no reverse hits -> nothing
        const std::vector<OrientedPrimerSearchHit> only_fwd = {{100, 0, F, 0}};
        require(assemble_single_primer_hits(primer, only_fwd, PrimerIdentity::Primer1, 50, 3000).empty(),
                "forward-only hits must give no product");

        std::cout << "single_primer_products: PASS (" << products.size() << " products)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
