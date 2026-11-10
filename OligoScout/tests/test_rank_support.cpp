#include <cstdint>
#include <iostream>
#include <string>

#include "primerpair/rank_support.hpp"

namespace {

bool check(
    bool condition,
    const std::string& name
) {
    if (condition) {
        std::cout << "[PASS] " << name << '\n';
        return true;
    }

    std::cerr << "[FAIL] " << name << '\n';
    return false;
}

}  // namespace

int main() {
    int failures = 0;

    /*
     * 300 bazlık BWT-benzeri test dizisi:
     *
     * İlk 128 baz:   A
     * Sonraki 128:   C
     * Son 44 baz:    G
     */
    const std::string sequence =
        std::string(128, 'A') +
        std::string(128, 'C') +
        std::string(44, 'G');

    primerpair::CheckpointRank rank;

    rank.build(sequence, 128);

    if (!check(
            rank.checkpoint_count() == 3,
            "Checkpoint count"
        )) {
        ++failures;
    }

    if (!check(
            rank.rank(sequence, 'A', 0) == 0,
            "Rank A at position 0"
        )) {
        ++failures;
    }

    if (!check(
            rank.rank(sequence, 'A', 128) == 128,
            "Rank A at first checkpoint"
        )) {
        ++failures;
    }

    if (!check(
            rank.rank(sequence, 'A', 129) == 128,
            "Rank A after first checkpoint"
        )) {
        ++failures;
    }

    if (!check(
            rank.rank(sequence, 'C', 128) == 0,
            "Rank C before C region"
        )) {
        ++failures;
    }

    if (!check(
            rank.rank(sequence, 'C', 256) == 128,
            "Rank C at second checkpoint"
        )) {
        ++failures;
    }

    if (!check(
            rank.rank(sequence, 'G', 256) == 0,
            "Rank G before final region"
        )) {
        ++failures;
    }

    if (!check(
            rank.rank(sequence, 'G', 300) == 44,
            "Rank G at sequence end"
        )) {
        ++failures;
    }

    if (!check(
            rank.memory_bytes() ==
                rank.checkpoint_count() *
                sizeof(std::array<std::uint32_t, 6>),
            "Checkpoint memory calculation"
        )) {
        ++failures;
    }

    if (failures != 0) {
        std::cerr
            << failures
            << " rank test(s) failed.\n";

        return 1;
    }

    std::cout
        << "All checkpoint-rank tests passed.\n";

    return 0;
}
