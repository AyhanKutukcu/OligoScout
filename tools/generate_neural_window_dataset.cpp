#include "primerpair/ip_bwt_index.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

std::string make_reference(
    const std::size_t length
) {
    if (length < 1000) {
        throw std::invalid_argument(
            "Reference must contain at least 1000 bases."
        );
    }


    constexpr char alphabet[] = {
        'A',
        'C',
        'G',
        'T'
    };


    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;


    std::string reference;

    reference.reserve(
        length
    );


    /*
     * Deterministic pseudo-random background.
     */
    for (
        std::size_t i = 0;
        i < length;
        ++i
    ) {
        state =
            state *
            6364136223846793005ULL
            +
            1442695040888963407ULL;


        const std::size_t index =
            static_cast<std::size_t>(
                (
                    state >>
                    32
                )
                &
                3ULL
            );


        reference.push_back(
            alphabet[index]
        );
    }


    /*
     * Deliberately inject repetitive and
     * low-complexity regions.
     *
     * A neural search-window model must learn
     * both unique and repetitive sequence regimes.
     */
    const std::string motif =
        "ACGTACGTACGTACGT"
        "AAAAAAAACCCCCCCC"
        "GGGGGGGGTTTTTTTT"
        "ACACACACGTGTGTGT";


    for (
        std::size_t block_start = 8192;
        block_start + 512 < reference.size();
        block_start += 16384
    ) {
        for (
            std::size_t offset = 0;
            offset < 512;
            ++offset
        ) {
            reference.at(
                block_start + offset
            ) =
                motif.at(
                    offset %
                    motif.size()
                );
        }
    }


    return reference;
}


char mutate_base(
    const char base
) {
    switch (base) {

        case 'A':
            return 'C';

        case 'C':
            return 'G';

        case 'G':
            return 'T';

        case 'T':
            return 'A';

        default:
            throw std::runtime_error(
                "Unexpected nucleotide."
            );
    }
}


double gc_fraction(
    const std::string_view sequence
) {
    std::size_t gc = 0;


    for (
        const char base :
        sequence
    ) {
        if (
            base == 'G'
            ||
            base == 'C'
        ) {
            ++gc;
        }
    }


    return
        static_cast<double>(
            gc
        )
        /
        static_cast<double>(
            sequence.size()
        );
}


double shannon_entropy(
    const std::string_view sequence
) {
    std::array<std::size_t, 4>
        counts{};


    for (
        const char base :
        sequence
    ) {
        switch (base) {

            case 'A':
                ++counts[0];
                break;

            case 'C':
                ++counts[1];
                break;

            case 'G':
                ++counts[2];
                break;

            case 'T':
                ++counts[3];
                break;

            default:
                throw std::runtime_error(
                    "Unexpected nucleotide."
                );
        }
    }


    double entropy = 0.0;


    for (
        const std::size_t count :
        counts
    ) {
        if (
            count == 0
        ) {
            continue;
        }


        const double probability =
            static_cast<double>(
                count
            )
            /
            static_cast<double>(
                sequence.size()
            );


        entropy -=
            probability *
            std::log2(
                probability
            );
    }


    return entropy;
}


std::size_t maximum_homopolymer_run(
    const std::string_view sequence
) {
    if (
        sequence.empty()
    ) {
        return 0;
    }


    std::size_t maximum = 1;
    std::size_t current = 1;


    for (
        std::size_t i = 1;
        i < sequence.size();
        ++i
    ) {
        if (
            sequence.at(i)
            ==
            sequence.at(
                i - 1
            )
        ) {
            ++current;

            if (
                current >
                maximum
            ) {
                maximum =
                    current;
            }

        } else {

            current = 1;
        }
    }


    return maximum;
}


std::size_t distinct_base_count(
    const std::string_view sequence
) {
    bool a = false;
    bool c = false;
    bool g = false;
    bool t = false;


    for (
        const char base :
        sequence
    ) {
        switch (base) {

            case 'A':
                a = true;
                break;

            case 'C':
                c = true;
                break;

            case 'G':
                g = true;
                break;

            case 'T':
                t = true;
                break;

            default:
                break;
        }
    }


    return
        static_cast<std::size_t>(a)
        +
        static_cast<std::size_t>(c)
        +
        static_cast<std::size_t>(g)
        +
        static_cast<std::size_t>(t);
}


std::uint64_t splitmix64(
    std::uint64_t value
) {
    value +=
        0x9E3779B97F4A7C15ULL;

    value =
        (
            value ^
            (
                value >>
                30U
            )
        )
        *
        0xBF58476D1CE4E5B9ULL;

    value =
        (
            value ^
            (
                value >>
                27U
            )
        )
        *
        0x94D049BB133111EBULL;

    return
        value ^
        (
            value >>
            31U
        );
}


const char* dataset_split(
    const std::size_t source_start
) {
    /*
     * Split by genomic block rather than by
     * individual query.
     *
     * Nearby overlapping k-mers therefore stay
     * together and cannot trivially leak from
     * training into validation/test.
     */
    const std::uint64_t block =
        static_cast<std::uint64_t>(
            source_start /
            4096
        );


    const std::uint64_t bucket =
        splitmix64(
            block
        )
        %
        10ULL;


    if (
        bucket <
        8ULL
    ) {
        return "train";
    }


    if (
        bucket ==
        8ULL
    ) {
        return "validation";
    }


    return "test";
}

}  // namespace


int main(
    const int argc,
    char** argv
) {
    using namespace primerpair;

    try {
        if (
            argc !=
            5
        ) {
            std::cerr
                << "Usage:\n"
                << "  "
                << argv[0]
                << " <reference_length>"
                << " <sample_count>"
                << " <k>"
                << " <output.tsv>\n";

            return 2;
        }


        const std::size_t reference_length =
            static_cast<std::size_t>(
                std::stoull(
                    argv[1]
                )
            );


        const std::size_t sample_count =
            static_cast<std::size_t>(
                std::stoull(
                    argv[2]
                )
            );


        const std::size_t k =
            static_cast<std::size_t>(
                std::stoull(
                    argv[3]
                )
            );


        const std::string output_path =
            argv[4];


        if (
            k == 0
            ||
            k >
            IPBWTIndex::max_chunk_length
        ) {
            throw std::invalid_argument(
                "Training k must be between 1 and 21."
            );
        }


        if (
            reference_length <=
            k
        ) {
            throw std::invalid_argument(
                "Reference is shorter than k."
            );
        }


        if (
            sample_count == 0
        ) {
            throw std::invalid_argument(
                "Sample count cannot be zero."
            );
        }


        std::cout
            << "reference_length\t"
            << reference_length
            << '\n';

        std::cout
            << "sample_count\t"
            << sample_count
            << '\n';

        std::cout
            << "k\t"
            << k
            << '\n';


        const std::string reference =
            make_reference(
                reference_length
            );


        std::cout
            << "building_ip_bwt\tYES\n";


        IPBWTIndex index(
            reference,
            k
        );


        std::cout
            << "ip_bwt_rows\t"
            << index.row_count()
            << '\n';

        std::cout
            << "ip_bwt_storage_bytes\t"
            << index.compact_storage_bytes()
            << '\n';


        std::ofstream output(
            output_path
        );


        if (
            !output
        ) {
            throw std::runtime_error(
                "Could not open output dataset."
            );
        }


        output
            << "sample_id"
            << '\t'
            << "split"
            << '\t'
            << "kmer"
            << '\t'
            << "source_start"
            << '\t'
            << "mutated"
            << '\t'
            << "present"
            << '\t'
            << "sa_lower"
            << '\t'
            << "sa_upper"
            << '\t'
            << "sa_width"
            << '\t'
            << "lower_norm"
            << '\t'
            << "upper_norm"
            << '\t'
            << "center_norm"
            << '\t'
            << "log2_width_plus1"
            << '\t'
            << "gc_fraction"
            << '\t'
            << "entropy"
            << '\t'
            << "max_homopolymer"
            << '\t'
            << "distinct_bases"
            << '\n';


        output
            << std::setprecision(
                10
            );


        std::uint64_t random_state =
            0xD1B54A32D192ED03ULL;


        std::size_t actual_present = 0;
        std::size_t actual_absent = 0;

        std::size_t train_count = 0;
        std::size_t validation_count = 0;
        std::size_t test_count = 0;


        for (
            std::size_t sample = 0;
            sample < sample_count;
            ++sample
        ) {
            random_state =
                random_state *
                2862933555777941757ULL
                +
                3037000493ULL;


            const std::size_t maximum_start =
                reference.size() -
                k;


            const std::size_t source_start =
                static_cast<std::size_t>(
                    random_state %
                    static_cast<std::uint64_t>(
                        maximum_start + 1
                    )
                );


            std::string query =
                reference.substr(
                    source_start,
                    k
                );


            const bool mutated =
                (
                    sample &
                    std::size_t{1}
                )
                != 0;


            if (
                mutated
            ) {
                random_state =
                    random_state *
                    2862933555777941757ULL
                    +
                    3037000493ULL;


                const std::size_t mutation_position =
                    static_cast<std::size_t>(
                        random_state %
                        static_cast<std::uint64_t>(
                            k
                        )
                    );


                query.at(
                    mutation_position
                ) =
                    mutate_base(
                        query.at(
                            mutation_position
                        )
                    );
            }


            const Interval interval =
                index.exact_search(
                    query
                );


            const std::uint64_t width =
                interval.size();


            const bool present =
                width != 0;


            if (
                present
            ) {
                ++actual_present;
            } else {
                ++actual_absent;
            }


            const double row_count =
                static_cast<double>(
                    index.row_count()
                );


            const double lower_norm =
                static_cast<double>(
                    interval.begin
                )
                /
                row_count;


            const double upper_norm =
                static_cast<double>(
                    interval.end
                )
                /
                row_count;


            const double center_norm =
                (
                    static_cast<double>(
                        interval.begin
                    )
                    +
                    static_cast<double>(
                        interval.end
                    )
                )
                /
                (
                    2.0 *
                    row_count
                );


            const double log_width =
                std::log2(
                    static_cast<double>(
                        width
                    )
                    +
                    1.0
                );


            const char* split =
                dataset_split(
                    source_start
                );


            if (
                std::string_view(split) ==
                "train"
            ) {
                ++train_count;

            } else if (
                std::string_view(split) ==
                "validation"
            ) {
                ++validation_count;

            } else {

                ++test_count;
            }


            output
                << sample
                << '\t'
                << split
                << '\t'
                << query
                << '\t'
                << source_start
                << '\t'
                << (
                    mutated
                    ? 1
                    : 0
                )
                << '\t'
                << (
                    present
                    ? 1
                    : 0
                )
                << '\t'
                << interval.begin
                << '\t'
                << interval.end
                << '\t'
                << width
                << '\t'
                << lower_norm
                << '\t'
                << upper_norm
                << '\t'
                << center_norm
                << '\t'
                << log_width
                << '\t'
                << gc_fraction(
                    query
                )
                << '\t'
                << shannon_entropy(
                    query
                )
                << '\t'
                << maximum_homopolymer_run(
                    query
                )
                << '\t'
                << distinct_base_count(
                    query
                )
                << '\n';
        }


        output.close();


        if (
            !output
        ) {
            throw std::runtime_error(
                "Dataset output failed during writing."
            );
        }


        std::cout
            << "actual_present\t"
            << actual_present
            << '\n';

        std::cout
            << "actual_absent\t"
            << actual_absent
            << '\n';

        std::cout
            << "train_samples\t"
            << train_count
            << '\n';

        std::cout
            << "validation_samples\t"
            << validation_count
            << '\n';

        std::cout
            << "test_samples\t"
            << test_count
            << '\n';

        std::cout
            << "output\t"
            << output_path
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
