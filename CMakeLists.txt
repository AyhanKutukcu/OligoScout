cmake_minimum_required(VERSION 3.20)

project(
    OligoScout
    VERSION 0.1.0
    LANGUAGES CXX
)

option(
    PRIMERPAIR_ENABLE_AVX2
    "Build the optional AVX2 backend"
    OFF
)

option(
    PRIMERPAIR_ENABLE_SANITIZERS
    "Enable AddressSanitizer and UndefinedBehaviorSanitizer"
    OFF
)

# ------------------------------------------------------------
# Sanitizer flags
#
# INTERFACE hedefi sayesinde:
# 1. primerpair_core sanitizer ile derlenir.
# 2. Son executable sanitizer runtime ile linklenir.
# ------------------------------------------------------------

add_library(
    primerpair_sanitizers
    INTERFACE
)

if(PRIMERPAIR_ENABLE_SANITIZERS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(
            primerpair_sanitizers
            INTERFACE
                -O1
                -g
                -fsanitize=address,undefined
                -fno-omit-frame-pointer
        )

        target_link_options(
            primerpair_sanitizers
            INTERFACE
                -fsanitize=address,undefined
        )
    else()
        message(
            FATAL_ERROR
            "Sanitizers currently require GCC or Clang."
        )
    endif()
endif()

# ------------------------------------------------------------
# Core library
# ------------------------------------------------------------

add_library(
    primerpair_core
    STATIC
    src/fm_index.cpp
    src/suffix_array_builder.cpp
    src/bidirectional_fm_index.cpp
    src/ip_bwt_index.cpp
    src/batched_anchor_lookup.cpp
    src/batched_candidate_search.cpp
    src/hybrid_batched_primer_search.cpp
    src/hybrid_strand_aware_primer_search.cpp
    src/hybrid_primer_pair_search.cpp
    src/multiplex_primer_search.cpp
    src/global_multiplex_cross_join.cpp
    src/multiplex_primer_search_v2.cpp
    src/adaptive_multiplex_primer_search.cpp
    src/ip_bwt_rmi.cpp
    src/anchor_search.cpp
    src/approximate_anchor_search.cpp
    src/search_difficulty.cpp
    src/search_strategy.cpp
    src/packed_reference.cpp
    src/anchor_candidate_search.cpp
    src/single_primer_search.cpp
    src/strand_aware_primer_search.cpp
    src/primer_pair_search.cpp
    src/primer_pair_assembler.cpp
    src/mismatch_features.cpp
    src/biological_mismatch_profile.cpp
    src/thermodynamic_backend.cpp
    src/binding_site_thermodynamics.cpp
    src/combined_biological_features.cpp
    src/biological_risk_scoring.cpp
    src/pcr_chemistry_profile.cpp
    src/report_model.cpp
    src/report_serialization.cpp
    src/html_report.cpp
    src/annotated_report_pipeline.cpp
    src/primer_binding_site.cpp
    src/search_profile.cpp
    src/sensitive_primer_search.cpp
    src/sensitive_candidate_search.cpp
    src/sensitive_cost_estimator.cpp
    src/sensitive_adaptive_search.cpp
    src/sensitive_pair_constrained_search.cpp
    src/rank_scalar.cpp
    src/packed_bwt.cpp
)

target_include_directories(
    primerpair_core
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_compile_features(
    primerpair_core
    PUBLIC
        cxx_std_20
)

if(MSVC)
    target_compile_options(
        primerpair_core
        PRIVATE
            /W4
            /permissive-
    )
else()
    target_compile_options(
        primerpair_core
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
    )
endif()

# Sanitizer compile seçeneklerini core kütüphanesine uygular.
target_link_libraries(
    primerpair_core
    PRIVATE
        primerpair_sanitizers
)

# ------------------------------------------------------------
# Optional AVX2 backend
# ------------------------------------------------------------

if(PRIMERPAIR_ENABLE_AVX2)
    add_library(
        primerpair_avx2
        OBJECT
        src/rank_avx2.cpp
    )

    target_include_directories(
        primerpair_avx2
        PRIVATE
            ${PROJECT_SOURCE_DIR}/include
    )

    target_compile_features(
        primerpair_avx2
        PRIVATE
            cxx_std_20
    )

    if(MSVC)
        target_compile_options(
            primerpair_avx2
            PRIVATE
                /arch:AVX2
        )
    else()
        target_compile_options(
            primerpair_avx2
            PRIVATE
                -O3
                -mavx2
        )
    endif()

    target_sources(
        primerpair_core
        PRIVATE
            $<TARGET_OBJECTS:primerpair_avx2>
    )

    target_compile_definitions(
        primerpair_core
        PUBLIC
            PRIMERPAIR_HAVE_AVX2_BACKEND=1
    )
endif()

# ------------------------------------------------------------
# CLI executable
# ------------------------------------------------------------

add_executable(
    primerpair-bifm
    src/main.cpp
)

target_link_libraries(
    primerpair-bifm
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

set_target_properties(
    primerpair-bifm
    PROPERTIES
        OUTPUT_NAME oligoscout
)

enable_testing()
# ------------------------------------------------------------
# Exact FM-index tests
# ------------------------------------------------------------

add_executable(
    test_fm_index
    tests/test_fm_index.cpp
)

target_link_libraries(
    test_fm_index
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME exact_fm_index
    COMMAND test_fm_index
)
# ------------------------------------------------------------
# Checkpoint-rank tests
# ------------------------------------------------------------

add_executable(
    test_rank_support
    tests/test_rank_support.cpp
)

target_link_libraries(
    test_rank_support
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME checkpoint_rank
    COMMAND test_rank_support
)

# ------------------------------------------------------------
# Packed-BWT tests
# ------------------------------------------------------------

add_executable(
    test_packed_bwt
    tests/test_packed_bwt.cpp
)

target_link_libraries(
    test_packed_bwt
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME packed_bwt
    COMMAND test_packed_bwt
)

# ------------------------------------------------------------
# Packed checkpoint-rank equivalence tests
# ------------------------------------------------------------

add_executable(
    test_packed_rank_support
    tests/test_packed_rank_support.cpp
)

target_link_libraries(
    test_packed_rank_support
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME packed_checkpoint_rank
    COMMAND test_packed_rank_support
)
# ------------------------------------------------------------
# Bidirectional FM-index tests
# ------------------------------------------------------------

add_executable(
    test_bidirectional_fm_index
    tests/test_bidirectional_fm_index.cpp
)

target_link_libraries(
    test_bidirectional_fm_index
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME bidirectional_fm_index
    COMMAND test_bidirectional_fm_index
)

# ------------------------------------------------------------
# FM-index microbenchmark
# ------------------------------------------------------------

add_executable(
    benchmark_fm_index
    tests/benchmark_fm_index.cpp
)

target_link_libraries(
    benchmark_fm_index
    PRIVATE
        primerpair_core
)
# ------------------------------------------------------------
# FM-index build benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_index_build
    tests/benchmark_index_build.cpp
)

target_link_libraries(
    benchmark_index_build
    PRIVATE
        primerpair_core
)
# ------------------------------------------------------------
# Query scaling benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_query_scaling
    tests/benchmark_query_scaling.cpp
)

target_link_libraries(
    benchmark_query_scaling
    PRIVATE
        primerpair_core
)



# ------------------------------------------------------------
# Bidirectional FM-index extension benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_bidirectional_extension
    tests/benchmark_bidirectional_extension.cpp
)

target_link_libraries(
    benchmark_bidirectional_extension
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Exact 3-prime anchor search tests
# ------------------------------------------------------------

add_executable(
    test_anchor_search
    tests/test_anchor_search.cpp
)

target_link_libraries(
    test_anchor_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME anchor_search
    COMMAND test_anchor_search
)



# ------------------------------------------------------------
# 3-prime-aware approximate anchor tests
# ------------------------------------------------------------

add_executable(
    test_approximate_anchor_search
    tests/test_approximate_anchor_search.cpp
)

target_link_libraries(
    test_approximate_anchor_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME approximate_anchor_search
    COMMAND test_approximate_anchor_search
)



# ------------------------------------------------------------
# Approximate 3-prime anchor benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_approximate_anchor
    tests/benchmark_approximate_anchor.cpp
)

target_link_libraries(
    benchmark_approximate_anchor
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Real-FASTA approximate anchor benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_approximate_anchor_fasta
    tests/benchmark_approximate_anchor_fasta.cpp
)

target_link_libraries(
    benchmark_approximate_anchor_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Anchor difficulty profiler
# ------------------------------------------------------------

add_executable(
    profile_anchor_difficulty
    tests/profile_anchor_difficulty.cpp
)

target_link_libraries(
    profile_anchor_difficulty
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Search-difficulty bin benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_difficulty_bins
    tests/benchmark_difficulty_bins.cpp
)

target_link_libraries(
    benchmark_difficulty_bins
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Search difficulty tests
# ------------------------------------------------------------

add_executable(
    test_search_difficulty
    tests/test_search_difficulty.cpp
)

target_link_libraries(
    test_search_difficulty
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME search_difficulty
    COMMAND test_search_difficulty
)



# ------------------------------------------------------------
# Search strategy router tests
# ------------------------------------------------------------

add_executable(
    test_search_strategy
    tests/test_search_strategy.cpp
)

target_link_libraries(
    test_search_strategy
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME search_strategy
    COMMAND test_search_strategy
)



# ------------------------------------------------------------
# Packed reference tests
# ------------------------------------------------------------

add_executable(
    test_packed_reference
    tests/test_packed_reference.cpp
)

target_link_libraries(
    test_packed_reference
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME packed_reference
    COMMAND test_packed_reference
)



# ------------------------------------------------------------
# Anchor candidate verification tests
# ------------------------------------------------------------

add_executable(
    test_anchor_candidate_search
    tests/test_anchor_candidate_search.cpp
)

target_link_libraries(
    test_anchor_candidate_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME anchor_candidate_search
    COMMAND test_anchor_candidate_search
)



# ------------------------------------------------------------
# Candidate verification vs BiFM branching benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_candidate_vs_branching
    tests/benchmark_candidate_vs_branching.cpp
)

target_link_libraries(
    benchmark_candidate_vs_branching
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Candidate vs branching end-to-end benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_candidate_vs_branching_e2e
    tests/benchmark_candidate_vs_branching_e2e.cpp
)

target_link_libraries(
    benchmark_candidate_vs_branching_e2e
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Single-primer routed search tests
# ------------------------------------------------------------

add_executable(
    test_single_primer_search
    tests/test_single_primer_search.cpp
)

target_link_libraries(
    test_single_primer_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME single_primer_search
    COMMAND test_single_primer_search
)



# ------------------------------------------------------------
# Routed single-primer vs always-branching benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_routed_vs_branching
    tests/benchmark_routed_vs_branching.cpp
)

target_link_libraries(
    benchmark_routed_vs_branching
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Routed strategy stability benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_routed_stability
    tests/benchmark_routed_stability.cpp
)

target_link_libraries(
    benchmark_routed_stability
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Strand-aware primer search tests
# ------------------------------------------------------------

add_executable(
    test_strand_aware_primer_search
    tests/test_strand_aware_primer_search.cpp
)

target_link_libraries(
    test_strand_aware_primer_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME strand_aware_primer_search
    COMMAND test_strand_aware_primer_search
)



# ------------------------------------------------------------
# Primer-pair constrained search tests
# ------------------------------------------------------------

add_executable(
    test_primer_pair_search
    tests/test_primer_pair_search.cpp
)

target_link_libraries(
    test_primer_pair_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME primer_pair_search
    COMMAND test_primer_pair_search
)



# ------------------------------------------------------------
# Real-FASTA primer-pair constraint benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_pair_constraint_fasta
    tests/benchmark_pair_constraint_fasta.cpp
)

target_link_libraries(
    benchmark_pair_constraint_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Lower-bound vs sweep-line pairing kernel benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_pairing_kernel
    tests/benchmark_pairing_kernel.cpp
)

target_link_libraries(
    benchmark_pairing_kernel
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Mismatch positional feature tests
# ------------------------------------------------------------

add_executable(
    test_mismatch_features
    tests/test_mismatch_features.cpp
)

target_link_libraries(
    test_mismatch_features
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME mismatch_features
    COMMAND test_mismatch_features
)



# ------------------------------------------------------------
# Real-FASTA STRICT mismatch invariant validation
# ------------------------------------------------------------

add_executable(
    validate_strict_invariant_fasta
    tests/validate_strict_invariant_fasta.cpp
)

target_link_libraries(
    validate_strict_invariant_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Search profile façade tests
# ------------------------------------------------------------

add_executable(
    test_search_profile
    tests/test_search_profile.cpp
)

target_link_libraries(
    test_search_profile
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME search_profile
    COMMAND test_search_profile
)



# ------------------------------------------------------------
# SENSITIVE full-primer Hamming reference search tests
# ------------------------------------------------------------

add_executable(
    test_sensitive_primer_search
    tests/test_sensitive_primer_search.cpp
)

target_link_libraries(
    test_sensitive_primer_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME sensitive_primer_search
    COMMAND test_sensitive_primer_search
)



# ------------------------------------------------------------
# Real-FASTA STRICT vs SENSITIVE benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_strict_vs_sensitive_fasta
    tests/benchmark_strict_vs_sensitive_fasta.cpp
)

target_link_libraries(
    benchmark_strict_vs_sensitive_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Lossless SENSITIVE half-seed candidate tests
# ------------------------------------------------------------

add_executable(
    test_sensitive_candidate_search
    tests/test_sensitive_candidate_search.cpp
)

target_link_libraries(
    test_sensitive_candidate_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME sensitive_candidate_search
    COMMAND test_sensitive_candidate_search
)



# ------------------------------------------------------------
# Real-FASTA SENSITIVE reference vs candidate benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_sensitive_candidate_fasta
    tests/benchmark_sensitive_candidate_fasta.cpp
)

target_link_libraries(
    benchmark_sensitive_candidate_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# SENSITIVE candidate cost-estimator tests
# ------------------------------------------------------------

add_executable(
    test_sensitive_cost_estimator
    tests/test_sensitive_cost_estimator.cpp
)

target_link_libraries(
    test_sensitive_cost_estimator
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME sensitive_cost_estimator
    COMMAND test_sensitive_cost_estimator
)



# ------------------------------------------------------------
# Real-FASTA SENSITIVE adaptive-router training benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_sensitive_router_fasta
    tests/benchmark_sensitive_router_fasta.cpp
)

target_link_libraries(
    benchmark_sensitive_router_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Adaptive SENSITIVE backend tests
# ------------------------------------------------------------

add_executable(
    test_sensitive_adaptive_search
    tests/test_sensitive_adaptive_search.cpp
)

target_link_libraries(
    test_sensitive_adaptive_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME sensitive_adaptive_search
    COMMAND test_sensitive_adaptive_search
)



# ------------------------------------------------------------
# Independent chr22 adaptive SENSITIVE validation
# ------------------------------------------------------------

add_executable(
    benchmark_sensitive_adaptive_validation_fasta
    tests/benchmark_sensitive_adaptive_validation_fasta.cpp
)

target_link_libraries(
    benchmark_sensitive_adaptive_validation_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Search-profile SENSITIVE integration tests
# ------------------------------------------------------------

add_executable(
    test_search_profile_sensitive_integration
    tests/test_search_profile_sensitive_integration.cpp
)

target_link_libraries(
    test_search_profile_sensitive_integration
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME search_profile_sensitive_integration
    COMMAND test_search_profile_sensitive_integration
)



# ------------------------------------------------------------
# STRICT / SENSITIVE primer-pair profile integration
# ------------------------------------------------------------

add_executable(
    test_primer_pair_profile
    tests/test_primer_pair_profile.cpp
)

target_link_libraries(
    test_primer_pair_profile
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME primer_pair_profile
    COMMAND test_primer_pair_profile
)



# ------------------------------------------------------------
# STRICT vs SENSITIVE real-FASTA primer-pair benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_strict_vs_sensitive_pair_fasta
    tests/benchmark_strict_vs_sensitive_pair_fasta.cpp
)

target_link_libraries(
    benchmark_strict_vs_sensitive_pair_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Primer-pair runtime breakdown benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_pair_runtime_breakdown_fasta
    tests/benchmark_pair_runtime_breakdown_fasta.cpp
)

target_link_libraries(
    benchmark_pair_runtime_breakdown_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# ------------------------------------------------------------
# Pair-constrained SENSITIVE search
# ------------------------------------------------------------

add_executable(
    test_sensitive_pair_constrained_search
    tests/test_sensitive_pair_constrained_search.cpp
)

target_link_libraries(
    test_sensitive_pair_constrained_search
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)

add_test(
    NAME sensitive_pair_constrained_search
    COMMAND test_sensitive_pair_constrained_search
)



# ------------------------------------------------------------
# Global-global vs pair-constrained SENSITIVE benchmark
# ------------------------------------------------------------

add_executable(
    benchmark_sensitive_pair_constrained_fasta
    tests/benchmark_sensitive_pair_constrained_fasta.cpp
)

target_link_libraries(
    benchmark_sensitive_pair_constrained_fasta
    PRIVATE
        primerpair_core
        primerpair_sanitizers
)



# Anchor-selection oracle benchmark.
# Benchmark only: intentionally not registered with CTest.
add_executable(
    benchmark_anchor_oracle_fasta
    tests/benchmark_anchor_oracle_fasta.cpp
)

target_link_libraries(
    benchmark_anchor_oracle_fasta
    PRIVATE
        primerpair_core
)


# Suffix-array sampling-rate benchmark.
# Benchmark only; not registered with CTest.
add_executable(
    benchmark_sa_sampling_fasta
    tests/benchmark_sa_sampling_fasta.cpp
)

target_link_libraries(
    benchmark_sa_sampling_fasta
    PRIVATE
        primerpair_core
)


add_executable(
    test_sensitive_pair_constrained_lengths
    tests/test_sensitive_pair_constrained_lengths.cpp
)

target_link_libraries(
    test_sensitive_pair_constrained_lengths
    PRIVATE
        primerpair_core
)

add_test(
    NAME sensitive_pair_constrained_lengths
    COMMAND test_sensitive_pair_constrained_lengths
)


# 18–35 nt k3 Candidate-vs-Exhaustive router benchmark.
# Benchmark only; intentionally not registered with CTest.
add_executable(
    benchmark_sensitive_length_router_fasta
    tests/benchmark_sensitive_length_router_fasta.cpp
)

target_link_libraries(
    benchmark_sensitive_length_router_fasta
    PRIVATE
        primerpair_core
)


# Repeat-rich real-genome stress benchmark for the
# 23..35 direct-Candidate adaptive routing policy.
# Benchmark only; not registered with CTest.
add_executable(
    benchmark_repeat_rich_router_fasta
    tests/benchmark_repeat_rich_router_fasta.cpp
)

target_link_libraries(
    benchmark_repeat_rich_router_fasta
    PRIVATE
        primerpair_core
)


# Whole-genome chromosome/shard orchestration.
target_sources(
    primerpair_core
    PRIVATE
        src/genome_search.cpp
    src/ppfm_io.cpp
    src/ppfm_manifest.cpp
    src/ppfm_shard_cache.cpp
    src/persistent_genome_search.cpp
)


add_executable(
    test_genome_search
    tests/test_genome_search.cpp
)

target_link_libraries(
    test_genome_search
    PRIVATE
        primerpair_core
)

add_test(
    NAME genome_search
    COMMAND test_genome_search
)


# Reusable single/multi-record FASTA parser.
target_sources(
    primerpair_core
    PRIVATE
        src/fasta_reader.cpp
)


add_executable(
    test_fasta_reader
    tests/test_fasta_reader.cpp
)

target_link_libraries(
    test_fasta_reader
    PRIVATE
        primerpair_core
)

add_test(
    NAME fasta_reader
    COMMAND test_fasta_reader
)


# GRCh38 primary chromosome metadata resolver.
target_sources(
    primerpair_core
    PRIVATE
        src/grch38_reference.cpp
)


add_executable(
    test_grch38_reference
    tests/test_grch38_reference.cpp
)

target_link_libraries(
    test_grch38_reference
    PRIVATE
        primerpair_core
)

add_test(
    NAME grch38_reference
    COMMAND test_grch38_reference
)


# Selective streaming GRCh38 chromosome loader.
target_sources(
    primerpair_core
    PRIVATE
        src/grch38_genome_loader.cpp
)


add_executable(
    test_grch38_genome_loader
    tests/test_grch38_genome_loader.cpp
)

target_link_libraries(
    test_grch38_genome_loader
    PRIVATE
        primerpair_core
)

add_test(
    NAME grch38_genome_loader
    COMMAND test_grch38_genome_loader
)


# Real GRCh38 chr21+chr22 in-memory shard benchmark.
# Benchmark only; not registered with CTest.
add_executable(
    benchmark_grch38_two_shards
    tests/benchmark_grch38_two_shards.cpp
)

target_link_libraries(
    benchmark_grch38_two_shards
    PRIVATE
        primerpair_core
)


add_executable(
    test_ppfm_io
    tests/test_ppfm_io.cpp
)

target_link_libraries(
    test_ppfm_io
    PRIVATE
        primerpair_core
)

add_test(
    NAME ppfm_io
    COMMAND test_ppfm_io
)


add_executable(
    build_grch38_ppfm
    tools/build_grch38_ppfm.cpp
)

target_link_libraries(
    build_grch38_ppfm
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_ppfm_load
    tools/benchmark_ppfm_load.cpp
)

target_link_libraries(
    benchmark_ppfm_load
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_ppfm_two_shards_pair
    tools/benchmark_ppfm_two_shards_pair.cpp
)

target_link_libraries(
    benchmark_ppfm_two_shards_pair
    PRIVATE
        primerpair_core
)


add_executable(
    test_genome_search_ppfm
    tests/test_genome_search_ppfm.cpp
)

target_link_libraries(
    test_genome_search_ppfm
    PRIVATE
        primerpair_core
)

add_test(
    NAME genome_search_ppfm
    COMMAND test_genome_search_ppfm
)


add_executable(
    benchmark_genome_ppfm_two_shards
    tools/benchmark_genome_ppfm_two_shards.cpp
)

target_link_libraries(
    benchmark_genome_ppfm_two_shards
    PRIVATE
        primerpair_core
)


add_executable(
    test_ppfm_manifest
    tests/test_ppfm_manifest.cpp
)

target_link_libraries(
    test_ppfm_manifest
    PRIVATE
        primerpair_core
)

add_test(
    NAME ppfm_manifest
    COMMAND test_ppfm_manifest
)

set_tests_properties(
    ppfm_manifest
    PROPERTIES
        WORKING_DIRECTORY
            "${CMAKE_SOURCE_DIR}"
)


add_executable(
    test_ppfm_shard_cache
    tests/test_ppfm_shard_cache.cpp
)

target_link_libraries(
    test_ppfm_shard_cache
    PRIVATE
        primerpair_core
)

add_test(
    NAME ppfm_shard_cache
    COMMAND test_ppfm_shard_cache
)


add_executable(
    test_persistent_genome_search
    tests/test_persistent_genome_search.cpp
)

target_link_libraries(
    test_persistent_genome_search
    PRIVATE
        primerpair_core
)

add_test(
    NAME persistent_genome_search
    COMMAND test_persistent_genome_search
)


add_executable(
    benchmark_grch38_persistent_pair
    tools/benchmark_grch38_persistent_pair.cpp
)

target_link_libraries(
    benchmark_grch38_persistent_pair
    PRIVATE
        primerpair_core
)


add_executable(
    test_primer_binding_site
    tests/test_primer_binding_site.cpp
)

target_link_libraries(
    test_primer_binding_site
    PRIVATE
        primerpair_core
)

add_test(
    NAME primer_binding_site
    COMMAND test_primer_binding_site
)


add_executable(
    report_grch38_binding_sites
    tools/report_grch38_binding_sites.cpp
)

target_link_libraries(
    report_grch38_binding_sites
    PRIVATE
        primerpair_core
)


add_executable(
    test_ip_bwt_index
    tests/test_ip_bwt_index.cpp
)

target_link_libraries(
    test_ip_bwt_index
    PRIVATE
        primerpair_core
)

add_test(
    NAME ip_bwt_index
    COMMAND test_ip_bwt_index
)


add_executable(
    test_suffix_array_builder
    tests/test_suffix_array_builder.cpp
)

target_link_libraries(
    test_suffix_array_builder
    PRIVATE
        primerpair_core
)

add_test(
    NAME suffix_array_builder
    COMMAND test_suffix_array_builder
)


add_executable(
    test_ip_bwt_rmi
    tests/test_ip_bwt_rmi.cpp
)

target_link_libraries(
    test_ip_bwt_rmi
    PRIVATE
        primerpair_core
)

add_test(
    NAME ip_bwt_rmi
    COMMAND test_ip_bwt_rmi
)


add_executable(
    test_ip_bwt_rmi_primer_lengths
    tests/test_ip_bwt_rmi_primer_lengths.cpp
)

target_link_libraries(
    test_ip_bwt_rmi_primer_lengths
    PRIVATE
        primerpair_core
)

add_test(
    NAME ip_bwt_rmi_primer_lengths
    COMMAND test_ip_bwt_rmi_primer_lengths
)


add_executable(
    diagnose_fm_backward_search
    tools/diagnose_fm_backward_search.cpp
)

target_link_libraries(
    diagnose_fm_backward_search
    PRIVATE
        primerpair_core
)


add_executable(
    diagnose_ip_bwt_partial
    tools/diagnose_ip_bwt_partial.cpp
)

target_link_libraries(
    diagnose_ip_bwt_partial
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_ip_bwt_rmi
    tools/benchmark_ip_bwt_rmi.cpp
)

target_link_libraries(
    benchmark_ip_bwt_rmi
    PRIVATE
        primerpair_core
)


add_executable(
    test_ip_bwt_batch4
    tests/test_ip_bwt_batch4.cpp
)

target_link_libraries(
    test_ip_bwt_batch4
    PRIVATE
        primerpair_core
)

add_test(
    NAME ip_bwt_batch4
    COMMAND test_ip_bwt_batch4
)


add_executable(
    generate_neural_window_dataset
    tools/generate_neural_window_dataset.cpp
)

target_link_libraries(
    generate_neural_window_dataset
    PRIVATE
        primerpair_core
)

target_sources(
    primerpair_core
    PRIVATE
        src/neural_window_predictor.cpp
)

add_executable(
    test_neural_window_predictor
    tests/test_neural_window_predictor.cpp
)

target_link_libraries(
    test_neural_window_predictor
    PRIVATE
        primerpair_core
)

add_test(
    NAME neural_window_predictor
    COMMAND test_neural_window_predictor
)

add_executable(
    test_ip_bwt_neural_window
    tests/test_ip_bwt_neural_window.cpp
)

target_link_libraries(
    test_ip_bwt_neural_window
    PRIVATE
        primerpair_core
)

add_test(
    NAME ip_bwt_neural_window
    COMMAND test_ip_bwt_neural_window
)

add_executable(
    benchmark_neural_certified_ip_bwt
    tools/benchmark_neural_certified_ip_bwt.cpp
)

target_link_libraries(
    benchmark_neural_certified_ip_bwt
    PRIVATE
        primerpair_core
)

add_executable(
    test_neural_window_batch8
    tests/test_neural_window_batch8.cpp
)

target_link_libraries(
    test_neural_window_batch8
    PRIVATE
        primerpair_core
)

add_test(
    NAME neural_window_batch8
    COMMAND test_neural_window_batch8
)

add_executable(
    benchmark_neural_batch8
    tools/benchmark_neural_batch8.cpp
)

target_link_libraries(
    benchmark_neural_batch8
    PRIVATE
        primerpair_core
)


# BEGIN NEURAL WINDOW AVX2 FLAGS

option(
    PRIMERPAIR_ENABLE_NEURAL_AVX2
    "Enable AVX2 batch neural-window inference"
    OFF
)

if(PRIMERPAIR_ENABLE_NEURAL_AVX2)

    if(
        CMAKE_CXX_COMPILER_ID
        MATCHES
        "GNU|Clang"
    )

        # Only the neural predictor translation unit
        # receives AVX2 ISA instructions.
        #
        # The rest of primerpair_core remains on the
        # normal release ISA so the binary-IPBWT
        # baseline is not silently changed.
        set_source_files_properties(
            src/neural_window_predictor.cpp
            PROPERTIES
                COMPILE_OPTIONS "-mavx2"
        )

        target_compile_options(
            test_neural_window_batch8
            PRIVATE
                -mavx2
        )

        target_compile_options(
            benchmark_neural_batch8
            PRIVATE
                -mavx2
        )

    else()

        message(
            FATAL_ERROR
            "Neural AVX2 benchmark currently supports "
            "GNU/Clang compilers."
        )

    endif()

endif()

# END NEURAL WINDOW AVX2 FLAGS

add_executable(
    test_ip_bwt_neural_window_batch8
    tests/test_ip_bwt_neural_window_batch8.cpp
)

target_link_libraries(
    test_ip_bwt_neural_window_batch8
    PRIVATE
        primerpair_core
)

add_test(
    NAME ip_bwt_neural_window_batch8
    COMMAND test_ip_bwt_neural_window_batch8
)


add_executable(
    benchmark_cdf_actual
    tools/benchmark_cdf_actual.cpp
)

target_link_libraries(
    benchmark_cdf_actual
    PRIVATE
        primerpair_core
)

if(
    CMAKE_CXX_COMPILER_ID
    MATCHES
    "GNU|Clang"
)
    target_compile_options(
        benchmark_cdf_actual
        PRIVATE
            -mavx2
    )
endif()


add_executable(
    test_ip_bwt_many
    tests/test_ip_bwt_many.cpp
)

target_link_libraries(
    test_ip_bwt_many
    PRIVATE
        primerpair_core
)

add_test(
    NAME ip_bwt_many
    COMMAND test_ip_bwt_many
)


add_executable(
    benchmark_ip_bwt_many
    tools/benchmark_ip_bwt_many.cpp
)

target_link_libraries(
    benchmark_ip_bwt_many
    PRIVATE
        primerpair_core
)


add_executable(
    test_ip_bwt_bifm_anchor_equivalence
    tests/test_ip_bwt_bifm_anchor_equivalence.cpp
)

target_link_libraries(
    test_ip_bwt_bifm_anchor_equivalence
    PRIVATE
        primerpair_core
)

add_test(
    NAME ip_bwt_bifm_anchor_equivalence
    COMMAND test_ip_bwt_bifm_anchor_equivalence
)


add_executable(
    test_batched_anchor_lookup
    tests/test_batched_anchor_lookup.cpp
)

target_link_libraries(
    test_batched_anchor_lookup
    PRIVATE
        primerpair_core
)

add_test(
    NAME batched_anchor_lookup
    COMMAND test_batched_anchor_lookup
)


add_executable(
    test_anchor_candidate_ipbwt_equivalence
    tests/test_anchor_candidate_ipbwt_equivalence.cpp
)

target_link_libraries(
    test_anchor_candidate_ipbwt_equivalence
    PRIVATE
        primerpair_core
)

add_test(
    NAME anchor_candidate_ipbwt_equivalence
    COMMAND test_anchor_candidate_ipbwt_equivalence
)


add_executable(
    test_batched_candidate_search
    tests/test_batched_candidate_search.cpp
)

target_link_libraries(
    test_batched_candidate_search
    PRIVATE
        primerpair_core
)

add_test(
    NAME batched_candidate_search
    COMMAND test_batched_candidate_search
)


add_executable(
    benchmark_batched_candidate
    tools/benchmark_batched_candidate.cpp
)

target_link_libraries(
    benchmark_batched_candidate
    PRIVATE
        primerpair_core
)


add_executable(
    test_hybrid_batched_primer_search
    tests/test_hybrid_batched_primer_search.cpp
)

target_link_libraries(
    test_hybrid_batched_primer_search
    PRIVATE
        primerpair_core
)

add_test(
    NAME hybrid_batched_primer_search
    COMMAND test_hybrid_batched_primer_search
)


add_executable(
    benchmark_hybrid_batched_primer
    tools/benchmark_hybrid_batched_primer.cpp
)

target_link_libraries(
    benchmark_hybrid_batched_primer
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_hybrid_route_breakdown
    tools/benchmark_hybrid_route_breakdown.cpp
)

target_link_libraries(
    benchmark_hybrid_route_breakdown
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_branching_profile
    tools/benchmark_branching_profile.cpp
)

target_link_libraries(
    benchmark_branching_profile
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_branching_timing
    tools/benchmark_branching_timing.cpp
)

target_link_libraries(
    benchmark_branching_timing
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_sa_sample_rate
    tools/benchmark_sa_sample_rate.cpp
)

target_link_libraries(
    benchmark_sa_sample_rate
    PRIVATE
        primerpair_core
)


add_executable(
    test_rank_all
    tests/test_rank_all.cpp
)

target_link_libraries(
    test_rank_all
    PRIVATE
        primerpair_core
)

add_test(
    NAME rank_all
    COMMAND test_rank_all
)


add_executable(
    test_extend_left_all
    tests/test_extend_left_all.cpp
)

target_link_libraries(
    test_extend_left_all
    PRIVATE
        primerpair_core
)

add_test(
    NAME extend_left_all
    COMMAND test_extend_left_all
)


add_executable(
    test_locate_unsorted
    tests/test_locate_unsorted.cpp
)

target_link_libraries(
    test_locate_unsorted
    PRIVATE
        primerpair_core
)

add_test(
    NAME locate_unsorted
    COMMAND test_locate_unsorted
)


add_executable(
    test_lf_symbol_id
    tests/test_lf_symbol_id.cpp
)

target_link_libraries(
    test_lf_symbol_id
    PRIVATE
        primerpair_core
)

add_test(
    NAME lf_symbol_id
    COMMAND test_lf_symbol_id
)


add_executable(
    test_hybrid_strand_aware_primer_search
    tests/test_hybrid_strand_aware_primer_search.cpp
)

target_link_libraries(
    test_hybrid_strand_aware_primer_search
    PRIVATE
        primerpair_core
)

add_test(
    NAME hybrid_strand_aware_primer_search
    COMMAND test_hybrid_strand_aware_primer_search
)


add_executable(
    test_hybrid_strand_equivalence
    tests/test_hybrid_strand_equivalence.cpp
)

target_link_libraries(
    test_hybrid_strand_equivalence
    PRIVATE
        primerpair_core
)

add_test(
    NAME hybrid_strand_equivalence
    COMMAND test_hybrid_strand_equivalence
)


add_executable(
    test_primer_pair_assembler_equivalence
    tests/test_primer_pair_assembler_equivalence.cpp
)

target_link_libraries(
    test_primer_pair_assembler_equivalence
    PRIVATE
        primerpair_core
)

add_test(
    NAME primer_pair_assembler_equivalence
    COMMAND test_primer_pair_assembler_equivalence
)


add_executable(
    test_hybrid_primer_pair_equivalence
    tests/test_hybrid_primer_pair_equivalence.cpp
)

target_link_libraries(
    test_hybrid_primer_pair_equivalence
    PRIVATE
        primerpair_core
)

add_test(
    NAME hybrid_primer_pair_equivalence
    COMMAND test_hybrid_primer_pair_equivalence
)


add_executable(
    test_multiplex_primer_search
    tests/test_multiplex_primer_search.cpp
)

target_link_libraries(
    test_multiplex_primer_search
    PRIVATE
        primerpair_core
)

add_test(
    NAME multiplex_primer_search
    COMMAND test_multiplex_primer_search
)


add_executable(
    test_multiplex_cross_pair_equivalence
    tests/test_multiplex_cross_pair_equivalence.cpp
)

target_link_libraries(
    test_multiplex_cross_pair_equivalence
    PRIVATE
        primerpair_core
)

add_test(
    NAME multiplex_cross_pair_equivalence
    COMMAND test_multiplex_cross_pair_equivalence
)


add_executable(
    test_multiplex_cross_join_cache
    tests/test_multiplex_cross_join_cache.cpp
)

target_link_libraries(
    test_multiplex_cross_join_cache
    PRIVATE
        primerpair_core
)

add_test(
    NAME multiplex_cross_join_cache
    COMMAND test_multiplex_cross_join_cache
)


add_executable(
    benchmark_multiplex_reuse
    tools/benchmark_multiplex_reuse.cpp
)

target_link_libraries(
    benchmark_multiplex_reuse
    PRIVATE
        primerpair_core
)


add_executable(
    test_global_multiplex_cross_join
    tests/test_global_multiplex_cross_join.cpp
)

target_link_libraries(
    test_global_multiplex_cross_join
    PRIVATE
        primerpair_core
)

add_test(
    NAME global_multiplex_cross_join
    COMMAND test_global_multiplex_cross_join
)


add_executable(
    test_multiplex_v2_equivalence
    tests/test_multiplex_v2_equivalence.cpp
)

target_link_libraries(
    test_multiplex_v2_equivalence
    PRIVATE
        primerpair_core
)

add_test(
    NAME multiplex_v2_equivalence
    COMMAND test_multiplex_v2_equivalence
)


add_executable(
    test_multiplex_v2_repeated_primer_equivalence
    tests/test_multiplex_v2_repeated_primer_equivalence.cpp
)

target_link_libraries(
    test_multiplex_v2_repeated_primer_equivalence
    PRIVATE
        primerpair_core
)

add_test(
    NAME multiplex_v2_repeated_primer_equivalence
    COMMAND test_multiplex_v2_repeated_primer_equivalence
)


add_executable(
    benchmark_multiplex_v2_scaling
    tools/benchmark_multiplex_v2_scaling.cpp
)

target_link_libraries(
    benchmark_multiplex_v2_scaling
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_grch38_chr22_multiplex_v2
    tools/benchmark_grch38_chr22_multiplex_v2.cpp
)

target_link_libraries(
    benchmark_grch38_chr22_multiplex_v2
    PRIVATE
        primerpair_core
)


add_executable(
    benchmark_grch38_shard_multiplex_v2
    tools/benchmark_grch38_shard_multiplex_v2.cpp
)

target_link_libraries(
    benchmark_grch38_shard_multiplex_v2
    PRIVATE
        primerpair_core
)


# Persistent whole-genome Multiplex V2 production engine.
target_sources(
    primerpair_core
    PRIVATE
        src/persistent_multiplex_primer_search_v2.cpp
)


add_executable(
    validate_persistent_multiplex_v2_grch38
    tools/validate_persistent_multiplex_v2_grch38.cpp
)

target_link_libraries(
    validate_persistent_multiplex_v2_grch38
    PRIVATE
        primerpair_core
)


# Test #64:
# Persistent whole-genome Multiplex V2.
#
# Verifies:
# - chromosome-local assembly
# - no cross-chromosome PCR products
# - persistent PPFM loading
# - capacity-one LRU eviction
# - production V1 local oracle equivalence
add_executable(
    test_persistent_multiplex_v2
    tests/test_persistent_multiplex_v2.cpp
)

target_link_libraries(
    test_persistent_multiplex_v2
    PRIVATE
        primerpair_core
)

add_test(
    NAME persistent_multiplex_v2
    COMMAND test_persistent_multiplex_v2
)


add_executable(
    profile_persistent_multiplex_v2_memory
    tools/profile_persistent_multiplex_v2_memory.cpp
)

target_link_libraries(
    profile_persistent_multiplex_v2_memory
    PRIVATE
        primerpair_core
)


# Test #65:
# PPFM cache evict-before-load memory policy.
add_executable(
    test_ppfm_shard_cache_evict_before_load
    tests/test_ppfm_shard_cache_evict_before_load.cpp
)

target_link_libraries(
    test_ppfm_shard_cache_evict_before_load
    PRIVATE
        primerpair_core
)

add_test(
    NAME ppfm_shard_cache_evict_before_load
    COMMAND test_ppfm_shard_cache_evict_before_load
)


add_executable(
    profile_persistent_multiplex_v2_repeated
    tools/profile_persistent_multiplex_v2_repeated.cpp
)

target_link_libraries(
    profile_persistent_multiplex_v2_repeated
    PRIVATE
        primerpair_core
)


# Test #66:
# Persistent Multiplex V2 shard-major search_many().
add_executable(
    test_persistent_multiplex_v2_search_many
    tests/test_persistent_multiplex_v2_search_many.cpp
)

target_link_libraries(
    test_persistent_multiplex_v2_search_many
    PRIVATE
        primerpair_core
)

add_test(
    NAME persistent_multiplex_v2_search_many
    COMMAND test_persistent_multiplex_v2_search_many
)


add_executable(
    profile_persistent_multiplex_v2_search_many
    tools/profile_persistent_multiplex_v2_search_many.cpp
)

target_link_libraries(
    profile_persistent_multiplex_v2_search_many
    PRIVATE
        primerpair_core
)


# Test #67:
# Empirical Production V1 / Global Sweep V2
# multiplex crossover strategy router.
add_executable(
    test_multiplex_strategy_router
    tests/test_multiplex_strategy_router.cpp
)

target_link_libraries(
    test_multiplex_strategy_router
    PRIVATE
        primerpair_core
)

add_test(
    NAME multiplex_strategy_router
    COMMAND test_multiplex_strategy_router
)


# Test #68:
# Adaptive production multiplex execution.
#
# Validates:
# - conservative 12-pair production crossover
# - 11-pair Production V1 route
# - 12-pair Global Sweep V2 route
# - exact selected-backend result equivalence
# - custom router override
add_executable(
    test_adaptive_multiplex_primer_search
    tests/test_adaptive_multiplex_primer_search.cpp
)

target_link_libraries(
    test_adaptive_multiplex_primer_search
    PRIVATE
        primerpair_core
)

add_test(
    NAME adaptive_multiplex_primer_search
    COMMAND test_adaptive_multiplex_primer_search
)


# ----------------------------------------------------------
# Biological mismatch positional profile — Test #69
# ----------------------------------------------------------
add_executable(
    test_biological_mismatch_profile
    tests/test_biological_mismatch_profile.cpp
)

target_link_libraries(
    test_biological_mismatch_profile
    PRIVATE
    primerpair_core
)

add_test(
    NAME biological_mismatch_profile
    COMMAND test_biological_mismatch_profile
)


# ----------------------------------------------------------
# Primer3 Thermodynamic Backend V1 — Test #70
# ----------------------------------------------------------
add_executable(
    test_thermodynamic_backend
    tests/test_thermodynamic_backend.cpp
)

target_link_libraries(
    test_thermodynamic_backend
    PRIVATE
    primerpair_core
)

add_test(
    NAME thermodynamic_backend
    COMMAND test_thermodynamic_backend
)

set_tests_properties(
    thermodynamic_backend
    PROPERTIES
    ENVIRONMENT
    "PRIMERPAIR_NTTHAL=${CMAKE_SOURCE_DIR}/third_party/primer3/src/ntthal;PRIMERPAIR_OLIGOTM=${CMAKE_SOURCE_DIR}/third_party/primer3/src/oligotm;PRIMERPAIR_THERMO_CONFIG=${CMAKE_SOURCE_DIR}/third_party/primer3/src/primer3_config"
)


# ----------------------------------------------------------
# Binding-Site Thermodynamics V1 — Test #71
# ----------------------------------------------------------
add_executable(
    test_binding_site_thermodynamics
    tests/test_binding_site_thermodynamics.cpp
)

target_link_libraries(
    test_binding_site_thermodynamics
    PRIVATE
    primerpair_core
)

add_test(
    NAME binding_site_thermodynamics
    COMMAND test_binding_site_thermodynamics
)

set_tests_properties(
    binding_site_thermodynamics
    PROPERTIES
    ENVIRONMENT
    "PRIMERPAIR_NTTHAL=${CMAKE_SOURCE_DIR}/third_party/primer3/src/ntthal;PRIMERPAIR_OLIGOTM=${CMAKE_SOURCE_DIR}/third_party/primer3/src/oligotm;PRIMERPAIR_THERMO_CONFIG=${CMAKE_SOURCE_DIR}/third_party/primer3/src/primer3_config"
)


# ----------------------------------------------------------
# Combined Biological Feature Vector V1 — Test #72
# ----------------------------------------------------------
add_executable(
    test_combined_biological_features
    tests/test_combined_biological_features.cpp
)

target_link_libraries(
    test_combined_biological_features
    PRIVATE
    primerpair_core
)

add_test(
    NAME combined_biological_features
    COMMAND test_combined_biological_features
)

set_tests_properties(
    combined_biological_features
    PROPERTIES
    ENVIRONMENT
    "PRIMERPAIR_NTTHAL=${CMAKE_SOURCE_DIR}/third_party/primer3/src/ntthal;PRIMERPAIR_OLIGOTM=${CMAKE_SOURCE_DIR}/third_party/primer3/src/oligotm;PRIMERPAIR_THERMO_CONFIG=${CMAKE_SOURCE_DIR}/third_party/primer3/src/primer3_config"
)


# ----------------------------------------------------------
# Biological Risk Scoring Interface V1 — Test #73
# ----------------------------------------------------------
add_executable(
    test_biological_risk_scoring
    tests/test_biological_risk_scoring.cpp
)

target_link_libraries(
    test_biological_risk_scoring
    PRIVATE
    primerpair_core
)

add_test(
    NAME biological_risk_scoring
    COMMAND test_biological_risk_scoring
)


# ----------------------------------------------------------
# PCR Chemistry Profile V1 — Test #74
# ----------------------------------------------------------
add_executable(
    test_pcr_chemistry_profile
    tests/test_pcr_chemistry_profile.cpp
)

target_link_libraries(
    test_pcr_chemistry_profile
    PRIVATE
    primerpair_core
)

add_test(
    NAME pcr_chemistry_profile
    COMMAND test_pcr_chemistry_profile
)


# ----------------------------------------------------------
# Structured Result / Report Model V1 — Test #75
# ----------------------------------------------------------
add_executable(
    test_report_model
    tests/test_report_model.cpp
)

target_link_libraries(
    test_report_model
    PRIVATE
    primerpair_core
)

add_test(
    NAME report_model
    COMMAND test_report_model
)


# ----------------------------------------------------------
# TSV + JSON Report Serialization V1 — Test #76
# ----------------------------------------------------------
add_executable(
    test_report_serialization
    tests/test_report_serialization.cpp
)

target_link_libraries(
    test_report_serialization
    PRIVATE
    primerpair_core
)

add_test(
    NAME report_serialization
    COMMAND test_report_serialization
)


# ----------------------------------------------------------
# Self-Contained HTML Report V1 — Test #77
# ----------------------------------------------------------
add_executable(
    test_html_report
    tests/test_html_report.cpp
)

target_link_libraries(
    test_html_report
    PRIVATE
    primerpair_core
)

add_test(
    NAME html_report
    COMMAND test_html_report
)


# ----------------------------------------------------------
# Annotated Report Pipeline V1 — Test #78
# ----------------------------------------------------------
add_executable(
    test_annotated_report_pipeline
    tests/test_annotated_report_pipeline.cpp
)

target_link_libraries(
    test_annotated_report_pipeline
    PRIVATE
    primerpair_core
)

add_test(
    NAME annotated_report_pipeline
    COMMAND test_annotated_report_pipeline
)

set_tests_properties(
    annotated_report_pipeline
    PROPERTIES
    ENVIRONMENT
    "PRIMERPAIR_NTTHAL=${CMAKE_SOURCE_DIR}/third_party/primer3/src/ntthal;PRIMERPAIR_OLIGOTM=${CMAKE_SOURCE_DIR}/third_party/primer3/src/oligotm;PRIMERPAIR_THERMO_CONFIG=${CMAKE_SOURCE_DIR}/third_party/primer3/src/primer3_config"
)


# ----------------------------------------------------------
# Whole-GRCh38 Annotated Report Integration Tool
# ----------------------------------------------------------
add_executable(
    integrate_grch38_annotated_report
    tools/integrate_grch38_annotated_report.cpp
)

target_link_libraries(
    integrate_grch38_annotated_report
    PRIVATE
    primerpair_core
)

# ------------------------------------------------------------
# Test #79B full PrimerPair product export tool
# ------------------------------------------------------------

add_executable(
    export_test79_primerpair_products
    tools/export_test79_primerpair_products.cpp
)

target_link_libraries(
    export_test79_primerpair_products
    PRIVATE
    primerpair_core
)
