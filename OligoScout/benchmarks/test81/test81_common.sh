#!/usr/bin/env bash
set -euo pipefail

TEST81_ARTIFACT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
TEST81_REPO="${TEST81_REPO:-/home/nilay/primerpair-search}"
TEST81_E_ROOT="${TEST81_OUTPUT_ROOT:-/mnt/e/OligoScout_Benchmarks/test81_publication}"

REFERENCE="$TEST81_REPO/data/indexes/external_benchmarks/mfeprimer_4.5.1/grch38_p14_primary24/grch38_p14_primary24.fa"
REFERENCE_FAI="${REFERENCE}.fai"
MANIFEST="$TEST81_REPO/results/benchmarks/grch38_ppfm_batch/manifest.sa8.tsv"
PPFM_DIR="$TEST81_REPO/data/indexes/ppfm"
LEGACY_PANEL="$TEST81_REPO/data/benchmarks/test79_inputs/panel32_pairs.tsv"
CORE_LIBRARY="$TEST81_REPO/build/release/libprimerpair_core.a"
PRIMER3="$TEST81_REPO/third_party/primer3/src/primer3_core"
JELLYFISH="$(command -v jellyfish || true)"
BLASTN="$(command -v blastn || true)"
BLAST_DB="$TEST81_REPO/data/indexes/external_benchmarks/blast_2.16.0/grch38_p14_primary24/grch38_p14_primary24"
BOWTIE="$TEST81_REPO/data/tools/external_benchmarks/test80_multi_aligner/bowtie-1.3.1/bowtie"
BOWTIE_INDEX="$TEST81_REPO/data/indexes/external_benchmarks/test80_multi_aligner/bowtie1_1.3.1/grch38_primary24"
BWA="$TEST81_REPO/data/tools/external_benchmarks/test80_multi_aligner/bwa-0.7.19/bwa"
BWA_INDEX="$TEST81_REPO/data/indexes/external_benchmarks/test80_multi_aligner/bwa_0.7.19/grch38_primary24"
PIPELINE="$TEST81_ARTIFACT_DIR/test81_pipeline.py"
PANEL_GENERATOR="$TEST81_ARTIFACT_DIR/test81_generate_panel.py"

banner() {
    printf '\n============================================================\n%s\n============================================================\n' "$1"
}

require_file() {
    [[ -f "$1" ]] || { echo "ERROR: required file missing: $1" >&2; exit 1; }
}

test81_preflight() {
    banner "TEST #81 — E: DRIVE PREFLIGHT"
    [[ -d /mnt/e ]] || {
        echo "ERROR: E: is not mounted as /mnt/e in WSL." >&2
        echo "Open PowerShell and verify the E: drive, then restart WSL." >&2
        exit 1
    }
    mkdir -p "$TEST81_E_ROOT"
    local probe="$TEST81_E_ROOT/.write_probe_$$"
    printf 'test81\n' > "$probe"
    rm -f -- "$probe"
    local free_kib
    free_kib="$(df -Pk /mnt/e | awk 'NR==2 {print $4}')"
    (( free_kib >= 40 * 1024 * 1024 )) || {
        echo "ERROR: at least 40 GiB free is required on E:; available KiB: $free_kib" >&2
        exit 1
    }
    for path in "$REFERENCE" "$REFERENCE_FAI" "$MANIFEST" "$LEGACY_PANEL" \
        "$CORE_LIBRARY" "$PRIMER3" "$BOWTIE" "$BWA" "$PIPELINE" \
        "$PANEL_GENERATOR" "$TEST81_ARTIFACT_DIR/src/test81_anchor_oracle.cpp" \
        "$TEST81_ARTIFACT_DIR/src/benchmark_test81_panel.cpp"; do
        require_file "$path"
    done
    [[ -n "$JELLYFISH" ]] || { echo "ERROR: jellyfish not found" >&2; exit 1; }
    [[ -n "$BLASTN" ]] || { echo "ERROR: blastn not found" >&2; exit 1; }
    require_file "${BLAST_DB}.nsq"
    require_file "${BOWTIE_INDEX}.1.ebwt"
    require_file "${BWA_INDEX}.bwt"
    [[ "$(find "$PPFM_DIR" -maxdepth 1 -type f -name '*.ppfm' | wc -l)" -eq 24 ]] || {
        echo "ERROR: expected 24 PPFM shards in $PPFM_DIR" >&2
        exit 1
    }
    echo -e "E_DRIVE_FREE_GIB\t$((free_kib / 1024 / 1024))"
    echo -e "OUTPUT_ROOT\t$TEST81_E_ROOT"
    echo -e "PREFLIGHT\tPASS"
}

test81_initialize_run() {
    mkdir -p "$TEST81_E_ROOT/runs" "$TEST81_E_ROOT/shared"
    local current_file="$TEST81_E_ROOT/CURRENT_RUN.txt"
    if [[ -n "${TEST81_RUN_TAG:-}" ]]; then
        RUN_TAG="$TEST81_RUN_TAG"
    elif [[ -s "$current_file" ]]; then
        RUN_TAG="$(tr -d '\r\n' < "$current_file")"
    else
        RUN_TAG="test81_publication_$(date +%Y%m%d_%H%M%S)"
    fi
    [[ "$RUN_TAG" =~ ^test81_publication_[0-9]{8}_[0-9]{6}$ ]] || {
        echo "ERROR: invalid TEST81_RUN_TAG: $RUN_TAG" >&2
        exit 1
    }
    RUN_DIR="$TEST81_E_ROOT/runs/$RUN_TAG"
    SHARED_DIR="$TEST81_E_ROOT/shared"
    BIN_DIR="$RUN_DIR/bin"
    INPUT_DIR="$RUN_DIR/inputs"
    LOG_DIR="$RUN_DIR/logs"
    TMP_DIR="$RUN_DIR/tmp"
    RESULTS_DIR="$RUN_DIR/results"
    mkdir -p "$BIN_DIR" "$INPUT_DIR" "$LOG_DIR" "$TMP_DIR" "$RESULTS_DIR"
    printf '%s\n' "$RUN_TAG" > "$current_file"
    export TMPDIR="$TMP_DIR"
    echo -e "RUN_TAG\t$RUN_TAG"
    printf 'WINDOWS_OUTPUT\t%s\n' "E:\\OligoScout_Benchmarks\\test81_publication\\runs\\$RUN_TAG"
    echo -e "WSL_OUTPUT\t$RUN_DIR"
    echo -e "RESUME\tTEST81_RUN_TAG=$RUN_TAG bash '$TEST81_ARTIFACT_DIR/setup_e_drive_and_pilot.sh'"
}

test81_record_environment() {
    if [[ -f "$RUN_DIR/.environment_complete" ]]; then return; fi
    banner "RECORD IMMUTABLE ENVIRONMENT"
    {
        echo -e "run_tag\t$RUN_TAG"
        echo -e "date_utc\t$(date -u +%FT%TZ)"
        echo -e "hostname\t$(hostname)"
        echo -e "kernel\t$(uname -srmo)"
        echo -e "processor\t$(lscpu | awk -F: '/Model name/ {sub(/^[ \t]+/,"",$2); print $2; exit}')"
        echo -e "logical_cpus\t$(nproc)"
        echo -e "memory_kib\t$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
        echo -e "output_filesystem\t$(findmnt -T /mnt/e -no FSTYPE,OPTIONS)"
        echo -e "reference\t$REFERENCE"
        echo -e "reference_sha256\t$(sha256sum "$REFERENCE" | awk '{print $1}')"
        echo -e "manifest_sha256\t$(sha256sum "$MANIFEST" | awk '{print $1}')"
        echo -e "git_commit\t$(git -C "$TEST81_REPO" rev-parse HEAD)"
        echo -e "git_status_porcelain_count\t$(git -C "$TEST81_REPO" status --porcelain | wc -l)"
        echo -e "compiler\t$(g++ --version | head -1)"
        echo -e "python\t$(python3 --version 2>&1)"
        echo -e "primer3\t$($PRIMER3 --about 2>/dev/null | head -1)"
        echo -e "jellyfish\t$($JELLYFISH --version 2>&1 | head -1)"
        echo -e "blastn\t$($BLASTN -version 2>&1 | head -1)"
        echo -e "bowtie\t$($BOWTIE --version 2>&1 | head -1)"
        echo -e "bwa\t$($BWA 2>&1 | awk '/Version:/ {print $2; exit}')"
        echo -e "search_threads\t1"
        echo -e "panel_generation_threads\t4"
        echo -e "anchor_length\t12"
        echo -e "max_mismatches_per_primer\t3"
        echo -e "amplicon_range_bp\t50..3000"
        echo -e "general_aligner_role\tpost-filtered baseline, not primer-native software"
        echo -e "timing_caveat\toutputs and temporary sort files are on WSL DrvFS mounted E:"
    } > "$RUN_DIR/environment.tsv"
    {
        cd "$TEST81_REPO"
        find include src -type f -print0 | sort -z | xargs -0 sha256sum
        sha256sum CMakeLists.txt
    } > "$RUN_DIR/core_source_before.sha256"
    cp -- "$TEST81_ARTIFACT_DIR/test81_pipeline.py" "$RUN_DIR/"
    cp -- "$TEST81_ARTIFACT_DIR/test81_generate_panel.py" "$RUN_DIR/"
    cp -- "$TEST81_ARTIFACT_DIR/test81_report.py" "$RUN_DIR/"
    cp -- "$TEST81_ARTIFACT_DIR/test81_common.sh" "$RUN_DIR/"
    cp -- "$TEST81_ARTIFACT_DIR/setup_e_drive_and_pilot.sh" "$RUN_DIR/"
    cp -- "$TEST81_ARTIFACT_DIR/launch_full_1024.sh" "$RUN_DIR/"
    cp -- "$TEST81_ARTIFACT_DIR/status_test81.sh" "$RUN_DIR/"
    cp -- "$TEST81_ARTIFACT_DIR/src/test81_anchor_oracle.cpp" "$RUN_DIR/"
    cp -- "$TEST81_ARTIFACT_DIR/src/benchmark_test81_panel.cpp" "$RUN_DIR/"
    touch "$RUN_DIR/.environment_complete"
    echo -e "ENVIRONMENT_RECORD\tCOMPLETE"
}

test81_build_helpers() {
    if [[ -f "$RUN_DIR/.helpers_complete" ]]; then return; fi
    banner "BUILD TEST #81 HELPERS ON E:"
    python3 -m py_compile "$PIPELINE" "$PANEL_GENERATOR"
    g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic \
        -I"$TEST81_REPO/include" \
        "$TEST81_ARTIFACT_DIR/src/test81_anchor_oracle.cpp" \
        -o "$BIN_DIR/test81_anchor_oracle"
    g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic \
        -I"$TEST81_REPO/include" \
        "$TEST81_ARTIFACT_DIR/src/benchmark_test81_panel.cpp" \
        "$CORE_LIBRARY" -pthread -o "$BIN_DIR/benchmark_test81_panel"
    sha256sum "$BIN_DIR/test81_anchor_oracle" \
        "$BIN_DIR/benchmark_test81_panel" > "$BIN_DIR/SHA256SUMS"
    touch "$RUN_DIR/.helpers_complete"
    echo -e "HELPERS\tCOMPILED"
}

test81_build_kmer_index() {
    KMER_DB="$SHARED_DIR/grch38_primary24.canonical.k12.jf"
    if [[ -f "$SHARED_DIR/.k12_complete" && -s "$KMER_DB" ]]; then
        echo -e "PRIMARY24_K12_INDEX\tREUSED"
        return
    fi
    banner "BUILD PRIMARY24 CANONICAL 12-MER INDEX ON E:"
    local temporary="$KMER_DB.partial.$$"
    /usr/bin/time -v -o "$SHARED_DIR/k12_index.time.txt" \
        "$JELLYFISH" count -m 12 -s 20M -t 4 -C \
        -o "$temporary" "$REFERENCE" \
        > "$SHARED_DIR/k12_index.stdout.txt" \
        2> "$SHARED_DIR/k12_index.stderr.txt"
    mv -- "$temporary" "$KMER_DB"
    sha256sum "$KMER_DB" > "$SHARED_DIR/k12_index.sha256"
    touch "$SHARED_DIR/.k12_complete"
    echo -e "PRIMARY24_K12_INDEX\tCOMPLETE"
}

test81_generate_panel() {
    PANEL1024="$INPUT_DIR/test81_panel_1024.tsv"
    PANEL128="$INPUT_DIR/test81_panel_128.tsv"
    if [[ ! -f "$RUN_DIR/.panel_complete" ]]; then
        banner "DESIGN 1,024 PRIMER PAIRS WITH FOUR DIFFICULTY STRATA"
        if [[ "$(($(wc -l < "$PANEL1024" 2>/dev/null || echo 0) - 1))" -ne 1024 ]] || \
           [[ ! -s "$INPUT_DIR/test81_panel_1024.metadata.json" ]]; then
            /usr/bin/time -v -o "$LOG_DIR/panel_generation.time.txt" \
                python3 "$PANEL_GENERATOR" \
                    --reference "$REFERENCE" --fai "$REFERENCE_FAI" \
                    --primer3 "$PRIMER3" --jellyfish "$JELLYFISH" \
                    --kmer-db "$SHARED_DIR/grch38_primary24.canonical.k12.jf" \
                    --output "$PANEL1024" \
                    --metadata "$INPUT_DIR/test81_panel_1024.metadata.json" \
                    --work-dir "$TMP_DIR/panel_generation" \
                    --pairs 1024 --candidate-windows 6000 --seed 810024 \
                > "$LOG_DIR/panel_generation.stdout.txt" \
                2> "$LOG_DIR/panel_generation.stderr.txt"
        else
            echo -e "PANEL_GENERATION\tRECOVERED_COMPLETED_OUTPUT"
        fi
        head -n 129 "$PANEL1024" > "$PANEL128"
        [[ "$(($(wc -l < "$PANEL1024") - 1))" -eq 1024 ]]
        [[ "$(($(wc -l < "$PANEL128") - 1))" -eq 128 ]]
        sha256sum "$PANEL1024" "$PANEL128" > "$INPUT_DIR/panel_SHA256SUMS"
        touch "$RUN_DIR/.panel_complete"
    fi
    echo -e "PANEL_1024\t$PANEL1024"
    echo -e "PANEL_128\t$PANEL128"
}

sort_bindings() {
    local input="$1" output="$2" temporary="${2}.partial"
    { head -n 1 "$input"; tail -n +2 "$input" | LC_ALL=C sort \
        --parallel=2 -S 1536M -T "$TMP_DIR" -t $'\t' \
        -k1,1n -k4,4 -k2,2 -k5,5 -k6,6n -k7,7n -k9,9n -k10,10 -u; \
    } > "$temporary"
    mv -- "$temporary" "$output"
}

sort_products() {
    local input="$1" output="$2" temporary="${2}.partial"
    { head -n 1 "$input"; tail -n +2 "$input" | LC_ALL=C sort \
        --parallel=2 -S 1536M -T "$TMP_DIR" -t $'\t' \
        -k1,1n -k2,2 -k7,7n -k8,8n -k3,3 -k4,4 \
        -k10,10n -k11,11n -k17,17 -k18,18 -u; \
    } > "$temporary"
    mv -- "$temporary" "$output"
}

prepare_queries() {
    local panel="$1" directory="$2"
    if [[ -f "$directory/.queries_complete" ]]; then return; fi
    mkdir -p "$directory/queries" "$directory/times"
    python3 "$PIPELINE" make-queries --panel "$panel" \
        --fasta "$directory/queries/primers.fa" \
        --fastq "$directory/queries/primers.fq" \
        > "$directory/queries/make_queries.log"
    touch "$directory/.queries_complete"
}

finalize_external_method() {
    local method="$1" panel="$2" directory="$3" normalized="$4"
    if [[ ! -s "$directory/bindings.tsv" ]]; then
        /usr/bin/time -v -o "$directory/times/sort_bindings.time.txt" \
            bash -c 'source "$1"; TMP_DIR="$2"; sort_bindings "$3" "$4"' \
            _ "$TEST81_ARTIFACT_DIR/test81_common.sh" "$TMP_DIR" \
            "$normalized" "$directory/bindings.tsv"
    fi
    if [[ ! -s "$directory/products.unsorted.tsv" && ! -s "$directory/products.tsv" ]]; then
        /usr/bin/time -v -o "$directory/times/assemble.time.txt" \
            python3 "$PIPELINE" assemble --bindings "$directory/bindings.tsv" \
                --output "$directory/products.unsorted.tsv" \
                --min-amplicon 50 --max-amplicon 3000 \
            > "$directory/assemble.log" 2>&1
    fi
    if [[ ! -s "$directory/products.tsv" ]]; then
        /usr/bin/time -v -o "$directory/times/sort_products.time.txt" \
            bash -c 'source "$1"; TMP_DIR="$2"; sort_products "$3" "$4"' \
            _ "$TEST81_ARTIFACT_DIR/test81_common.sh" "$TMP_DIR" \
            "$directory/products.unsorted.tsv" "$directory/products.tsv"
    fi
    rm -f -- "$normalized" "$directory/products.unsorted.tsv"
    sha256sum "$directory/bindings.tsv" "$directory/products.tsv" \
        > "$directory/SHA256SUMS"
    touch "$directory/.complete"
    echo -e "METHOD_COMPLETE\t$method"
}

run_oligoscout() {
    local panel="$1" directory="$2"
    [[ -f "$directory/.complete" ]] && return
    banner "OLIGOSCOUT — DIRECT ANCHOR-AWARE SEARCH"
    mkdir -p "$directory/times"
    if [[ ! -f "$directory/.search_complete" ]]; then
        /usr/bin/time -v -o "$directory/times/search.time.txt" \
            "$BIN_DIR/benchmark_test81_panel" \
                "$MANIFEST" "$PPFM_DIR" "$panel" \
                "$directory/bindings.unsorted.tsv" \
                "$directory/products.unsorted.tsv" 12 3 50 3000 \
            > "$directory/search.stdout.txt" 2> "$directory/search.stderr.txt"
        touch "$directory/.search_complete"
    fi
    if [[ ! -s "$directory/bindings.tsv" ]]; then
        /usr/bin/time -v -o "$directory/times/sort_bindings.time.txt" \
            bash -c 'source "$1"; TMP_DIR="$2"; sort_bindings "$3" "$4"' \
            _ "$TEST81_ARTIFACT_DIR/test81_common.sh" "$TMP_DIR" \
            "$directory/bindings.unsorted.tsv" "$directory/bindings.tsv"
    fi
    if [[ ! -s "$directory/products.tsv" ]]; then
        /usr/bin/time -v -o "$directory/times/sort_products.time.txt" \
            bash -c 'source "$1"; TMP_DIR="$2"; sort_products "$3" "$4"' \
            _ "$TEST81_ARTIFACT_DIR/test81_common.sh" "$TMP_DIR" \
            "$directory/products.unsorted.tsv" "$directory/products.tsv"
    fi
    rm -f -- "$directory/bindings.unsorted.tsv" "$directory/products.unsorted.tsv"
    sha256sum "$directory/bindings.tsv" "$directory/products.tsv" > "$directory/SHA256SUMS"
    touch "$directory/.complete"
}

run_oracle() {
    local panel="$1" directory="$2"
    [[ -f "$directory/.complete" ]] && return
    banner "INDEPENDENT FULL-REFERENCE ORACLE"
    mkdir -p "$directory/times"
    if [[ ! -f "$directory/.search_complete" ]]; then
        /usr/bin/time -v -o "$directory/times/search.time.txt" \
            "$BIN_DIR/test81_anchor_oracle" "$REFERENCE" "$panel" \
                "$directory/bindings.unsorted.tsv" 12 3 \
            > "$directory/search.stdout.txt" 2> "$directory/search.stderr.txt"
        touch "$directory/.search_complete"
    fi
    if [[ ! -s "$directory/bindings.tsv" ]]; then
        /usr/bin/time -v -o "$directory/times/sort_bindings.time.txt" \
            bash -c 'source "$1"; TMP_DIR="$2"; sort_bindings "$3" "$4"' \
            _ "$TEST81_ARTIFACT_DIR/test81_common.sh" "$TMP_DIR" \
            "$directory/bindings.unsorted.tsv" "$directory/bindings.tsv"
    fi
    if [[ ! -s "$directory/products.unsorted.tsv" && ! -s "$directory/products.tsv" ]]; then
        /usr/bin/time -v -o "$directory/times/assemble.time.txt" \
            python3 "$PIPELINE" assemble --bindings "$directory/bindings.tsv" \
                --output "$directory/products.unsorted.tsv" \
            > "$directory/assemble.log" 2>&1
    fi
    if [[ ! -s "$directory/products.tsv" ]]; then
        /usr/bin/time -v -o "$directory/times/sort_products.time.txt" \
            bash -c 'source "$1"; TMP_DIR="$2"; sort_products "$3" "$4"' \
            _ "$TEST81_ARTIFACT_DIR/test81_common.sh" "$TMP_DIR" \
            "$directory/products.unsorted.tsv" "$directory/products.tsv"
    fi
    rm -f -- "$directory/bindings.unsorted.tsv" "$directory/products.unsorted.tsv"
    sha256sum "$directory/bindings.tsv" "$directory/products.tsv" > "$directory/SHA256SUMS"
    touch "$directory/.complete"
}

run_bowtie1() {
    local panel="$1" directory="$2"
    [[ -f "$directory/.complete" ]] && return
    banner "BOWTIE 1 — EXHAUSTIVE UNGAPPED BASELINE"
    prepare_queries "$panel" "$directory"
    if [[ ! -f "$directory/.search_complete" ]]; then
        set -o pipefail
        /usr/bin/time -v -o "$directory/times/search.time.txt" \
            "$BOWTIE" -S -f -v 3 -a --no-unal --seed 0 -p 1 \
                "$BOWTIE_INDEX" "$directory/queries/primers.fa" \
            2> "$directory/search.stderr.txt" | gzip -1 > "$directory/raw.sam.gz"
        touch "$directory/.search_complete"
    fi
    if [[ ! -f "$directory/.normalize_complete" ]]; then
        /usr/bin/time -v -o "$directory/times/normalize.time.txt" \
            python3 "$PIPELINE" normalize-sam --panel "$panel" \
                --reference "$REFERENCE" --fai "$REFERENCE_FAI" \
                --input "$directory/raw.sam.gz" \
                --output "$directory/bindings.normalized.tsv" \
            > "$directory/normalize.log" 2>&1
        touch "$directory/.normalize_complete"
    fi
    finalize_external_method bowtie1 "$panel" "$directory" \
        "$directory/bindings.normalized.tsv"
}

run_bwa() {
    local panel="$1" directory="$2"
    [[ -f "$directory/.complete" ]] && return
    banner "BWA-ALN — SEED-DISABLED ALL-HIT BASELINE"
    prepare_queries "$panel" "$directory"
    if [[ ! -f "$directory/.search_complete" ]]; then
        /usr/bin/time -v -o "$directory/times/search.time.txt" \
            bash -c 'set -euo pipefail
            bwa="$1"; index="$2"; fastq="$3"; sai="$4"; raw="$5"; err="$6"
            "$bwa" aln -n 3 -o 0 -e -1 -l 25 -M 1000 -N -t 1 \
                "$index" "$fastq" > "$sai" 2> "$err"
            "$bwa" samse -n 1000000 "$index" "$sai" "$fastq" \
                2>> "$err" | gzip -1 > "$raw"' \
        _ "$BWA" "$BWA_INDEX" "$directory/queries/primers.fq" \
        "$directory/search.sai" "$directory/raw.sam.gz" \
        "$directory/search.stderr.txt"
        rm -f -- "$directory/search.sai"
        touch "$directory/.search_complete"
    fi
    if [[ ! -f "$directory/.normalize_complete" ]]; then
        /usr/bin/time -v -o "$directory/times/normalize.time.txt" \
            python3 "$PIPELINE" normalize-sam --panel "$panel" \
                --reference "$REFERENCE" --fai "$REFERENCE_FAI" \
                --input "$directory/raw.sam.gz" \
                --output "$directory/bindings.normalized.tsv" --include-xa \
            > "$directory/normalize.log" 2>&1
        touch "$directory/.normalize_complete"
    fi
    finalize_external_method bwa_aln "$panel" "$directory" \
        "$directory/bindings.normalized.tsv"
}

run_blast() {
    local panel="$1" directory="$2"
    [[ -f "$directory/.complete" ]] && return
    banner "BLASTN-SHORT — E-VALUE 1100 CONTROL BASELINE"
    prepare_queries "$panel" "$directory"
    mkdir -p "$directory/query_chunks" "$directory/raw_chunks" "$directory/times/blast_chunks"
    if [[ ! -f "$directory/.blast_search_complete" ]]; then
        if ! compgen -G "$directory/query_chunks/queries_*.fa" > /dev/null; then
            python3 "$PIPELINE" split-fasta \
                --input "$directory/queries/primers.fa" \
                --output-dir "$directory/query_chunks" --records-per-chunk 64 \
                > "$directory/query_chunks/split.log"
        fi
        local expected_queries chunked_queries
        expected_queries="$(grep -c '^>' "$directory/queries/primers.fa")"
        chunked_queries="$(grep -h '^>' "$directory"/query_chunks/queries_*.fa | wc -l)"
        [[ "$chunked_queries" -eq "$expected_queries" ]] || {
            echo "ERROR: BLAST query chunks contain $chunked_queries queries; expected $expected_queries" >&2
            exit 1
        }
        local chunk stem
        for chunk in "$directory"/query_chunks/queries_*.fa; do
            stem="$(basename "$chunk" .fa)"
            [[ -s "$directory/raw_chunks/$stem.tsv.gz" ]] && continue
            set -o pipefail
            /usr/bin/time -v -o "$directory/times/blast_chunks/$stem.time.txt" \
                "$BLASTN" -task blastn-short -query "$chunk" -db "$BLAST_DB" \
                    -strand both -word_size 7 -dust no -soft_masking false \
                    -evalue 1100 -num_threads 1 -max_target_seqs 24 \
                    -outfmt '6 qseqid sseqid qstart qend sstart send sstrand length mismatch gaps qseq sseq evalue bitscore score' \
                2> "$directory/raw_chunks/$stem.stderr.txt" | \
                gzip -1 > "$directory/raw_chunks/$stem.tsv.gz"
        done
        find "$directory/raw_chunks" -maxdepth 1 -name 'queries_*.tsv.gz' \
            -print0 | sort -z | xargs -0 cat > "$directory/raw.tsv.gz"
        touch "$directory/.blast_search_complete"
    fi
    if [[ ! -f "$directory/.normalize_complete" ]]; then
        /usr/bin/time -v -o "$directory/times/normalize.time.txt" \
            python3 "$PIPELINE" normalize-blast --panel "$panel" \
                --reference "$REFERENCE" --fai "$REFERENCE_FAI" \
                --input "$directory/raw.tsv.gz" \
                --output "$directory/bindings.normalized.tsv" \
            > "$directory/normalize.log" 2>&1
        touch "$directory/.normalize_complete"
    fi
    finalize_external_method blastn_short_e1100 "$panel" "$directory" \
        "$directory/bindings.normalized.tsv"
}

compare_method() {
    local method="$1" panel="$2" stage_dir="$3"
    local oracle="$stage_dir/oracle" directory="$stage_dir/$method/comparison"
    mkdir -p "$directory"
    for kind in bindings products; do
        python3 "$PIPELINE" compare --panel "$panel" \
            --truth "$oracle/$kind.tsv" \
            --prediction "$stage_dir/$method/$kind.tsv" --kind "$kind" \
            --summary "$directory/$kind.summary.tsv" \
            --json "$directory/$kind.summary.json" \
            --per-pair "$directory/$kind.per_pair.tsv" \
            > "$directory/$kind.stdout.txt"
    done
}

test81_core_unchanged() {
    {
        cd "$TEST81_REPO"
        find include src -type f -print0 | sort -z | xargs -0 sha256sum
        sha256sum CMakeLists.txt
    } > "$RUN_DIR/core_source_after.sha256"
    cmp -s "$RUN_DIR/core_source_before.sha256" "$RUN_DIR/core_source_after.sha256" || {
        echo "ERROR: core source changed during Test #81" >&2
        exit 1
    }
    echo -e "CORE_SOURCE_CHANGED\tNO"
}
