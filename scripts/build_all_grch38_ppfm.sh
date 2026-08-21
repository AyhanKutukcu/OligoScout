#!/usr/bin/env bash

set -euo pipefail


SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &&
    pwd
)"

ROOT="$(
    cd -- "$SCRIPT_DIR/.." &&
    pwd
)"

FASTA="$ROOT/data/genomes/GRCh38.p14.fa"

BUILDER="$ROOT/build/release/build_grch38_ppfm"

OUTDIR="$ROOT/data/indexes/ppfm"

LOGDIR="$ROOT/results/benchmarks/grch38_ppfm_batch"

SA_RATE="${1:-8}"


mkdir -p \
    "$OUTDIR" \
    "$LOGDIR"


if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA not found:"
    echo "$FASTA"
    exit 1
fi


if [[ ! -x "$BUILDER" ]]; then
    echo "ERROR: builder executable not found:"
    echo "$BUILDER"
    exit 1
fi


chromosomes=()

for n in $(seq 1 22); do
    chromosomes+=("chr${n}")
done

chromosomes+=(
    "chrX"
    "chrY"
)


manifest_tmp="$LOGDIR/manifest.sa${SA_RATE}.tsv.tmp"

manifest="$LOGDIR/manifest.sa${SA_RATE}.tsv"


printf \
"chromosome\tfile_bytes\tsha256\tstatus\n" \
> "$manifest_tmp"


echo
echo "========================================"
echo "PrimerPair GRCh38 PPFM batch"
echo "SA rate: $SA_RATE"
echo "========================================"
echo


for chromosome in "${chromosomes[@]}"; do

    final="$OUTDIR/${chromosome}.sa${SA_RATE}.ppfm"

    part="${final}.part"

    log="$LOGDIR/build_${chromosome}_sa${SA_RATE}.tsv"

    time_log="$LOGDIR/build_${chromosome}_sa${SA_RATE}.time.txt"

    tmp_log="${log}.tmp"

    tmp_time="${time_log}.tmp"


    echo
    echo "----------------------------------------"
    echo "$chromosome"
    echo "----------------------------------------"


    if [[ -s "$final" ]]; then

        echo "Existing finalized index found."
        echo "Skipping build."

        status="EXISTING"

    else

        rm -f \
            "$part" \
            "$tmp_log" \
            "$tmp_time"


        echo "Building $chromosome ..."


        if ! /usr/bin/time -v \
            "$BUILDER" \
            "$FASTA" \
            "$chromosome" \
            "$part" \
            "$SA_RATE" \
            > "$tmp_log" \
            2> "$tmp_time"
        then
            echo "ERROR: build failed for $chromosome"

            rm -f "$part"

            exit 1
        fi


        if ! grep -q \
            $'^ALL_CHECKS\tYES$' \
            "$tmp_log"
        then
            echo "ERROR: ALL_CHECKS YES missing for $chromosome"

            rm -f "$part"

            exit 1
        fi


        if [[ ! -s "$part" ]]; then
            echo "ERROR: generated PPFM is empty."

            exit 1
        fi


        mv \
            "$part" \
            "$final"

        mv \
            "$tmp_log" \
            "$log"

        mv \
            "$tmp_time" \
            "$time_log"


        status="BUILT"


        echo "Build completed."

        grep -E \
            'chromosome|sequence_length|sa_rate|build_seconds|save_seconds|ppfm_bytes|ALL_CHECKS' \
            "$log"
    fi


    bytes=$(
        stat -c '%s' \
            "$final"
    )


    sha=$(
        sha256sum \
            "$final" \
        | awk '{print $1}'
    )


    printf \
        "%s\t%s\t%s\t%s\n" \
        "$chromosome" \
        "$bytes" \
        "$sha" \
        "$status" \
        >> "$manifest_tmp"


    echo "bytes=$bytes"
    echo "sha256=$sha"

done


mv \
    "$manifest_tmp" \
    "$manifest"


echo
echo "========================================"
echo "All canonical chromosomes completed."
echo "Manifest:"
echo "$manifest"
echo "========================================"
