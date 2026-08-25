#!/usr/bin/env bash
set -euo pipefail

ARTIFACT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_ROOT="${TEST81_OUTPUT_ROOT:-/mnt/e/OligoScout_Benchmarks/test81_publication}"
CURRENT="$OUTPUT_ROOT/CURRENT_RUN.txt"
[[ -s "$CURRENT" ]] || { echo "ERROR: CURRENT_RUN.txt missing" >&2; exit 1; }
RUN_TAG="${TEST81_RUN_TAG:-$(tr -d '\r\n' < "$CURRENT")}"
RUN="$OUTPUT_ROOT/runs/$RUN_TAG"
BLAST="$RUN/results/full1024/blastn_e1100"
ORIGINAL="$BLAST/query_chunks/queries_0031.fa"
RECOVERY="$BLAST/adaptive_chunk_recovery_queries_0031"

if pgrep -af 'launch_full_1024.sh|blastn' | grep -v grep >/dev/null; then
    echo "ERROR: a full benchmark or BLAST process is still running" >&2
    exit 1
fi
[[ -f "$RUN/.pilot_pass" ]] || { echo "ERROR: pilot did not pass" >&2; exit 1; }
[[ ! -f "$RUN/.full_complete" ]] || { echo "FULL RUN ALREADY COMPLETE"; exit 0; }
[[ -s "$ORIGINAL" ]] || { echo "ERROR: original queries_0031.fa missing" >&2; exit 1; }

if [[ -f "$RECOVERY/.complete" ]]; then
    echo -e "ADAPTIVE_CHUNK_RECOVERY\tALREADY_COMPLETE"
    echo -e "RUN_TAG\t$RUN_TAG"
    exit 0
fi

mkdir -p "$RECOVERY/original" "$RECOVERY/split"
cp -- "$ORIGINAL" "$RECOVERY/original/queries_0031.original64.fa"
[[ "$(grep -c '^>' "$RECOVERY/original/queries_0031.original64.fa")" -eq 64 ]] || {
    echo "ERROR: recovery source does not contain 64 queries" >&2
    exit 1
}

python3 "$ARTIFACT_DIR/test81_pipeline.py" split-fasta \
    --input "$RECOVERY/original/queries_0031.original64.fa" \
    --output-dir "$RECOVERY/split" --records-per-chunk 8 \
    > "$RECOVERY/split.log"

for offset in 0 1 2 3 4 5 6 7; do
    source_file="$(printf '%s/split/queries_%04d.fa' "$RECOVERY" "$offset")"
    destination_index="$((31 + offset))"
    destination="$(printf '%s/query_chunks/queries_%04d.fa' "$BLAST" "$destination_index")"
    [[ "$(grep -c '^>' "$source_file")" -eq 8 ]]
    cp -- "$source_file" "$destination"
done

TOTAL_QUERIES="$(grep -h '^>' "$BLAST"/query_chunks/queries_*.fa | wc -l)"
TOTAL_CHUNKS="$(find "$BLAST/query_chunks" -maxdepth 1 -type f -name 'queries_*.fa' | wc -l)"
[[ "$TOTAL_QUERIES" -eq 2048 ]] || {
    echo "ERROR: adaptive chunks contain $TOTAL_QUERIES queries; expected 2048" >&2
    exit 1
}
[[ "$TOTAL_CHUNKS" -eq 39 ]] || {
    echo "ERROR: adaptive chunk count is $TOTAL_CHUNKS; expected 39" >&2
    exit 1
}

{
    echo -e "event\tADAPTIVE_BLAST_QUERY_CHUNKING"
    echo -e "run_tag\t$RUN_TAG"
    echo -e "trigger\t64_QUERY_FINAL_CHUNK_TERMINATED_TWICE"
    echo -e "completed_original_chunks_preserved\t31"
    echo -e "recovered_original_chunk\tqueries_0031"
    echo -e "original_queries\t64"
    echo -e "replacement_chunks\t8"
    echo -e "queries_per_replacement_chunk\t8"
    echo -e "total_queries_after_recovery\t$TOTAL_QUERIES"
    echo -e "total_chunks_after_recovery\t$TOTAL_CHUNKS"
    echo -e "search_contract_changed\tNO"
    echo -e "scientific_query_set_changed\tNO"
    echo -e "resource_execution_schedule_changed\tYES"
} > "$RECOVERY/RECOVERY.tsv"

(
    cd "$RECOVERY"
    sha256sum RECOVERY.tsv original/queries_0031.original64.fa split/queries_*.fa \
        > SHA256SUMS
    sha256sum -c SHA256SUMS
)
touch "$RECOVERY/.complete"

echo -e "ADAPTIVE_CHUNK_RECOVERY\tCOMPLETE"
echo -e "PRESERVED_BLAST_CHUNKS\t31"
echo -e "REPLACEMENT_CHUNKS\t8_X_8_QUERIES"
echo -e "TOTAL_QUERY_COUNT\t$TOTAL_QUERIES"
echo -e "TOTAL_CHUNK_COUNT\t$TOTAL_CHUNKS"
echo -e "SCIENTIFIC_QUERY_SET_CHANGED\tNO"
echo
echo "Resume with:"
echo "TEST81_RUN_TAG=$RUN_TAG bash '$ARTIFACT_DIR/launch_full_1024.sh'"
