#!/usr/bin/env bash
set -euo pipefail

ROOT="${TEST81_OUTPUT_ROOT:-/mnt/e/OligoScout_Benchmarks/test81_publication}"
CURRENT="$ROOT/CURRENT_RUN.txt"
[[ -s "$CURRENT" ]] || { echo "TEST81_STATUS    NOT_STARTED"; exit 0; }
RUN_TAG="$(tr -d '\r\n' < "$CURRENT")"
RUN="$ROOT/runs/$RUN_TAG"

echo -e "RUN_TAG\t$RUN_TAG"
echo -e "RUN_DIRECTORY\t$RUN"
echo -e "USED_SPACE\t$(du -sh "$RUN" 2>/dev/null | awk '{print $1}')"
echo -e "E_DRIVE_FREE\t$(df -h /mnt/e | awk 'NR==2 {print $4}')"

FULL_DIR="$RUN/results/full1024"
FULL_METHODS_COMPLETE="$(find "$FULL_DIR" -mindepth 2 -maxdepth 2 -type f -name '.complete' 2>/dev/null | wc -l)"
BLAST_CHUNKS_COMPLETE="$(find "$FULL_DIR/blastn_e1100/raw_chunks" -maxdepth 1 -type f -name 'queries_*.tsv.gz' -size +0c 2>/dev/null | wc -l)"
BLAST_CHUNKS_TOTAL="$(find "$FULL_DIR/blastn_e1100/query_chunks" -maxdepth 1 -type f -name 'queries_*.fa' 2>/dev/null | wc -l)"

if [[ -d "$FULL_DIR" ]]; then
    echo -e "FULL_METHODS_COMPLETE\t$FULL_METHODS_COMPLETE/5"
    echo -e "FULL_BLAST_CHUNKS_COMPLETE\t$BLAST_CHUNKS_COMPLETE/$BLAST_CHUNKS_TOTAL"
fi

if [[ -f "$RUN/.full_complete" ]]; then
    echo -e "TEST81_STATUS\tFULL_1024_COMPLETE"
elif pgrep -af 'launch_full_1024.sh' | grep -v grep >/dev/null; then
    echo -e "TEST81_STATUS\tFULL_1024_RUNNING"
elif [[ "$FULL_METHODS_COMPLETE" -gt 0 || "$BLAST_CHUNKS_COMPLETE" -gt 0 ]]; then
    echo -e "TEST81_STATUS\tFULL_1024_PARTIAL_STOPPED;RERUN_FULL_COMMAND_TO_RESUME"
elif pgrep -af 'setup_e_drive_and_pilot.sh' | grep -v grep >/dev/null; then
    echo -e "TEST81_STATUS\tPILOT_128_RUNNING"
elif [[ -f "$RUN/.pilot_pass" ]]; then
    echo -e "TEST81_STATUS\tPILOT_128_PASS_FULL_NOT_STARTED"
else
    echo -e "TEST81_STATUS\tSTOPPED_OR_BETWEEN_STAGES;RERUN_SAME_COMMAND_TO_RESUME"
fi

find "$RUN/results" -type f \( -name search.stderr.txt -o -name normalize.log \) \
    -printf '%T@\t%p\n' 2>/dev/null | sort -n | tail -n 3 | cut -f2- | \
    sed 's/^/RECENT_LOG\t/' || true
