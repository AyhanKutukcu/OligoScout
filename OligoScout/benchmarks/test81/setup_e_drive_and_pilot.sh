#!/usr/bin/env bash
set -euo pipefail

ARTIFACT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$ARTIFACT_DIR/test81_common.sh"

test81_preflight
test81_initialize_run
test81_record_environment
test81_build_helpers
test81_build_kmer_index
test81_generate_panel

if [[ ! -f "$RUN_DIR/.regression_complete" ]]; then
    banner "OLIGOSCOUT RELEASE REGRESSION"
    /usr/bin/time -v -o "$LOG_DIR/ctest_release.time.txt" \
        ctest --test-dir "$TEST81_REPO/build/release" \
            --output-on-failure -j 1 \
        > "$LOG_DIR/ctest_release_78.txt" 2>&1
    grep -Eq '100% tests passed|78/78' "$LOG_DIR/ctest_release_78.txt" || {
        echo "ERROR: release regression did not report a complete pass" >&2
        exit 1
    }
    touch "$RUN_DIR/.regression_complete"
fi
echo -e "RELEASE_REGRESSION\tPASS"

SMOKE_DIR="$RESULTS_DIR/smoke32"
mkdir -p "$SMOKE_DIR/oligoscout"
if [[ ! -f "$SMOKE_DIR/.complete" ]]; then
    banner "LEGACY 32-PAIR SMOKE CONTROL"
    run_oligoscout "$LEGACY_PANEL" "$SMOKE_DIR/oligoscout"
    grep -q $'BINDING_COUNT\t675059' "$SMOKE_DIR/oligoscout/search.stdout.txt"
    grep -q $'PRODUCT_COUNT\t65467' "$SMOKE_DIR/oligoscout/search.stdout.txt"
    {
        echo -e "legacy_pairs\t32"
        echo -e "expected_bindings\t675059"
        echo -e "observed_bindings\t675059"
        echo -e "expected_products\t65467"
        echo -e "observed_products\t65467"
        echo -e "legacy_smoke_status\tPASS"
    } > "$SMOKE_DIR/summary.tsv"
    touch "$SMOKE_DIR/.complete"
fi
echo -e "LEGACY_32_PAIR_SMOKE\tPASS"

PILOT_DIR="$RESULTS_DIR/pilot128"
mkdir -p "$PILOT_DIR"/{oracle,oligoscout,bowtie1,bwa_aln,blastn_e1100}

run_oracle "$PANEL128" "$PILOT_DIR/oracle"
run_oligoscout "$PANEL128" "$PILOT_DIR/oligoscout"
run_bowtie1 "$PANEL128" "$PILOT_DIR/bowtie1"
run_bwa "$PANEL128" "$PILOT_DIR/bwa_aln"
run_blast "$PANEL128" "$PILOT_DIR/blastn_e1100"

banner "PILOT ACCURACY COMPARISONS"
for method in oracle oligoscout bowtie1 bwa_aln blastn_e1100; do
    compare_method "$method" "$PANEL128" "$PILOT_DIR"
done

mkdir -p "$PILOT_DIR/report"
python3 "$ARTIFACT_DIR/test81_report.py" \
    --stage-dir "$PILOT_DIR" --panel-size 128 \
    --output-dir "$PILOT_DIR/report" \
    > "$PILOT_DIR/report/report_generation.log"

python3 - "$PILOT_DIR" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
for method in ("oligoscout", "bowtie1", "bwa_aln"):
    for kind in ("bindings", "products"):
        path = root / method / "comparison" / f"{kind}.summary.json"
        data = json.loads(path.read_text())
        if not data["exact_set_equality"]:
            raise SystemExit(f"PILOT GATE FAIL: {method} {kind} is not oracle-identical")
for method in ("oligoscout", "bowtie1", "bwa_aln"):
    path = root / method / "comparison" / "products.summary.json"
    data = json.loads(path.read_text())
    if data["intended_products_recovered"] != 128:
        raise SystemExit(
            f"PILOT GATE FAIL: {method} intended recall is "
            f"{data['intended_products_recovered']}/128"
        )
print("PILOT_ORACLE_EQUALITY_GATE\tPASS")
print("PILOT_INTENDED_RECALL_GATE\t128/128_FOR_EXACT_BASELINES")
PY

test81_core_unchanged

{
    echo -e "test\tTEST81_PUBLICATION_PILOT"
    echo -e "run_tag\t$RUN_TAG"
    echo -e "panel_pairs\t128"
    echo -e "primer_sequences\t256"
    echo -e "difficulty_strata\t4_BALANCED"
    echo -e "reference\tGRCh38.p14_PRIMARY24"
    echo -e "anchor_contract\tEXACT_12_NT_3PRIME"
    echo -e "mismatch_contract\tAT_MOST_3_PER_PRIMER"
    echo -e "amplicon_contract\t50_TO_3000_BP"
    echo -e "oracle\tINDEPENDENT_FULL_REFERENCE_SCAN"
    echo -e "oligoscout_exact_vs_oracle\tYES"
    echo -e "bowtie1_exact_vs_oracle\tYES"
    echo -e "bwa_aln_exact_vs_oracle\tYES"
    echo -e "blast_gate\tINTENDED_RECALL_ONLY;FULL_SET_DISCREPANCIES_REPORTED"
    echo -e "intended_recall\tREPORTED_PER_METHOD"
    echo -e "core_source_changed\tNO"
    echo -e "full_1024_launched\tNO"
    echo -e "pilot_status\tPASS"
} > "$PILOT_DIR/PILOT_SUMMARY.tsv"

sha256sum "$PILOT_DIR/PILOT_SUMMARY.tsv" \
    "$PILOT_DIR/report/test81_results.tsv" \
    "$PILOT_DIR/report/test81_report.md" \
    "$PILOT_DIR/report/test81_product_f1.svg" \
    "$PILOT_DIR/report/test81_search_runtime.svg" \
    > "$PILOT_DIR/REPORT_SHA256SUMS"
touch "$RUN_DIR/.pilot_pass"

banner "TEST #81 PILOT COMPLETE"
cat "$PILOT_DIR/PILOT_SUMMARY.tsv"
echo
echo -e "PILOT_REPORT\t$PILOT_DIR/report/test81_report.md"
printf 'WINDOWS_PILOT_REPORT\t%s\n' "E:\\OligoScout_Benchmarks\\test81_publication\\runs\\$RUN_TAG\\results\\pilot128\\report\\test81_report.md"
echo -e "FULL_RUN_STATUS\tNOT_LAUNCHED"
echo
echo "Pilot sonuçlarını gördükten sonra 1.024 çiftlik tam koşu komutu:"
echo "TEST81_RUN_TAG=$RUN_TAG bash '$ARTIFACT_DIR/launch_full_1024.sh'"
echo
echo "TEST81_PILOT_COMPLETE    YES"
