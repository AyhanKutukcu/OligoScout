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

[[ -f "$RUN_DIR/.pilot_pass" ]] || {
    echo "ERROR: 128-pair pilot has not passed for $RUN_TAG" >&2
    echo "Run this first:" >&2
    echo "TEST81_RUN_TAG=$RUN_TAG bash '$ARTIFACT_DIR/setup_e_drive_and_pilot.sh'" >&2
    exit 1
}

FULL_DIR="$RESULTS_DIR/full1024"
mkdir -p "$FULL_DIR"/{oracle,oligoscout,bowtie1,bwa_aln,blastn_e1100}

run_oracle "$PANEL1024" "$FULL_DIR/oracle"
run_oligoscout "$PANEL1024" "$FULL_DIR/oligoscout"
run_bowtie1 "$PANEL1024" "$FULL_DIR/bowtie1"
run_bwa "$PANEL1024" "$FULL_DIR/bwa_aln"
run_blast "$PANEL1024" "$FULL_DIR/blastn_e1100"

banner "FULL 1,024-PAIR ACCURACY COMPARISONS"
for method in oracle oligoscout bowtie1 bwa_aln blastn_e1100; do
    compare_method "$method" "$PANEL1024" "$FULL_DIR"
done

mkdir -p "$FULL_DIR/report"
python3 "$ARTIFACT_DIR/test81_report.py" \
    --stage-dir "$FULL_DIR" --panel-size 1024 \
    --output-dir "$FULL_DIR/report" \
    > "$FULL_DIR/report/report_generation.log"

python3 - "$FULL_DIR" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
oligo_bindings = json.loads(
    (root / "oligoscout/comparison/bindings.summary.json").read_text()
)
oligo_products = json.loads(
    (root / "oligoscout/comparison/products.summary.json").read_text()
)
if not oligo_bindings["exact_set_equality"] or not oligo_products["exact_set_equality"]:
    raise SystemExit("FULL VALIDATION FAIL: OligoScout differs from independent oracle")
print("FULL_OLIGOSCOUT_ORACLE_EQUALITY\tPASS")
print("FULL_INTENDED_RECALL\tREPORTED_PER_METHOD")
PY

test81_core_unchanged

{
    echo -e "test\tTEST81_PUBLICATION_FULL"
    echo -e "run_tag\t$RUN_TAG"
    echo -e "panel_pairs\t1024"
    echo -e "primer_sequences\t2048"
    echo -e "difficulty_strata\t4_X_256_BALANCED"
    echo -e "nested_balanced_subsets\t16_64_128_256_1024_PAIRS"
    echo -e "reference\tGRCh38.p14_PRIMARY24"
    echo -e "anchor_contract\tEXACT_12_NT_3PRIME"
    echo -e "mismatch_contract\tAT_MOST_3_PER_PRIMER"
    echo -e "amplicon_contract\t50_TO_3000_BP"
    echo -e "oracle\tINDEPENDENT_FULL_REFERENCE_SCAN"
    echo -e "oligoscout_exact_vs_oracle\tYES"
    echo -e "intended_recall\tREPORTED_PER_METHOD"
    echo -e "bootstrap_unit\tPRIMER_PAIR"
    echo -e "bootstrap_replicates\t10000"
    echo -e "core_source_changed\tNO"
    echo -e "timing_scope\tTHIS_MACHINE_WITH_OUTPUTS_ON_E_DRIVE_DRVFS"
    echo -e "full_status\tCOMPLETE"
} > "$FULL_DIR/FINAL_SUMMARY.tsv"

sha256sum "$FULL_DIR/FINAL_SUMMARY.tsv" \
    "$FULL_DIR/report/test81_results.tsv" \
    "$FULL_DIR/report/test81_report.md" \
    "$FULL_DIR/report/test81_product_f1.svg" \
    "$FULL_DIR/report/test81_search_runtime.svg" \
    > "$FULL_DIR/REPORT_SHA256SUMS"
touch "$RUN_DIR/.full_complete"

banner "TEST #81 FULL BENCHMARK COMPLETE"
cat "$FULL_DIR/FINAL_SUMMARY.tsv"
echo
echo -e "RESULT_TABLE\t$FULL_DIR/report/test81_results.tsv"
echo -e "REPORT\t$FULL_DIR/report/test81_report.md"
echo -e "ACCURACY_FIGURE\t$FULL_DIR/report/test81_product_f1.svg"
echo -e "RUNTIME_FIGURE\t$FULL_DIR/report/test81_search_runtime.svg"
printf 'WINDOWS_RESULT_DIRECTORY\t%s\n' "E:\\OligoScout_Benchmarks\\test81_publication\\runs\\$RUN_TAG\\results\\full1024"
echo
echo "TEST81_FULL_COMPLETE    YES"
