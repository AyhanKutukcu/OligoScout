#!/usr/bin/env bash
# Dizi/kompozisyon/cift/amplikon kontrolleri (SEQ, CMP, SS-02, PAIR, AMP)
set -Eeuo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PY_DIR="$ROOT/reference_py"
cd "$PY_DIR"
DESIGN="${DESIGN:-$PY_DIR/out/design_ctdna}"
TARGET="${TARGET:-$ROOT/hedef_chr22.fa}"
PROFILE="${PROFILE:-ctdna}"
OUT="${OUT:-$DESIGN}"
if [[ ! -s "$DESIGN/design_pairs.tsv" ]]; then
    echo "Tasarim ciktisi yok: $DESIGN/design_pairs.tsv"
    echo "Once: RUN_DESIGN.cmd <hedef.fa> <profil>"
    exit 1
fi
ARGS=(--design "$DESIGN" --target "$TARGET" --profile "$PROFILE" --out "$OUT"
      --na "${NA:-50}" --mg "${MG:-1.5}" --dntp "${DNTP:-0.6}"
      --primer-nm "${PRIMER_NM:-200}" --ta "${TA:-60}")
if [[ -n "${PRODUCT_BP:-}" ]]; then ARGS+=(--product-bp "$PRODUCT_BP"); fi
if [[ -n "${GENE:-}" ]]; then ARGS+=(--gene "$GENE"); fi
if [[ -n "${PROBES:-}" ]]; then ARGS+=(--probes "$PROBES"); fi
python3 cli/sequence_qc.py "${ARGS[@]}"
