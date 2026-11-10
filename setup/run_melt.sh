#!/usr/bin/env bash
# =============================================================================
# Amplikon erime egrisi ve QPCR-02 (tek erime piki)
#
#   DESIGN=/yol/out/design_qpcr_sybr TARGET=/yol/hedef_chr22.fa \
#       bash setup/run_melt.sh
#
# Ayarlar: OUT, NA, MG, DNTP, T_MIN, T_MAX, STEP, MIN_SEPARATION,
#          TM_RESOLUTION, SEQUENCES (id<TAB>sequence dosyasi)
#
# Model: Zimm-Bragg tipi isbirlikli helis-kangal. sigma KALIBRE EDILMEMIS bir
# model sabitidir, dongu entropisi henuz yoktur (tablolar pending).
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PY_DIR="$ROOT/reference_py"
cd "$PY_DIR"

DESIGN="${DESIGN:-$PY_DIR/out/design_qpcr_sybr}"
TARGET="${TARGET:-$ROOT/hedef_chr22.fa}"
OUT="${OUT:-$PY_DIR/out/melt_$(basename "$DESIGN")}"

if [[ ! -s "$DESIGN/design_pairs.tsv" ]]; then
    echo "Tasarim ciktisi yok: $DESIGN/design_pairs.tsv"
    echo "Once: RUN_DESIGN.cmd <hedef.fa> <profil>"
    exit 1
fi

ARGS=(--design "$DESIGN" --target "$TARGET" --out "$OUT"
      --na "${NA:-50}" --mg "${MG:-1.5}" --dntp "${DNTP:-0.6}"
      --t-min "${T_MIN:-50}" --t-max "${T_MAX:-98}" --step "${STEP:-0.25}"
      --min-separation "${MIN_SEPARATION:-2.0}"
      --tm-resolution "${TM_RESOLUTION:-1.0}")
if [[ -n "${SEQUENCES:-}" ]]; then ARGS+=(--sequences "$SEQUENCES"); fi

mkdir -p "$OUT"
python3 cli/amplicon_melt.py "${ARGS[@]}"
