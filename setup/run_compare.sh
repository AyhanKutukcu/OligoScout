#!/usr/bin/env bash
# =============================================================================
# Iki profilin karsilastirma sekli (profile_comparison.png + .tsv)
#
#   bash setup/run_compare.sh
#   PANEL_A=out/panel_design_qpcr_sybr PANEL_B=out/panel_design_ctdna \
#       bash setup/run_compare.sh
#
# Girdi, RUN_PANEL.cmd'nin yazdigi iki tarama klasorudur. Varsayilan:
#   A = reference_py/out/panel_design_qpcr_sybr   (qPCR profili)
#   B = reference_py/out/panel_design_ctdna       (ctDNA profili)
#
# DIKKAT: Bu adim matplotlib + numpy gerektirir; referans katmanin geri kalani
# yalniz standart kutuphane kullanir. Eksikse asagida kurulum satiri yazilir.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PY_DIR="$ROOT/reference_py"
cd "$PY_DIR"

PANEL_A="${PANEL_A:-$PY_DIR/out/panel_design_qpcr_sybr}"
PANEL_B="${PANEL_B:-$PY_DIR/out/panel_design_ctdna}"
NAME_A="${NAME_A:-qpcr_sybr}"
NAME_B="${NAME_B:-ctdna}"
LABEL_A="${LABEL_A:-qPCR profile\n(70-150 bp)}"
LABEL_B="${LABEL_B:-ctDNA profile\n(60-100 bp)}"
OUT="${OUT:-$PY_DIR/out}"

for d in "$PANEL_A" "$PANEL_B"; do
    if [[ ! -s "$d/panel_qc.tsv" ]]; then
        echo "Tarama klasoru eksik: $d/panel_qc.tsv"
        echo "Once o profili tarayin:"
        echo "  RUN_PANEL.cmd reference_py\\out\\<tasarim_klasoru>\\design_panel.tsv"
        exit 1
    fi
done

if ! python3 -c "import matplotlib, numpy" 2>/dev/null; then
    echo "matplotlib/numpy bulunamadi — sekil cizilemez."
    echo "WSL icinde bir kez kurun:"
    echo "  sudo apt-get install -y python3-matplotlib python3-numpy"
    echo "  (ya da: python3 -m pip install --user matplotlib numpy)"
    exit 1
fi

python3 cli/compare_profiles.py \
    --panel "$NAME_A=$PANEL_A" --panel "$NAME_B=$PANEL_B" \
    --label "$NAME_A=$LABEL_A" --label "$NAME_B=$LABEL_B" \
    --out "$OUT"
