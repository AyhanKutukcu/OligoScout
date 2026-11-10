#!/usr/bin/env bash
# =============================================================================
# BIS-04: bisulfit iplik ozgullugu — donusturulmus genom taramasi
#
#   bash setup/run_bisulfite_panel.sh                 depodaki chr22 paneli
#   PANEL=/yol/panel.tsv bash setup/run_bisulfite_panel.sh   kendi paneliniz
#
# Primerleri HEM CT HEM GA donusturulmus genom indeksinde arar (RUN_BISULFITE_
# INDEX.cmd ile kurulan), sonra cli/bisulfite_qc.py ile cift basina BIS-04
# kararini uretir. Yeniden DERLEME YOK: ayni oligoscout_panel ikilisi env ile
# farkli indeks dizinine yonlendirilir.
#
# ON KOSUL: data/indexes/bisulfite/{CT,GA}/manifest.sa8.tsv mevcut olmali.
# Ayarlar: ANCHOR, MAXMM, MINAMP, MAXAMP, INTENDED (CT|GA), SA_RATE
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
REPO="$ROOT/OligoScout"
PRESET="${PRESET:-release}"
BUILDDIR="$REPO/build/$PRESET"
SA_RATE="${SA_RATE:-8}"
BISDIR="$REPO/data/indexes/bisulfite"

DEFAULT_PANEL="$REPO/data/benchmarks/grch38_chr22_panel_32.tsv"
PANEL="${PANEL:-$DEFAULT_PANEL}"
OUT="${OUT:-$ROOT/reference_py/out/bisulfite_panel}"
ANCHOR="${ANCHOR:-12}"; MAXMM="${MAXMM:-3}"
MINAMP="${MINAMP:-50}"; MAXAMP="${MAXAMP:-3000}"
INTENDED="${INTENDED:-CT}"

[[ -x "$BUILDDIR/oligoscout_panel" ]] || {
    echo "oligoscout_panel yok; once: bash setup/run_panel.sh (derler)"; exit 1; }
for conv in CT GA; do
    [[ -f "$BISDIR/$conv/manifest.sa${SA_RATE}.tsv" ]] || {
        echo "Bisulfit indeksi yok: $BISDIR/$conv (once RUN_BISULFITE_INDEX.cmd)";
        exit 1; }
done
[[ -f "$PANEL" ]] || { echo "panel yok: $PANEL"; exit 1; }
mkdir -p "$OUT/CT" "$OUT/GA"

for conv in CT GA; do
    echo
    echo "===== $conv donusturulmus genomda arama ====="
    "$BUILDDIR/oligoscout_panel" \
        "$BISDIR/$conv/manifest.sa${SA_RATE}.tsv" "$BISDIR/$conv" \
        "$PANEL" "$OUT/$conv/panel" \
        --anchor "$ANCHOR" --max-mismatch "$MAXMM" \
        --min-amplicon "$MINAMP" --max-amplicon "$MAXAMP" --cache 1
done

echo
echo "===== BIS-04 iplik ozgullugu birlestirmesi ====="
cd "$ROOT/reference_py"
python3 cli/bisulfite_qc.py \
    --ct-prefix "$OUT/CT/panel" --ga-prefix "$OUT/GA/panel" \
    --out "$OUT" --intended "$INTENDED" --panel "$PANEL"

echo
echo "Tamamlandi. BIS-04: $OUT/bisulfite_qc.tsv"
