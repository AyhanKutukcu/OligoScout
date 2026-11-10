#!/usr/bin/env bash
# =============================================================================
# RT-03 izoform kapsami — GENCODE ekson anotasyonundan (transkriptom FM-index YOK)
#
#   PAIRS=/yol/pairs.tsv bash setup/run_isoform.sh
#
# pairs.tsv sutunlari (GENOM koordinatlari, 0-tabanli yari-acik):
#   pair_id  gene  chrom  fwd_start  fwd_end  rev_start  rev_end
#
# Ekson indeksi (transcript_exons.tsv) yoksa GENCODE GTF'den uretir
# (data/annotations icinde *.annotation.gtf.gz bulunmali — RUN_ANNOTATIONS ile
# ayni kaynak). Ayarlar: EXONS, OUT, ANNDIR
# =============================================================================
set -Eeuo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
ANNDIR="${ANNDIR:-$ROOT/OligoScout/data/annotations}"
EXONS="${EXONS:-$ANNDIR/transcript_exons.tsv}"
OUT="${OUT:-$ROOT/reference_py/out/isoform}"

[[ -n "${PAIRS:-}" ]] || { echo "PAIRS verilmedi (cift koordinat TSV'si)"; exit 1; }
[[ -f "$PAIRS" ]] || { echo "pairs yok: $PAIRS"; exit 1; }

if [[ ! -s "$EXONS" ]]; then
    echo "transcript_exons.tsv yok; GENCODE'dan uretiliyor ..."
    python3 "$ROOT/reference_py/cli/build_annotations.py" --annotations "$ANNDIR"
fi
[[ -s "$EXONS" ]] || { echo "ekson indeksi uretilemedi (GENCODE GTF yok mu?)"; exit 1; }

python3 "$ROOT/reference_py/cli/isoform_qc.py" \
    --pairs "$PAIRS" --exons "$EXONS" --out "$OUT"
echo
echo "Tamamlandi. RT-03: $OUT/isoform_qc.tsv"
