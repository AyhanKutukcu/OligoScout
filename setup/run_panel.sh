#!/usr/bin/env bash
# =============================================================================
# Panel boru hatti: derle -> tum genomda ara -> termodinamik QC
#
#   bash setup/run_panel.sh                      depodaki chr22 32-cift paneli
#   PANEL=/yol/panel.tsv bash setup/run_panel.sh kendi paneliniz
#
# Ayarlar: TA, NA, MG, DNTP, PRIMER_NM, TM_WINDOW, ANCHOR, MAXMM,
#          MINAMP, MAXAMP, CACHE, FRAGMENT_BP (ctDNA icin 166)
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
REPO="${REPO:-$ROOT/OligoScout}"
PRESET="${PRESET:-release}"
BUILDDIR="$REPO/build/$PRESET"
SA_RATE="${SA_RATE:-8}"

DEFAULT_PANEL="$REPO/data/benchmarks/grch38_chr22_panel_32.tsv"
PANEL="${PANEL:-$DEFAULT_PANEL}"
MANIFEST="${MANIFEST:-$REPO/results/benchmarks/grch38_ppfm_batch/manifest.sa${SA_RATE}.tsv}"
INDEXDIR="${INDEXDIR:-$REPO/data/indexes/ppfm}"

# Her panel kendi klasorune yazar; farkli kosular birbirini ezmez.
# design_panel.tsv adi tum tasarimlarda ayni oldugu icin, bu durumda klasor
# adi panelin BULUNDUGU dizinden turetilir (orn. out/design_ctdna ->
# out/panel_design_ctdna). Aksi halde profiller birbirinin sonucunu silerdi.
if [[ "$PANEL" == "$DEFAULT_PANEL" ]]; then
    OUT="${OUT:-$ROOT/reference_py/out/panel}"
elif [[ "$(basename "$PANEL")" == "design_panel.tsv" ]]; then
    OUT="${OUT:-$ROOT/reference_py/out/panel_$(basename "$(dirname "$PANEL")")}"
else
    OUT="${OUT:-$ROOT/reference_py/out/panel_$(basename "${PANEL%.*}")}"
fi

ANCHOR="${ANCHOR:-12}"
MAXMM="${MAXMM:-3}"
MINAMP="${MINAMP:-50}"
MAXAMP="${MAXAMP:-3000}"
CACHE="${CACHE:-1}"

FRAGMENT_BP="${FRAGMENT_BP:-0}"   # cfDNA/ctDNA icin 166 verin
SKIP_SCAN="${SKIP_SCAN:-0}"       # 1 ise derleme+tarama atlanir (QC'yi tekrarlar)
TA="${TA:-60}"; NA="${NA:-50}"; MG="${MG:-1.5}"; DNTP="${DNTP:-0.6}"
PRIMER_NM="${PRIMER_NM:-200}"; TM_WINDOW="${TM_WINDOW:-60:66}"

c_cyn=$'\033[36m'; c_grn=$'\033[32m'; c_off=$'\033[0m'
head1(){ printf '\n%s===== %s =====%s\n' "$c_cyn" "$*" "$c_off"; }

[[ -f "$MANIFEST" ]] || { echo "manifest yok: $MANIFEST (once RUN_SETUP.cmd)"; exit 1; }
[[ -f "$PANEL" ]]    || { echo "panel yok: $PANEL"; exit 1; }
mkdir -p "$OUT"

if [[ "$SKIP_SCAN" == "1" ]]; then
    # Tarama cikti dosyalari zaten varsa yeniden uretmeye gerek yok; yalniz
    # QC/multipleks asamalarini yeni ayarlarla tekrarlamak icin.
    if [[ ! -s "$OUT/panel.bindings.tsv" ]]; then
        echo "SKIP_SCAN=1 verildi ama $OUT/panel.bindings.tsv yok."
        echo "Once taramayi bir kez calistirin."
        exit 1
    fi
    head1 "1-2 — Tarama atlandi (mevcut kayitlar kullaniliyor)"
    echo "kaynak: $OUT/panel.bindings.tsv"
else
head1 "1 — oligoscout_panel derleniyor"
if [[ -f "$SCRIPT_DIR/cmake_path.env" ]]; then
    source "$SCRIPT_DIR/cmake_path.env"
fi
( cd "$REPO" && cmake --preset "$PRESET" >/dev/null )
( cd "$REPO" && cmake --build --preset "$PRESET" \
      --target oligoscout_panel -j"$(nproc)" )

head1 "2 — Tum genomda panel aramasi"
"$BUILDDIR/oligoscout_panel" \
    "$MANIFEST" "$INDEXDIR" "$PANEL" "$OUT/panel" \
    --anchor "$ANCHOR" --max-mismatch "$MAXMM" \
    --min-amplicon "$MINAMP" --max-amplicon "$MAXAMP" \
    --cache "$CACHE"
fi

head1 "3 — Termodinamik QC birlestirmesi"
cd "$ROOT/reference_py"
python3 cli/panel_qc.py --prefix "$OUT/panel" --out "$OUT" --panel "$PANEL" \
    --ta "$TA" --na "$NA" --mg "$MG" --dntp "$DNTP" \
    --primer-nm "$PRIMER_NM" --tm-window "$TM_WINDOW" --anchor "$ANCHOR"

head1 "4 — Multipleks capraz reaktivite (MPX-01..MPX-03)"
# Yeni genom taramasi YOK: capraz urunler 2. fazin baglanma kayitlarindan
# kurulur, boylece aramanin eksiksizlik garantisi korunur.
python3 cli/multiplex_qc.py --prefix "$OUT/panel" --out "$OUT/multiplex" \
    --panel "$PANEL" --min-amplicon "$MINAMP" --max-amplicon "$MAXAMP" \
    --na "$NA" --mg "$MG" --dntp "$DNTP" --primer-nm "$PRIMER_NM" \
    --fragment-bp "$FRAGMENT_BP"

# Panel bir tasarim ciktisiysa (design_panel.tsv), siralamayi ozgullukle
# yeniden duzenle. Zincirin son halkasi; ek ayar gerektirmez.
DESIGN_DIR="$(dirname "$PANEL")"
if [[ "$(basename "$PANEL")" == "design_panel.tsv" \
      && -s "$DESIGN_DIR/design_pairs.tsv" ]]; then
    head1 "5 — Ozgullukle yeniden siralama"
    # Varyant taramasi daha once yapildiysa (RUN_VARIANT.cmd), onun karari da
    # nihai siralamaya girer. Yoksa adim eskisi gibi yalniz ozgullukle calisir.
    VAR_DIR="$ROOT/reference_py/out/variant_$(basename "$DESIGN_DIR")"
    RERANK_ARGS=(--design "$DESIGN_DIR" --qc "$OUT")
    if [[ -s "$VAR_DIR/variant_qc.tsv" ]]; then
        RERANK_ARGS+=(--variants "$VAR_DIR")
        echo "varyant katmani   : $VAR_DIR"
    else
        echo "varyant katmani   : yok (RUN_VARIANT.cmd ile eklenebilir)"
    fi
    python3 cli/rerank.py "${RERANK_ARGS[@]}"
fi

printf '\n%sTamamlandi. Ciktilar: %s%s\n' "$c_grn" "$OUT" "$c_off"
