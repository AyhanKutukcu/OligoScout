#!/usr/bin/env bash
# =============================================================================
# Varyant katmani (VAR-01..VAR-08): gnomAD v4 ile alel dusmesi taramasi
#
#   PANEL=/yol/out/panel_design_ctdna DESIGN=/yol/out/design_ctdna \
#       bash setup/run_variant.sh
#
# Ayarlar: SCOPE (intended|products), AF_COMMON, AF_WARN, THREE_PRIME_NT,
#          OFFLINE=1 (ag'a cikma, onbellegi kullan), CACHE, OUT
#
# Cekilen her bolge onbellege yazilir (varsayilan reference_py/out/
# variant_cache.jsonl); ayni bolge ikinci kez sorgulanmaz.
# Internet yoksa OFFLINE=1 ile onbellekten calisir.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PY_DIR="$ROOT/reference_py"
cd "$PY_DIR"

PANEL="${PANEL:-$PY_DIR/out/panel_design_ctdna}"
DESIGN="${DESIGN:-$PY_DIR/out/design_ctdna}"
SCOPE="${SCOPE:-intended}"
CACHE="${CACHE:-$PY_DIR/out/variant_cache.jsonl}"
OUT="${OUT:-$PY_DIR/out/variant_$(basename "$DESIGN")}"

if [[ ! -s "$PANEL/panel.products.tsv" ]]; then
    echo "Tarama ciktisi yok: $PANEL/panel.products.tsv"
    echo "Once: RUN_PANEL.cmd reference_py\\out\\<tasarim>\\design_panel.tsv"
    exit 1
fi

ARGS=(--panel "$PANEL" --design "$DESIGN" --out "$OUT" --scope "$SCOPE"
      --cache "$CACHE")
if [[ "${OFFLINE:-}" == "1" ]];     then ARGS+=(--offline); fi
if [[ -n "${AF_COMMON:-}" ]];       then ARGS+=(--af-common "$AF_COMMON"); fi
if [[ -n "${AF_WARN:-}" ]];         then ARGS+=(--af-warn "$AF_WARN"); fi
if [[ -n "${THREE_PRIME_NT:-}" ]];  then ARGS+=(--three-prime-nt "$THREE_PRIME_NT"); fi

mkdir -p "$OUT"
python3 cli/variant_qc.py "${ARGS[@]}"

# Varyant karari nihai siralamaya islenir; iki katman ayri tablolarda kalmaz.
if [[ -s "$DESIGN/design_pairs.tsv" && -s "$PANEL/panel_qc.tsv" ]]; then
    printf '\n===== Varyantla yeniden siralama =====\n'
    python3 cli/rerank.py --design "$DESIGN" --qc "$PANEL" --variants "$OUT"
fi
