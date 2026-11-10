#!/usr/bin/env bash
# =============================================================================
# Tekrarlanabilir uctan uca regresyon (test paketi + butun CLI zinciri)
#
#   bash setup/run_regresyon.sh
#
# Ayarlar: OUT, TIMEOUT, LATEST (son panelin selected_panel.tsv yolu), HIZLI=1
#
# Cikis kodu: herhangi bir adim FAIL ise 1.
# Tek rapor: $OUT/REGRESYON_RAPORU.md
# =============================================================================
set -Euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PY_DIR="$ROOT/reference_py"
OUT="${OUT:-$PY_DIR/out/regresyon}"

cd "$PY_DIR"

ARGS=(--out "$OUT" --timeout "${TIMEOUT:-900}")
if [[ -n "${LATEST:-}" ]]; then ARGS+=(--latest "$LATEST"); fi
if [[ "${HIZLI:-}" == "1" ]]; then ARGS+=(--hizli); fi

python3 cli/regresyon.py "${ARGS[@]}"
code=$?

echo
echo "rapor: $OUT/REGRESYON_RAPORU.md"
exit $code
