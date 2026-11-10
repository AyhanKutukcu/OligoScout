#!/usr/bin/env bash
# =============================================================================
# Tek komutluk kontrol: test paketi + parametre capraz denetimi +
# terminal tablo denetimi.
#
#   bash setup/run_kontrol.sh
#
# Ayarlar: OUT, TA, NA, MG, DNTP, PRIMER_NM
#
# Hicbir adim digerini durdurmaz; her adimin cikis kodu rapora yazilir.
# Tek dosya uretir: $OUT/KONTROL_RAPORU.md  -- bu dosyayi geri gonderin.
# =============================================================================
set -Euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PY_DIR="$ROOT/reference_py"
OUT="${OUT:-$PY_DIR/out/kontrol}"
REPORT="$OUT/KONTROL_RAPORU.md"

mkdir -p "$OUT"
cd "$PY_DIR"

: > "$REPORT"
{
    echo "# OligoScout kontrol raporu"
    echo
    echo '```text'
    echo "python   : $(python3 -V 2>&1)"
    echo "calisma  : $PY_DIR"
    echo '```'
} >> "$REPORT"

run_step() {
    local title="$1"; shift
    local log="$OUT/$(echo "$title" | tr -c 'A-Za-z0-9' '_').log"
    echo "== $title"
    set +e
    "$@" > "$log" 2>&1
    local code=$?
    set -e
    {
        echo
        echo "## $title"
        echo
        echo "- komut: \`$*\`"
        echo "- cikis kodu: $code"
        echo
        echo '```text'
        tail -n 60 "$log"
        echo '```'
    } >> "$REPORT"
    echo "   cikis kodu $code  ($log)"
    return 0
}

run_step "1 test paketi" python3 -m unittest discover -s tests
run_step "2 parametre capraz denetimi" python3 cli/param_crosscheck.py
run_step "3 terminal tablo denetimi" python3 cli/terminal_table_audit.py \
    --ta "${TA:-60}" --na "${NA:-50}" --mg "${MG:-1.5}" \
    --dntp "${DNTP:-0.6}" --primer-nm "${PRIMER_NM:-200}"

{
    echo
    echo "## Capraz denetimin sapan satirlari"
    echo
    echo '```text'
    awk -F'\t' 'NR==1 || $NF != "MATCH"' out/param_crosscheck/crosscheck.tsv \
        | cut -f1,3,4,5,6,7,9 || true
    echo '```'
    echo
    echo "## Terminal tablo denetimi (ozet)"
    echo
    cat out/terminal_audit/terminal_audit_summary.md 2>/dev/null || echo "yok"
} >> "$REPORT"

echo
echo "rapor: $REPORT"
