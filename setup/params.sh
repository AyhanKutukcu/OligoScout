#!/usr/bin/env bash
# =============================================================================
# Parametre aktarim is akisi
#
#   bash setup/params.sh template   params_in/ icine bos sablonlar
#   bash setup/params.sh check      params_in/*.tsv denetimi
#   bash setup/params.sh emit       params.py'ye yapistirilacak bloklar
#
# Cift aktarim: ikinci kisinin dosyasini params_in/<tablo>.2.tsv olarak
# kaydedin; check ikisini otomatik karsilastirir.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PY_DIR="$ROOT/reference_py"
IN="${IN:-$ROOT/params_in}"
OUT="${OUT:-$ROOT/params_out}"
TABLES=(mismatch dangling loop)

mkdir -p "$IN" "$OUT"
cd "$PY_DIR"

cmd="${1:-template}"

case "$cmd" in
  template)
    for t in "${TABLES[@]}"; do
        if [[ -s "$IN/$t.tsv" ]]; then
            echo "atlandi (zaten var): $IN/$t.tsv"
        else
            python3 cli/param_transcribe.py template "$t" "$IN/$t.tsv"
        fi
    done
    echo
    echo "Sablonlar: $IN"
    echo "Doldurduktan sonra: RUN_PARAM.cmd check"
    ;;

  check)
    status=0
    for t in "${TABLES[@]}"; do
        [[ -s "$IN/$t.tsv" ]] || { echo "--- $t: dosya yok, atlandi"; continue; }
        echo
        echo "===== $t ====="
        if [[ -s "$IN/$t.2.tsv" ]]; then
            python3 cli/param_transcribe.py check "$t" \
                "$IN/$t.tsv" "$IN/$t.2.tsv" || status=1
        else
            python3 cli/param_transcribe.py check "$t" "$IN/$t.tsv" || status=1
        fi
    done
    exit "$status"
    ;;

  emit)
    for t in "${TABLES[@]}"; do
        [[ -s "$IN/$t.tsv" ]] || continue
        python3 cli/param_transcribe.py emit "$t" "$IN/$t.tsv" \
            > "$OUT/$t.py"
        echo "yazildi: $OUT/$t.py"
    done
    echo
    echo "Bu bloklari oligoscout_thermo/params.py icindeki ilgili"
    echo "ParameterTable.data alanina yapistirip verification_status'u"
    echo "'transcribed_unverified' yapin."
    ;;

  *)
    echo "kullanim: bash setup/params.sh [template|check|emit]"
    exit 2
    ;;
esac
