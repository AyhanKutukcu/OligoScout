#!/usr/bin/env bash
# =============================================================================
# Modul 1 — NN termodinamik motoru (Python referans uygulamasi)
#
# Tek script, sirayla:
#   1) parametre kayit defteri (tablolar, durumlari, SHA-256)
#   2) dogrulama paketi (analitik ozdeslikler, degismezlikler, oracle)
#   3) demo panel uzerinde QC raporu + erime egrileri
#   4) depodaki chr22 32-cift paneli uzerinde ayni rapor (varsa)
#
# Kullanim:
#   bash reference_py/RUN_THERMO.sh
#   TA=58 MG=3.0 PRIMER_NM=400 bash reference_py/RUN_THERMO.sh
# =============================================================================
set -Eeuo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$HERE/.." && pwd)"
PANEL="${PANEL:-$ROOT/OligoScout/data/benchmarks/grch38_chr22_panel_32.tsv}"
OUT="${OUT:-$HERE/out}"

TA="${TA:-60}"
NA="${NA:-50}"
MG="${MG:-1.5}"
DNTP="${DNTP:-0.6}"
PRIMER_NM="${PRIMER_NM:-200}"
TM_WINDOW="${TM_WINDOW:-58:62}"

PY="${PY:-python3}"
command -v "$PY" >/dev/null 2>&1 || { echo "python3 bulunamadi"; exit 1; }

c_cyn=$'\033[36m'; c_grn=$'\033[32m'; c_off=$'\033[0m'
head1(){ printf '\n%s===== %s =====%s\n' "$c_cyn" "$*" "$c_off"; }

cd "$HERE"
mkdir -p "$OUT"

head1 "1/4  Parametre kayit defteri"
"$PY" - <<'PYEOF'
import sys, os
sys.path.insert(0, os.getcwd())
from oligoscout_thermo import params
print(f"{'tablo':<38}{'durum':<26}{'kayit':>6}  sha256")
for row in params.registry_summary():
    print(f"{row['table_id']:<38}{row['verification_status']:<26}"
          f"{row['entries']:>6}  {row['sha256'][:16]}...")
print(f"\nbirlesik digest: {params.registry_digest()}")
PYEOF

head1 "2/4  Dogrulama paketi"
"$PY" -m unittest discover -s tests -t . -v 2>&1 | tail -n 40

head1 "3/4  Demo panel raporu"
"$PY" cli/thermo_report.py --out "$OUT/demo" \
    --ta "$TA" --na "$NA" --mg "$MG" --dntp "$DNTP" \
    --primer-nm "$PRIMER_NM" --tm-window "$TM_WINDOW"

head1 "4/4  Depo paneli (chr22, 32 cift)"
# Bu panel yayindaki 24-nt standart PCR primerleridir (374-424 bp urun),
# qPCR paneli degil. Bu yuzden profil penceresi standart PCR icin acilir.
PANEL_TM_WINDOW="${PANEL_TM_WINDOW:-60:66}"
if [[ -f "$PANEL" ]]; then
    "$PY" cli/thermo_report.py --panel "$PANEL" --out "$OUT/chr22_panel" \
        --ta "$TA" --na "$NA" --mg "$MG" --dntp "$DNTP" \
        --primer-nm "$PRIMER_NM" --tm-window "$PANEL_TM_WINDOW"
else
    echo "panel bulunamadi, atlandi: $PANEL"
fi

printf '\n%sTamamlandi. Ciktilar: %s%s\n' "$c_grn" "$OUT" "$c_off"
