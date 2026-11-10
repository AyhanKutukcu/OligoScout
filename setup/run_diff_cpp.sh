#!/usr/bin/env bash
# =============================================================================
# Python <-> C++ diferansiyel testi (native_thermo cekirdegi)
#
#   bash setup/run_diff_cpp.sh
#
# Ayarlar: PRESET, SEQUENCES, OUT, NA, MG, DNTP, PRIMER_NM, TA, SALT_MODEL
#
# Adimlar:
#   1. oligoscout_thermo_dump hedefini derler (CMakeLists'te yeni hedef).
#   2. Dizi listesini C++ cekirdekten gecirir, TSV dokumu alir.
#   3. reference_py/cli/diff_cpp.py ayni dizileri Python'dan gecirip
#      sutun sutun karsilastirir.
#
# Cikis kodu: cekirdek sutunlarda fark varsa 1. Dimer sutunlarindaki fark
# BEKLENIR: native_thermo.hpp'ye gore C++ dimer modeli yalniz kesintisiz
# Watson-Crick helisini tariyor; Python tarafinda tekli ic uyumsuzluk ve
# tek bazlik bulge da uretim yolunda. Rapor bunu ayri sinif olarak yazar.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
REPO="$ROOT/OligoScout"
PY_DIR="$ROOT/reference_py"

PRESET="${PRESET:-release}"
BUILD_DIR="$REPO/build/$PRESET"
OUT="${OUT:-$PY_DIR/out/diff_cpp}"
SEQUENCES="${SEQUENCES:-$OUT/diff_cpp_sequences.tsv}"
BINARY="$BUILD_DIR/oligoscout_thermo_dump"

mkdir -p "$OUT"

# --- 1. girdi listesi -------------------------------------------------------
if [[ ! -s "$SEQUENCES" ]]; then
    echo "[1/3] girdi listesi uretiliyor: $SEQUENCES"
    cd "$PY_DIR"
    python3 cli/diff_cpp_sequences.py --out "$SEQUENCES"
else
    echo "[1/3] girdi listesi hazir: $SEQUENCES"
fi

# --- 2. C++ hedefini derle --------------------------------------------------
echo "[2/3] oligoscout_thermo_dump derleniyor ($PRESET)"
# Uyumsuzluk/bulge tablolari params.py'den uretilir; once guncelligini dogrula.
python3 "$REPO/tools/gen_thermo_tables.py" --check || {
    echo "  tablo basligi params.py ile tutmuyor; yeniden uretiliyor"
    python3 "$REPO/tools/gen_thermo_tables.py"
}
if [[ ! -d "$BUILD_DIR" ]]; then
    echo "  yapilandirma yok: $BUILD_DIR"
    echo "  once: RUN_SETUP.cmd --only 5"
    exit 1
fi
cmake --build "$BUILD_DIR" --target oligoscout_thermo_dump -j "$(nproc)"

if [[ ! -x "$BINARY" ]]; then
    echo "  derleme sonrasi ikili bulunamadi: $BINARY"
    exit 1
fi

# --- 3. dokum ve karsilastirma ---------------------------------------------
echo "[3/3] diferansiyel karsilastirma"
cd "$PY_DIR"
python3 cli/diff_cpp.py \
    --sequences "$SEQUENCES" \
    --binary "$BINARY" \
    --out "$OUT" \
    --na "${NA:-50}" --mg "${MG:-1.5}" --dntp "${DNTP:-0.6}" \
    --primer-nm "${PRIMER_NM:-200}" --ta "${TA:-60}" \
    --salt-model "${SALT_MODEL:-owczarzy_2008}"
