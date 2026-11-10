#!/usr/bin/env bash
# =============================================================================
# Primer tasarimi (ozgulluk taramasi ONCESI)
#
#   TARGET=/yol/hedef.fa bash setup/run_design.sh
#   TARGET=ACGT...       bash setup/run_design.sh
#
#   PROFILE=ctdna TARGET=/yol/hedef.fa bash setup/run_design.sh
#
# Ayarlar: PROFILE, PRODUCT (90:160), TM (58:62), TM_OPT, LENGTH (18:27),
#          GC (40:60), TA, NA, MG, DNTP, PRIMER_NM, TOP, PROBE, REQUIRE_TH03
#
# PROBE=1 verilirse (ddpcr profilinde kendiliginden) her cift icin hidroliz
# probu secilir ve PRB-01..PRB-05 uygulanir. Prob icin TARGET bir FASTA
# dosyasi olmali; amplikon dizisi oradan kesilir.
#
# PROFILE verilirse kisitlar profilden gelir; ayrica verdiginiz bir ayar
# (orn. PRODUCT) profili ezer.
# Profiller: qpcr_sybr rt_qpcr ctdna ddpcr sanger ngs_amplicon
#            standard_pcr bisulfite
#
# Cikti design_panel.tsv dosyasini dogrudan RUN_PANEL.cmd'ye verebilirsiniz:
#   PANEL=<...>/design_panel.tsv bash setup/run_panel.sh
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
PY_DIR="$ROOT/reference_py"
OUT_DEFAULT="$PY_DIR/out/design"

TARGET="${TARGET:-}"
if [[ -z "$TARGET" ]]; then
    echo "TARGET verilmedi."
    echo "Ornek: TARGET=/mnt/c/.../hedef.fa bash setup/run_design.sh"
    exit 2
fi

PROFILE="${PROFILE:-}"
# Her profil kendi klasorune yazar; farkli profiller birbirini ezmez.
if [[ -n "$PROFILE" ]]; then OUT_DEFAULT="$PY_DIR/out/design_$PROFILE"; fi
OUT="${OUT:-$OUT_DEFAULT}"
# Profil yoksa klasik varsayilanlar devreye girer.
if [[ -z "$PROFILE" ]]; then
    PRODUCT="${PRODUCT:-90:160}"; TM="${TM:-58:62}"; TM_OPT="${TM_OPT:-60}"
    LENGTH="${LENGTH:-18:27}";    GC="${GC:-40:60}"
fi

# set -e altinda `[[ ... ]] && ...` kalibi, test basarisiz olunca betigi
# dusurur; bu yuzden acik if bloklari kullaniliyor.
ARGS=()
if [[ -n "$PROFILE" ]];     then ARGS+=(--profile "$PROFILE"); fi
if [[ -n "${PRODUCT:-}" ]]; then ARGS+=(--product "$PRODUCT"); fi
if [[ -n "${TM:-}" ]];      then ARGS+=(--tm "$TM"); fi
if [[ -n "${TM_OPT:-}" ]];  then ARGS+=(--tm-opt "$TM_OPT"); fi
if [[ -n "${LENGTH:-}" ]];  then ARGS+=(--length "$LENGTH"); fi
if [[ -n "${GC:-}" ]];      then ARGS+=(--gc "$GC"); fi
if [[ "${REQUIRE_TH03:-0}" == "1" ]]; then ARGS+=(--require-th03); fi

TA="${TA:-60}"; NA="${NA:-50}"; MG="${MG:-1.5}"; DNTP="${DNTP:-0.6}"
PRIMER_NM="${PRIMER_NM:-200}"; TOP="${TOP:-20}"

mkdir -p "$OUT"
cd "$PY_DIR"

python3 cli/design_primers.py \
    --target "$TARGET" --out "$OUT" \
    ${ARGS[@]+"${ARGS[@]}"} \
    --ta "$TA" --na "$NA" --mg "$MG" --dntp "$DNTP" \
    --primer-nm "$PRIMER_NM" --top "$TOP"

# ---------------------------------------------------------------------------
# Prob modulu (PRB-01..PRB-05). Genom taramasi gerektirmez.
# ---------------------------------------------------------------------------
PROBE="${PROBE:-}"
if [[ -z "$PROBE" && "$PROFILE" == "ddpcr" ]]; then PROBE=1; fi
if [[ "$PROBE" == "1" ]]; then
    if [[ ! -f "$TARGET" ]]; then
        echo
        echo "PROBE=1 verildi ama TARGET bir dosya degil; amplikon dizisi"
        echo "cikarilamadigi icin prob adimi atlandi."
    else
        printf '\n===== Prob modulu (PRB-01..PRB-05) =====\n'
        python3 cli/probe_qc.py \
            --design "$OUT" --target "$TARGET" --out "$OUT" \
            --profile "${PROFILE:-ddpcr}" \
            --ta "$TA" --na "$NA" --mg "$MG" --dntp "$DNTP" \
            --primer-nm "$PRIMER_NM" \
            ${PROBE_LENGTH:+--probe-length "$PROBE_LENGTH"} \
            ${TM_DELTA:+--tm-delta "$TM_DELTA"} \
            ${PROBES:+--probes "$PROBES"}
    fi
fi

echo
echo "Ozgulluk taramasi icin:"
REL="${OUT#"$ROOT/"}"
echo "  RUN_PANEL.cmd ${REL//\//\\}\\design_panel.tsv"
