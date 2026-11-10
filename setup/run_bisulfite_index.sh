#!/usr/bin/env bash
# =============================================================================
# Bisulfit-donusturulmus genom indeksleri (BIS-02/03/04 altyapisi)
#
#   bash setup/run_bisulfite_index.sh
#
# Bisulfit muamelesi metillenmemis sitozini urasile (dizilemede T) cevirir.
# Dort iplik olusur: donusturulmus + ve - iplikler ve bunlarin tamamlayici
# kopyalari. OligoScout'un bir ipligi ancak kendi donusturulmus genomunda
# ararsa dogru sonuc verir. Bu betik referanstan iki donusturulmus genom
# uretir (C->T ve G->A) ve her biri icin PPFM indeksi kurar.
#
# AGIR: tam genom donusumu + 2x24 kromozom-basina indeks kurulumu; ~12-17 GB
# ek disk ister (2 donusturulmus FASTA ~6 GB + 2x indeks ~6 GB). SURE: chr1
# tek is parcacigiyla ~40 dk; 24 kromozom sirali ~8 saat/donusum, iki donusum
# ~16 saat. JOBS ile paralel kurulursa cok kisalir (orn. JOBS=6 -> ~2-3 saat).
# JOBS her is ~birkac GB RAM kullanir; 126 GB makinede JOBS=6-8 guvenli.
#
# Faz 8'in build_all_grch38_ppfm.sh surucusuyle AYNI arayuz:
#   build_grch38_ppfm <FASTA> <chromosome> <output.ppfm> <sa_rate>
# Her kromozom icin tek .ppfm dosyasi + ALL_CHECKS dogrulamasi + manifest.
# Yeniden calistirilabilir: mevcut .ppfm/.fa dosyalari atlanir (kaldigi
# yerden devam eder) — bir calismayi durdurup JOBS ile yeniden baslatabilirsin.
#
# Ayarlar: PRESET (release), SA_RATE (8), JOBS (1), CHROMS (bos=24), OUTDIR
#
# NOT: build_grch38_ppfm kromozomu FASTA basligindan cozer ve yalniz NCBI
# RefSeq basliklarini (NC_...) kabul eder. Donusturulmus FASTA basliklari
# KORUNUR; yalnizca dizi harfleri donusur.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
REPO="$ROOT/OligoScout"
PRESET="${PRESET:-release}"
BUILDER="$REPO/build/$PRESET/build_grch38_ppfm"
GENOME="$REPO/data/genomes/GRCh38.p14.fa"
OUTDIR="${OUTDIR:-$REPO/data/indexes/bisulfite}"
SA_RATE="${SA_RATE:-8}"          # Faz 8 ana indeksle ayni (sa8)
JOBS="${JOBS:-1}"                # es zamanli kromozom kurulumu (paralellik)

[[ -s "$GENOME" ]] || { echo "Referans yok: $GENOME (Faz 7)"; exit 1; }
[[ -x "$BUILDER" ]] || {
    echo "build_grch38_ppfm yok; once: RUN_SETUP.cmd --only 5"; exit 1; }

# Kanonik kromozomlar — Faz 8 ana indeksle AYNI kume (24 shard; chrM haric,
# o Python tarafinda islenir). build_grch38_ppfm kromozomu AD ile cozer:
# arguman "chr1" iken FASTA basligindaki "chromosome 1" ile eslesir, bu
# nedenle donusturulmus FASTA basliklari KORUNMALIDIR (awk yalniz diziyi
# donusturur).
# CHROMS ile alt kume secilebilir (hizli duman testi), orn:
#   CHROMS="chrY" bash setup/run_bisulfite_index.sh   # tek kucuk kromozom
# Bos birakilirsa tam 24 kromozomluk kume kullanilir.
if [[ -n "${CHROMS:-}" ]]; then
    read -r -a chromosomes <<< "$CHROMS"
else
    chromosomes=()
    for n in $(seq 1 22); do chromosomes+=("chr${n}"); done
    chromosomes+=("chrX" "chrY")
fi

mkdir -p "$OUTDIR"

# --- 1) Iki donusturulmus genom (basliklar korunur, yalniz dizi harfleri) ---
CT="$OUTDIR/GRCh38.CT.fa"     # C->T (ust iplik bisulfit donusumu)
GA="$OUTDIR/GRCh38.GA.fa"     # G->A (alt iplik bisulfit donusumu)

if [[ ! -s "$CT" ]]; then
    echo "[1/3] C->T donusturulmus genom uretiliyor"
    awk '/^>/{print; next}{gsub(/C/,"T"); gsub(/c/,"t"); print}' "$GENOME" \
        > "$CT.part" && mv "$CT.part" "$CT"
else
    echo "[1/3] C->T genom mevcut, atlaniyor: $CT"
fi
if [[ ! -s "$GA" ]]; then
    echo "[2/3] G->A donusturulmus genom uretiliyor"
    awk '/^>/{print; next}{gsub(/G/,"A"); gsub(/g/,"a"); print}' "$GENOME" \
        > "$GA.part" && mv "$GA.part" "$GA"
else
    echo "[2/3] G->A genom mevcut, atlaniyor: $GA"
fi

# --- 2) Her donusum icin kromozom-basina PPFM (Faz 8 surucusuyle ayni) ------
echo "[3/3] PPFM indeksleri kuruluyor (sa_rate=$SA_RATE, JOBS=$JOBS)"

# Tek kromozom kurucu — xargs ile paralel cagrilabilmesi icin export edilir.
# Ortam: BUILDER, SA_RATE, CONV, FA, IDIR (export) uzerinden alinir.
build_chrom() {
    local chromosome="$1"
    local final="$IDIR/${chromosome}.sa${SA_RATE}.ppfm"
    local part="${final}.part"
    local log="$IDIR/build_${chromosome}.log"
    if [[ -s "$final" ]]; then
        echo "  $CONV/$chromosome: mevcut, atlaniyor"
        return 0
    fi
    echo "  $CONV/$chromosome: kuruluyor ..."
    rm -f "$part" "$log"
    # Pozisyonel arayuz: <FASTA> <chromosome> <output.ppfm> <sa_rate>
    if ! "$BUILDER" "$FA" "$chromosome" "$part" "$SA_RATE" > "$log" 2>&1; then
        echo "  HATA: $CONV/$chromosome kurulumu basarisiz. Log: $log"
        rm -f "$part"; return 1
    fi
    if ! grep -q $'^ALL_CHECKS\tYES$' "$log"; then
        echo "  HATA: $CONV/$chromosome ALL_CHECKS YES yok. Log: $log"
        rm -f "$part"; return 1
    fi
    [[ -s "$part" ]] || { echo "  HATA: bos PPFM $CONV/$chromosome"; return 1; }
    mv "$part" "$final"
    echo "  $CONV/$chromosome: BUILT"
}
export -f build_chrom
export BUILDER SA_RATE

build_one_conversion() {
    local conv="$1" fa="$2"
    local idir="$OUTDIR/$conv"
    mkdir -p "$idir"
    echo
    echo "======== $conv donusumu ($fa), JOBS=$JOBS ========"

    # Kromozomlari JOBS kadar es zamanli kur (xargs havuzu). Buyukten kucuge
    # siralamak kuyrugu dengeler. Her kurucu kendi .ppfm'ini yazar; hata
    # durumunda ilgili log'da kalir.
    export CONV="$conv" FA="$fa" IDIR="$idir"
    # || true: bir kurulum basarisiz olsa bile asagidaki eksik-denetimi
    # calissin ve hangi kromozomun basarisiz oldugunu bildirsin.
    printf '%s\n' "${chromosomes[@]}" \
        | xargs -P "$JOBS" -I{} bash -c 'build_chrom "$@"' _ {} || true

    # Havuz bitti — manifest'i final dosyalardan birlestir; eksik varsa hata.
    local manifest_tmp="$idir/manifest.sa${SA_RATE}.tsv.tmp"
    local manifest="$idir/manifest.sa${SA_RATE}.tsv"
    printf "chromosome\tfile_bytes\tsha256\tstatus\n" > "$manifest_tmp"
    local missing=0
    for chromosome in "${chromosomes[@]}"; do
        local final="$idir/${chromosome}.sa${SA_RATE}.ppfm"
        if [[ ! -s "$final" ]]; then
            echo "  HATA: $conv/$chromosome uretilemedi (log: $idir/build_${chromosome}.log)"
            missing=1; continue
        fi
        local bytes sha
        bytes=$(stat -c '%s' "$final")
        sha=$(sha256sum "$final" | awk '{print $1}')
        printf "%s\t%s\t%s\t%s\n" "$chromosome" "$bytes" "$sha" "BUILT" \
            >> "$manifest_tmp"
    done
    if [[ "$missing" -ne 0 ]]; then
        echo "  $conv: bazi kromozomlar eksik; manifest yazilmadi."
        rm -f "$manifest_tmp"; exit 1
    fi
    mv "$manifest_tmp" "$manifest"
    echo "  $conv manifest: $manifest"
}

build_one_conversion CT "$CT"
build_one_conversion GA "$GA"

echo
echo "Bitti. Donusturulmus indeksler:"
echo "  $OUTDIR/CT/chr*.sa${SA_RATE}.ppfm  (+ manifest.sa${SA_RATE}.tsv)"
echo "  $OUTDIR/GA/chr*.sa${SA_RATE}.ppfm  (+ manifest.sa${SA_RATE}.tsv)"
echo "Sonraki adim: panel_qc'ye bisulfit tarama yolunu baglamak (Python)."
