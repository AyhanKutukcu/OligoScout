#!/usr/bin/env bash
# =============================================================================
# OligoScout — tek adimli kurulum ve dogrulama scripti (WSL2 / Ubuntu)
#
# Bu tek script sirayla sunlari yapar:
#   FAZ 0  Ortam raporu (WSL, CPU, RAM, disk)
#   FAZ 1  Derleme bagimliliklari (apt)
#   FAZ 2  CMake >= 3.25 garantisi
#   FAZ 3  Depo + primer3 alt modulu (pinlenmis commit'ler)
#   FAZ 4  primer3 kutuphanesi (make)
#   FAZ 5  OligoScout derlemesi (Ninja, release preset)
#   FAZ 6  Hizli regresyon testleri (genom gerektirmeyenler)
#   FAZ 7  GRCh38.p14 referansi: indirme + normalizasyon + SHA-256 dogrulama
#   FAZ 8  24 kromozom icin PPFM indeksi (paralel) + manifest
#   FAZ 9  Tam test paketi + ozet rapor
#
# Kullanim:
#   bash setup/00_bootstrap.sh              # tum fazlar, kaldigi yerden
#   bash setup/00_bootstrap.sh --from 5     # 5. fazdan itibaren
#   bash setup/00_bootstrap.sh --only 8     # yalniz 8. faz
#   bash setup/00_bootstrap.sh --force      # tamamlanmis fazlari da yeniden calistir
#   bash setup/00_bootstrap.sh --list       # faz listesi
#
# Ayarlar ortam degiskeniyle ezilebilir, ornegin:
#   PRESET=release-avx2 INDEX_PAR=8 bash setup/00_bootstrap.sh --only 8
# =============================================================================

set -Eeuo pipefail

# ----------------------------- Ayarlar ---------------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
REPO="${REPO:-$ROOT/OligoScout}"
LOGDIR="${LOGDIR:-$ROOT/logs}"
STATEDIR="$LOGDIR/state"
DLDIR="${DLDIR:-$ROOT/downloads}"

REPO_URL="https://github.com/AyhanKutukcu/OligoScout.git"
REPO_COMMIT="${REPO_COMMIT:-45161ea2492f5086d122d9abb65bb24ff058b023}"
PRIMER3_COMMIT="${PRIMER3_COMMIT:-345cba9ca22fdb8f9f90da9b57948b81b4866330}"

# NCBI RefSeq GRCh38.p14. build_grch38_ppfm kromozomu FASTA basligindan
# cozumler (NC_ + "Homo sapiens chromosome N," + "GRCh38.p14 Primary Assembly"),
# bu yuzden kaynak NCBI olmak zorundadir. Ensembl/UCSC basliklari araci kirar.
REF_URL="${REF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz}"
REF_GZ="$DLDIR/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
FASTA="${FASTA:-$REPO/data/genomes/GRCh38.p14.fa}"
MANIFEST_REF="$REPO/data/benchmarks/reference/GRCh38.p14.primary24.manifest.tsv"

PRESET="${PRESET:-release}"              # release | release-avx2
BUILDDIR="$REPO/build/$PRESET"
SA_RATE="${SA_RATE:-8}"
INDEXDIR="${INDEXDIR:-$REPO/data/indexes/ppfm}"
JOBS="${JOBS:-$(nproc)}"
INDEX_PAR="${INDEX_PAR:-6}"              # es zamanli kromozom indeksi sayisi
CMAKE_MIN="3.25"

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 \
        chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

# ----------------------------- Yardimcilar -----------------------------------
RUN_TS="$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGDIR" "$STATEDIR" "$DLDIR"
MAIN_LOG="$LOGDIR/bootstrap_${RUN_TS}.log"

c_red=$'\033[31m'; c_grn=$'\033[32m'; c_ylw=$'\033[33m'; c_cyn=$'\033[36m'; c_off=$'\033[0m'

log()  { printf '%s %s\n' "$(date +%H:%M:%S)" "$*"; }
head1(){ printf '\n%s===== %s =====%s\n' "$c_cyn" "$*" "$c_off"; }
ok()   { printf '%s  OK  %s%s\n' "$c_grn" "$*" "$c_off"; }
warn() { printf '%s  !!  %s%s\n' "$c_ylw" "$*" "$c_off"; }
die()  { printf '%s  HATA  %s%s\n' "$c_red" "$*" "$c_off" >&2; exit 1; }

have()  { command -v "$1" >/dev/null 2>&1; }
done_f(){ echo "$STATEDIR/faz_${1}.done"; }
is_done(){ [[ -f "$(done_f "$1")" ]]; }
mark()  { date -Iseconds > "$(done_f "$1")"; }

ver_ge() { printf '%s\n%s\n' "$2" "$1" | sort -V -C; }   # $1 >= $2 ise 0

SUDO=""
need_sudo() {
    if [[ $EUID -eq 0 ]]; then SUDO=""; return 0; fi
    have sudo || die "sudo bulunamadi ve root degilsiniz."
    SUDO="sudo"
    log "Paket kurulumu icin sudo parolasi istenebilir."
    sudo -v
}

# ----------------------------- Argumanlar ------------------------------------
FROM=0; ONLY=""; FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)  FROM="$2"; shift 2 ;;
        --only)  ONLY="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        --list)
            grep -E '^#   FAZ' "$0" | sed 's/^#   //'
            exit 0 ;;
        -h|--help) sed -n '2,30p' "$0" | sed 's/^# \?//'; exit 0 ;;
        *) die "bilinmeyen arguman: $1" ;;
    esac
done

should_run() {
    local n="$1"
    if [[ -n "$ONLY" && "$ONLY" != "$n" ]]; then return 1; fi
    if [[ -z "$ONLY" && "$n" -lt "$FROM" ]]; then return 1; fi
    if [[ $FORCE -eq 0 ]] && is_done "$n"; then
        printf '%s  ATLANDI  FAZ %s daha once tamamlandi (%s)%s\n' \
               "$c_ylw" "$n" "$(cat "$(done_f "$n")")" "$c_off"
        return 1
    fi
    return 0
}

# =============================================================================
# FAZ 0 — Ortam raporu
# =============================================================================
faz0() {
    head1 "FAZ 0 — Ortam"
    if ! grep -qi microsoft /proc/version 2>/dev/null; then
        warn "WSL icinde gorunmuyorsunuz. Script yine de calisir ama yollar farkli olabilir."
    fi
    log "Dagitim : $(. /etc/os-release 2>/dev/null; echo "${PRETTY_NAME:-bilinmiyor}")"
    log "Cekirdek: $(uname -r)"
    log "CPU     : $(nproc) cekirdek"
    log "RAM     : $(awk '/MemTotal/{printf "%.1f GiB", $2/1048576}' /proc/meminfo)"
    log "ROOT    : $ROOT"
    log "REPO    : $REPO"
    df -h "$ROOT" | tail -n1 | awk '{printf "  disk    : %s toplam, %s bos (%s dolu)\n",$2,$4,$5}'

    local free_gb ram_gb
    free_gb=$(df -BG --output=avail "$ROOT" | tail -n1 | tr -dc '0-9')
    if (( free_gb < 40 )); then
        warn "Bos alan ${free_gb} GB. Genom + indeks icin ~40 GB onerilir."
    fi
    ram_gb=$(awk '/MemTotal/{printf "%d", $2/1048576}' /proc/meminfo)
    if (( ram_gb < 16 )); then
        warn "WSL'e ayrilan RAM ${ram_gb} GiB."
        warn "Windows tarafinda %USERPROFILE%\\.wslconfig dosyasina su satirlari ekleyip"
        warn "'wsl --shutdown' calistirin:  [wsl2] / memory=64GB / processors=24"
    fi
    mark 0
}

# =============================================================================
# FAZ 1 — Derleme bagimliliklari
# =============================================================================
faz1() {
    head1 "FAZ 1 — Paketler"
    local pkgs=(build-essential g++ make cmake ninja-build git curl ca-certificates \
                python3 pkg-config time coreutils)
    local missing=()
    local p
    for p in "${pkgs[@]}"; do
        dpkg -s "$p" >/dev/null 2>&1 || missing+=("$p")
    done
    if [[ ${#missing[@]} -eq 0 ]]; then
        ok "Tum paketler zaten kurulu."
    else
        log "Eksik: ${missing[*]}"
        need_sudo
        $SUDO apt-get update
        $SUDO DEBIAN_FRONTEND=noninteractive apt-get install -y "${missing[@]}"
        ok "Paketler kuruldu."
    fi
    g++ --version | head -n1
    mark 1
}

# =============================================================================
# FAZ 2 — CMake >= 3.25
# =============================================================================
faz2() {
    head1 "FAZ 2 — CMake"
    local v=""
    if have cmake; then v="$(cmake --version | head -n1 | awk '{print $3}')"; fi
    if [[ -n "$v" ]] && ver_ge "$v" "$CMAKE_MIN"; then
        ok "cmake $v (>= $CMAKE_MIN)"
    else
        warn "cmake ${v:-yok} yetersiz; yerel kuruluma gecis."
        local cv="3.31.6"
        local tgz="$DLDIR/cmake-${cv}-linux-x86_64.tar.gz"
        local dest="$ROOT/tools/cmake"
        mkdir -p "$dest"
        [[ -s "$tgz" ]] || curl -fL --retry 3 -o "$tgz" \
            "https://github.com/Kitware/CMake/releases/download/v${cv}/cmake-${cv}-linux-x86_64.tar.gz"
        tar -xzf "$tgz" -C "$dest" --strip-components=1
        export PATH="$dest/bin:$PATH"
        echo "export PATH=\"$dest/bin:\$PATH\"" > "$SCRIPT_DIR/cmake_path.env"
        ok "cmake $(cmake --version | head -n1 | awk '{print $3}') ($dest)"
    fi
    mark 2
}

# =============================================================================
# FAZ 3 — Depo ve alt modul
# =============================================================================
faz3() {
    head1 "FAZ 3 — Depo ve primer3 alt modulu"
    if [[ ! -d "$REPO" ]]; then
        log "Depo klonlaniyor..."
        git clone --recurse-submodules "$REPO_URL" "$REPO"
        git -C "$REPO" checkout --detach "$REPO_COMMIT"
    elif [[ ! -d "$REPO/.git" ]]; then
        # Dosyalar var ama git gecmisi yok (tarball kopyasi). Yerinde git'e cevrilir;
        # takip edilmeyen dosyalar (genom, indeks) silinmez.
        log "Git gecmisi yok; yerinde baglaniyor (dosyalar korunur)."
        git -C "$REPO" init -q
        git -C "$REPO" remote add origin "$REPO_URL" 2>/dev/null || true
        git -C "$REPO" fetch --depth 1 origin "$REPO_COMMIT"
        git -C "$REPO" reset --hard FETCH_HEAD
    else
        log "Mevcut depo pinlenmis commit'e hizalaniyor."
        git -C "$REPO" fetch --depth 1 origin "$REPO_COMMIT" || true
        git -C "$REPO" reset --hard "$REPO_COMMIT" 2>/dev/null || \
            git -C "$REPO" reset --hard FETCH_HEAD
    fi
    git -C "$REPO" submodule update --init --recursive
    if [[ -e "$REPO/third_party/primer3/.git" ]]; then
        git -C "$REPO/third_party/primer3" checkout --detach "$PRIMER3_COMMIT" 2>/dev/null || \
            warn "primer3 pinlenmis commit'e alinamadi; alt modul varsayilani kullanilacak."
    fi
    log "HEAD    : $(git -C "$REPO" rev-parse HEAD)"
    log "primer3 : $(git -C "$REPO/third_party/primer3" rev-parse HEAD 2>/dev/null || echo yok)"
    [[ -f "$MANIFEST_REF" ]] || die "Referans manifesti bulunamadi: $MANIFEST_REF"
    ok "Depo hazir."
    mark 3
}

# =============================================================================
# FAZ 4 — primer3
# =============================================================================
faz4() {
    head1 "FAZ 4 — primer3 (pinlenmis bagimlilik)"
    [[ -d "$REPO/third_party/primer3/src" ]] || die "primer3 alt modulu yok. Once FAZ 3."
    make -C "$REPO/third_party/primer3/src" -j"$JOBS"
    ok "primer3 derlendi."
    mark 4
}

# =============================================================================
# FAZ 5 — OligoScout derlemesi
# =============================================================================
faz5() {
    head1 "FAZ 5 — Derleme (preset: $PRESET, -j$JOBS)"
    [[ -f "$SCRIPT_DIR/cmake_path.env" ]] && source "$SCRIPT_DIR/cmake_path.env"
    ( cd "$REPO" && cmake --preset "$PRESET" )
    ( cd "$REPO" && cmake --build --preset "$PRESET" -j"$JOBS" )
    [[ -x "$BUILDDIR/oligoscout" ]]        || die "oligoscout uretilmedi."
    [[ -x "$BUILDDIR/build_grch38_ppfm" ]] || die "build_grch38_ppfm uretilmedi."
    ok "Calistirilabilir: $BUILDDIR/oligoscout , $BUILDDIR/build_grch38_ppfm"
    mark 5
}

# =============================================================================
# FAZ 6 — Hizli testler (genom gerektirmeyenler)
# =============================================================================
faz6() {
    head1 "FAZ 6 — Hizli regresyon testleri"
    [[ -f "$SCRIPT_DIR/cmake_path.env" ]] && source "$SCRIPT_DIR/cmake_path.env"
    local log="$LOGDIR/ctest_hizli_${RUN_TS}.log"
    # ppfm_manifest indeks ciktilarina bagimlidir; FAZ 8'den once calisamaz.
    if ctest --test-dir "$BUILDDIR" --output-on-failure -j"$JOBS" \
             -E 'persistent_genome|grch38|ppfm_batch|ppfm_manifest' 2>&1 | tee "$log" | tail -n 15
    then ok "Hizli testler gecti. Log: $log"
    else warn "Bazi hizli testler basarisiz. Log: $log"
    fi
    mark 6
}

# =============================================================================
# FAZ 7 — Referans genom: indir, normalize et, dogrula
# =============================================================================
faz7() {
    head1 "FAZ 7 — GRCh38.p14 referansi"
    mkdir -p "$(dirname "$FASTA")"

    local regen=0 rc=0
    if [[ -s "$FASTA" ]]; then
        log "Mevcut FASTA: $FASTA ($(du -h "$FASTA" | cut -f1))"
        log "Dogrulaniyor (uzunluk + SHA-256 + baslik uyumu)..."
        set +e
        python3 "$SCRIPT_DIR/normalize_reference.py" \
            --verify-only --output "$FASTA" --manifest "$MANIFEST_REF" \
            2>&1 | tee "$LOGDIR/reference_verify_${RUN_TS}.log"
        rc=${PIPESTATUS[0]}
        set -e
        case "$rc" in
            0) ok "Mevcut referans gecerli ve araca uyumlu." ; mark 7 ; return 0 ;;
            3) warn "Diziler dogru ama basliklar build_grch38_ppfm ile uyumsuz."
               warn "NCBI RefSeq kaynagindan yeniden uretilecek."
               regen=1 ;;
            *) die "Referans dogrulamasi basarisiz (cikis kodu $rc). Log: $LOGDIR" ;;
        esac
    else
        regen=1
    fi

    if (( regen )); then
        if [[ ! -s "$REF_GZ" ]]; then
            log "NCBI RefSeq GRCh38.p14 indiriliyor (~973 MB)..."
            curl -fL --retry 3 -C - -o "$REF_GZ" "$REF_URL"
        fi
        log "Normalizasyon (buyuk harf, IUPAC RYSWKMBDHV -> N, 24 birincil kromozom)..."
        log "Kaynak basliklari aynen korunur; arac kromozomu baslikten cozumler."
        python3 "$SCRIPT_DIR/normalize_reference.py" \
            --input "$REF_GZ" --output "$FASTA" --manifest "$MANIFEST_REF" \
            2>&1 | tee "$LOGDIR/reference_build_${RUN_TS}.log"
        python3 "$SCRIPT_DIR/normalize_reference.py" \
            --verify-only --output "$FASTA" --manifest "$MANIFEST_REF" \
            2>&1 | tee "$LOGDIR/reference_verify_${RUN_TS}.log"
    fi
    ok "Referans hazir: $FASTA"
    mark 7
}

# =============================================================================
# FAZ 8 — PPFM indeksleri (24 kromozom, paralel)
# =============================================================================
build_one() {
    local chrom="$1"
    local final="$INDEXDIR/${chrom}.sa${SA_RATE}.ppfm"
    local part="${final}.part"
    local olog="$LOGDIR/ppfm/build_${chrom}_sa${SA_RATE}.tsv"
    local tlog="$LOGDIR/ppfm/build_${chrom}_sa${SA_RATE}.time.txt"
    if [[ -s "$final" ]]; then echo "EXISTING $chrom"; return 0; fi
    rm -f "$part"
    if ! /usr/bin/time -v "$BUILDDIR/build_grch38_ppfm" \
            "$FASTA" "$chrom" "$part" "$SA_RATE" > "$olog" 2> "$tlog"; then
        echo "FAILED   $chrom (builder hatasi)"; rm -f "$part"; return 1
    fi
    if ! grep -q $'^ALL_CHECKS\tYES$' "$olog"; then
        echo "FAILED   $chrom (ALL_CHECKS YES yok)"; rm -f "$part"; return 1
    fi
    [[ -s "$part" ]] || { echo "FAILED   $chrom (bos indeks)"; return 1; }
    mv "$part" "$final"
    echo "BUILT    $chrom $(stat -c '%s' "$final") bayt"
}

faz8() {
    head1 "FAZ 8 — PPFM indeksleri (sa_rate=$SA_RATE, es zamanli=$INDEX_PAR)"
    [[ -s "$FASTA" ]]                      || die "FASTA yok. Once FAZ 7."
    [[ -x "$BUILDDIR/build_grch38_ppfm" ]] || die "Builder yok. Once FAZ 5."
    mkdir -p "$INDEXDIR" "$LOGDIR/ppfm"

    export -f build_one
    export INDEXDIR SA_RATE LOGDIR BUILDDIR FASTA

    printf '%s\n' "${CHROMS[@]}" \
        | xargs -P "$INDEX_PAR" -I{} bash -c 'build_one "$@"' _ {} \
        | tee "$LOGDIR/ppfm_build_${RUN_TS}.log"

    # Manifest, deponun kendi betiginin ve ppfm_manifest testinin bekledigi
    # yola ve bicime yazilir. Test calisma dizini depo kokudur ve basligin
    # tam olarak "chromosome\tfile_bytes\tsha256\tstatus" olmasini ister
    # (bkz. src/ppfm_manifest.cpp, tests/test_ppfm_manifest.cpp).
    local canon_dir="$REPO/results/benchmarks/grch38_ppfm_batch"
    local manifest="$canon_dir/manifest.sa${SA_RATE}.tsv"
    mkdir -p "$canon_dir"

    local missing=0 chrom f
    printf 'chromosome\tfile_bytes\tsha256\tstatus\n' > "$manifest"
    for chrom in "${CHROMS[@]}"; do
        f="$INDEXDIR/${chrom}.sa${SA_RATE}.ppfm"
        if [[ -s "$f" ]]; then
            printf '%s\t%s\t%s\t%s\n' "$chrom" "$(stat -c '%s' "$f")" \
                   "$(sha256sum "$f" | awk '{print $1}')" "BUILT" >> "$manifest"
        else
            warn "eksik: $chrom"; missing=$((missing+1))
        fi
    done
    (( missing > 0 )) && die "$missing kromozom indeksi uretilemedi. Log: $LOGDIR/ppfm/"
    cp "$manifest" "$INDEXDIR/manifest.sa${SA_RATE}.tsv"
    ok "24 indeks hazir. Manifest: $manifest"
    du -sh "$INDEXDIR" | awk '{print "  toplam indeks boyutu: "$1}'
    mark 8
}

# =============================================================================
# FAZ 9 — Tam test paketi ve ozet
# =============================================================================
faz9() {
    head1 "FAZ 9 — Tam test paketi"
    [[ -f "$SCRIPT_DIR/cmake_path.env" ]] && source "$SCRIPT_DIR/cmake_path.env"
    local log="$LOGDIR/ctest_tam_${RUN_TS}.log"
    if ctest --test-dir "$BUILDDIR" --output-on-failure -j"$JOBS" 2>&1 | tee "$log" | tail -n 20
    then ok "Tum testler gecti."
    else warn "Bazi testler basarisiz. Log: $log"
    fi

    head1 "OZET"
    printf '  depo commit      : %s\n' "$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo '?')"
    printf '  calistirilabilir : %s\n' "$BUILDDIR/oligoscout"
    printf '  referans FASTA   : %s (%s)\n' "$FASTA" "$(du -h "$FASTA" 2>/dev/null | cut -f1)"
    printf '  PPFM indeksleri  : %s (%s)\n' "$INDEXDIR" "$(du -sh "$INDEXDIR" 2>/dev/null | cut -f1)"
    printf '  loglar           : %s\n' "$LOGDIR"
    mark 9
}

# ----------------------------- Calistir --------------------------------------
{
    printf 'OligoScout bootstrap — %s\n' "$(date -Iseconds)"
    for n in 0 1 2 3 4 5 6 7 8 9; do
        if should_run "$n"; then "faz$n"; fi
    done
    printf '\n%sTamamlandi. Ana log: %s%s\n' "$c_grn" "$MAIN_LOG" "$c_off"
} 2>&1 | tee "$MAIN_LOG"
