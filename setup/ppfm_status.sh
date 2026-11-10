#!/usr/bin/env bash
# =============================================================================
# PPFM indeks uretiminin canli durumu.
#
#   bash setup/ppfm_status.sh          tek seferlik rapor
#   bash setup/ppfm_status.sh -w       10 saniyede bir yenile
#
# Builder ilerleme yazmaz; ozet satirlarini ve .part dosyasini ancak bitince
# yazar. Bu yuzden ilerlemeyi islem suresi, CPU ve RSS'ten okuyoruz.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
REPO="${REPO:-$ROOT/OligoScout}"
INDEXDIR="${INDEXDIR:-$REPO/data/indexes/ppfm}"
LOGDIR="${LOGDIR:-$ROOT/logs}"
SA_RATE="${SA_RATE:-8}"

report() {
    printf '\n===== %s =====\n' "$(date '+%F %T')"

    printf '\n-- calisan builder islemleri --\n'
    if pgrep -x build_grch38_ppfm >/dev/null 2>&1; then
        ps -o pid=,etime=,pcpu=,rss=,args= -C build_grch38_ppfm \
          | awk '{
                rss=$4/1048576;
                chrom="?";
                for(i=5;i<=NF;i++) if($i ~ /^chr[0-9XY]+$/) chrom=$i;
                printf "  %-8s pid %-7s sure %-10s cpu %5s%%  bellek %6.2f GB\n",
                       chrom,$1,$2,$3,rss
            }'
    else
        echo "  (calisan builder yok)"
    fi

    printf '\n-- tamamlanan indeksler --\n'
    local done_n=0 total_bytes=0 f
    shopt -s nullglob
    for f in "$INDEXDIR"/chr*.sa${SA_RATE}.ppfm; do
        done_n=$((done_n+1))
        total_bytes=$((total_bytes + $(stat -c '%s' "$f")))
        printf '  %-28s %8.2f GB\n' "$(basename "$f")" \
               "$(echo "$(stat -c '%s' "$f")" | awk '{print $1/1073741824}')"
    done
    for f in "$INDEXDIR"/chr*.part; do
        printf '  %-28s %8.2f GB  (yaziliyor)\n' "$(basename "$f")" \
               "$(echo "$(stat -c '%s' "$f")" | awk '{print $1/1073741824}')"
    done
    shopt -u nullglob
    printf '  toplam: %d / 24 kromozom, %.2f GB\n' \
           "$done_n" "$(echo "$total_bytes" | awk '{print $1/1073741824}')"

    printf '\n-- sistem --\n'
    free -g | awk '/Mem:/{printf "  bellek: %d GB kullanilan / %d GB toplam, %d GB bos\n",$3,$2,$7}'
    free -g | awk '/Swap:/{if($2>0) printf "  takas : %d GB / %d GB kullanimda\n",$3,$2}'
    df -h "$ROOT" | tail -n1 | awk '{printf "  disk  : %s bos (%s dolu)\n",$4,$5}'

    printf '\n-- son biten kromozomlar --\n'
    local log
    for log in $(ls -t "$LOGDIR"/ppfm/build_*_sa${SA_RATE}.tsv 2>/dev/null | head -n 3); do
        if [[ -s "$log" ]]; then
            printf '  %s\n' "$(basename "$log")"
            grep -E 'chromosome|sequence_length|build_seconds|ppfm_bytes|ALL_CHECKS' \
                 "$log" 2>/dev/null | sed 's/^/    /'
        fi
    done
    [[ -z "$(ls -A "$INDEXDIR" 2>/dev/null)" ]] && \
        echo "  (henuz hicbir kromozom bitmedi)"
}

if [[ "${1:-}" == "-w" ]]; then
    while true; do
        clear
        report
        sleep 10
    done
else
    report
fi
