#!/usr/bin/env bash

set -euo pipefail

K=12

if [[ $# -ne 2 ]]; then
    echo "Kullanım: $0 <jellyfish_index.jf> <primer_5to3>"
    exit 1
fi

INDEX_PATH="$1"
PRIMER=$(printf '%s' "$2" | tr '[:lower:]' '[:upper:]')

if [[ ! -f "$INDEX_PATH" ]]; then
    echo "Hata: İndeks bulunamadı: $INDEX_PATH" >&2
    exit 2
fi

if [[ ! "$PRIMER" =~ ^[ACGT]+$ ]]; then
    echo "Hata: Primer yalnızca A, C, G ve T içermelidir." >&2
    exit 3
fi

if (( ${#PRIMER} < K )); then
    echo "Hata: Primer en az $K baz olmalıdır." >&2
    exit 4
fi

SEED="${PRIMER: -K}"

QUERY_OUTPUT=$(
    jellyfish query \
        "$INDEX_PATH" \
        "$SEED"
)

COUNT=$(
    printf '%s\n' "$QUERY_OUTPUT" |
    awk 'NR == 1 {print $2}'
)

COUNT=${COUNT:-0}

if (( COUNT == 0 )); then
    DECISION="NO_EXACT_SEED"
    NEXT_STEP="SENSITIVE_FALLBACK"
elif (( COUNT <= 10 )); then
    DECISION="LOW_FREQUENCY_SEED"
    NEXT_STEP="FM_INDEX_SEARCH"
elif (( COUNT <= 100 )); then
    DECISION="MODERATE_FREQUENCY_SEED"
    NEXT_STEP="FM_INDEX_SEARCH"
else
    DECISION="HIGH_FREQUENCY_SEED"
    NEXT_STEP="ALTERNATIVE_ANCHOR_OR_FM_INDEX"
fi

printf 'primer\t%s\n' "$PRIMER"
printf 'primer_length\t%s\n' "${#PRIMER}"
printf 'seed_3prime\t%s\n' "$SEED"
printf 'seed_length\t%s\n' "$K"
printf 'canonical_seed_count\t%s\n' "$COUNT"
printf 'decision\t%s\n' "$DECISION"
printf 'next_step\t%s\n' "$NEXT_STEP"
