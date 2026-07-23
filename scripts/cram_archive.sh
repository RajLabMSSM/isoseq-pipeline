#!/bin/bash
set -uo pipefail
BAM="$1"; REF="$2"; THREADS="${3:-4}"
SAM=/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/samtools
[ -s "$BAM" ] || { echo "[skip] no BAM: $BAM"; exit 0; }
cram="${BAM%.bam}.cram"

nb=$("$SAM" view -c "$BAM") || { echo "[FAIL bam-read] $BAM"; exit 1; }
if ! { [ -s "$cram" ] && "$SAM" quickcheck "$cram" 2>/dev/null \
       && [ "$("$SAM" view -c -T "$REF" "$cram")" = "$nb" ]; }; then
    "$SAM" view -@ "$THREADS" -T "$REF" -C -o "$cram" "$BAM" || { echo "[FAIL encode] $BAM"; exit 1; }
    "$SAM" index "$cram" || { echo "[FAIL index] $cram"; exit 1; }
fi
nc=$("$SAM" view -c -T "$REF" "$cram")
[ "$nb" = "$nc" ] || { echo "[FAIL count] $BAM bam=$nb cram=$nc"; exit 1; }
"$SAM" quickcheck "$cram" || { echo "[FAIL quickcheck] $cram"; exit 1; }
rm -f "$BAM" "$BAM.bai" "${BAM%.bam}.bai"
echo "[ok] $cram ($nc reads; $(du -h "$cram" | cut -f1))"
