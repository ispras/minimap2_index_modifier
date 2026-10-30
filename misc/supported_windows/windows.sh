#!/bin/bash
# Forced realignment windows from an aligned BAM: samtools, mawk and sort.
# usage: windows.sh BAM PANEL FAI MQ OUTDIR [samtools_threads]
# writes OUTDIR/mq$MQ.bed, mq$MQ.events.jsonl, support.tsv, params.json
set -euo pipefail
export LC_ALL=C
bam=$1 panel=$2 fai=$3 mq=$4 out=$5 threads=${6:-4}
here=$(dirname "$0")
n=2 w=200 maxlen=1000
mkdir "$out"
samtools view -@ "$threads" -d lt "$bam" | mawk -v MQ="$mq" -f "$here/support.awk" > "$out/support.tsv"
mawk -v N="$n" -v VALIDATE=1 -v FAI="$fai" -v SUP="$out/support.tsv" -f "$here/events.awk" "$fai" "$out/support.tsv" "$panel" > "$out/events.tsv"
sort -t$'\t' -k1,1n -k2,2n -k3,3n "$out/events.tsv" \
    | mawk -v W="$w" -v MAXLEN="$maxlen" -v FAI="$fai" -f "$here/geometry.awk" > "$out/mq$mq.bed"
sort -t$'\t' -k4,4 -k2,2n -k5,5 -k6,6 "$out/events.tsv" \
    | mawk -F'\t' '{printf "{\"chrom\": \"%s\", \"pos\": %d, \"ref\": \"%s\", \"alt\": \"%s\"}\n",$4,$2,$5,$6}' > "$out/mq$mq.events.jsonl"
rm "$out/events.tsv"
printf '{"mq": %d, "n": %d, "w": %d, "max_window": %d}\n' "$mq" "$n" "$w" "$maxlen" > "$out/params.json"
