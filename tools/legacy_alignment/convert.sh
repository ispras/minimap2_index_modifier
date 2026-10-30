#!/usr/bin/env bash
# Convert a hap1-coordinate BAM (fixmate + sort) to GRCh38: leviosam2 lift, bamfix,
# mate completion, remap of deferred pairs to GRCh38, reconcile, merge, fixmate, sort, markdup.
# usage: convert.sh -a BAM -c CLFT -x YAML -f GRCh38.fa -i GRCh38.mmi -r READ_GROUP -o OUTDIR [-t N] [-m SORT_MEM]
# writes OUTDIR/final.bam, final.bam.bai, duplicate_metrics.txt
set -euo pipefail

T=8 MEM=768M
while getopts "a:c:x:f:i:r:o:t:m:" o; do
  case "$o" in
    a) BAM=$OPTARG ;; c) CLFT=$OPTARG ;; x) YAML=$OPTARG ;; f) FASTA=$OPTARG ;;
    i) MMI=$OPTARG ;; r) RG=$OPTARG ;; o) OUT=$OPTARG ;; t) T=$OPTARG ;; m) MEM=$OPTARG ;;
    *) sed -n 4p "$0" >&2; exit 2 ;;
  esac
done
here=$(cd "$(dirname "$0")" && pwd)
bamfix=${BAMFIX:-bamfix}
st=$((T > 1 ? T - 1 : 1))
mkdir "$OUT"
cd "$OUT"

leviosam2 lift -t "$T" -T 512 -G 500 -C "$CLFT" -x "$YAML" -m -f "$FASTA" \
    -S lifted -S clipped_frac:0.1 -O bam -a "$BAM" -p lift
for part in committed deferred; do
  "$bamfix" "lift-$part.bam" "$part-ok.bam" "$part-bad.bam"
done
samtools cat -o deferred.bam deferred-ok.bam deferred-bad.bam committed-bad.bam
samtools sort -@ "$st" -m "$MEM" -T committed-sorttmp -o committed-sorted.bam committed-ok.bam
bash "$here/collx.sh" -a committed-sorted.bam -b deferred.bam -p paired -t 1 -m 256M

minimap2 -x sr -R "$RG" -t "$T" "$MMI" paired-deferred-R1.fq.gz paired-deferred-R2.fq.gz -L -a \
  | samtools view -b -F 0x900 -o remapped-primary.bam -
samtools view -b -F 0x900 -o source-primary.bam paired-deferred.bam
leviosam2 reconcile -m -s source:source-primary.bam -s target:remapped-primary.bam -o reconciled.bam

samtools merge -@ "$st" -c -u -o merged.bam paired-committed.bam reconciled.bam
samtools collate -@ "$st" -u -O merged.bam collate-tmp > collated.bam
samtools fixmate -@ "$st" -m -u collated.bam fixmate.bam
samtools sort -@ "$st" -m "$MEM" -T final-sorttmp -o sorted.bam fixmate.bam
samtools markdup -@ "$st" -d 100 -f duplicate_metrics.txt sorted.bam final.bam
samtools index -@ "$st" final.bam
