#!/usr/bin/env bash
# Complete mate pairs between a committed BAM (-a) and a deferred BAM (-b) after a lift:
# deferred reads whose mate is only in committed get that mate moved to deferred.
# Outputs <prefix>-committed.bam, <prefix>-deferred.bam (name-sorted), paired FASTQ of the
# deferred reads and <prefix>.collx.stats (key=value counts).
set -u
set -o pipefail

T=24
TMPD=""
A=""; B=""; P=""
SORTMEM="16G"
while getopts "a:b:p:t:T:m:" o; do
  case "$o" in
    a) A=$OPTARG ;;
    b) B=$OPTARG ;;
    p) P=$OPTARG ;;
    t) T=$OPTARG ;;
    T) TMPD=$OPTARG ;;
    m) SORTMEM=$OPTARG ;;
    *) echo "usage: $0 -a <bam> -b <deferred bam> -p <prefix> [-t N] [-T dir] [-m sort_mem]" >&2; exit 2 ;;
  esac
done
[ -n "$A" ] && [ -n "$B" ] && [ -n "$P" ] || {
  echo "[collx] -a, -b and -p are required" >&2; exit 2; }
[ -s "$A" ] || { echo "[collx] ERROR: missing or empty $A" >&2; exit 1; }
[ -s "$B" ] || { echo "[collx] ERROR: missing or empty $B" >&2; exit 1; }

TMPD=${TMPD:-$P.collxtmp}
for target in "$TMPD" "$P.collx.stats" "$P-committed.bam" "$P-deferred.bam" "$P-deferred-R1.fq.gz" "$P-deferred-R2.fq.gz"; do
  [ ! -e "$target" ] && [ ! -L "$target" ] || { echo "Occupied output: $target" >&2; exit 1; }
done
mkdir -p "$TMPD" || exit 1
S=$TMPD
ST=$P.collx.stats
: > "$ST"

t0=$(date +%s)
say () { echo "[collx $(date +%H:%M:%S) +$(( $(date +%s) - t0 ))s] $*" >&2; }
stat_kv () { echo "$1=$2" >> "$ST"; }
die () { echo "[collx] ERROR: $*" >&2; exit 1; }

say "input A=$A"
say "input B=$B"
say "output prefix $P, threads $T, temp dir $S"

HA=$(samtools view -H "$A" | grep -c '^@SQ') || die "cannot read header of A"
HB=$(samtools view -H "$B" | grep -c '^@SQ') || die "cannot read header of B"
MA=$(samtools view -H "$A" | grep '^@SQ' | md5sum | cut -d' ' -f1)
MB=$(samtools view -H "$B" | grep '^@SQ' | md5sum | cut -d' ' -f1)
stat_kv sq_lines_a "$HA"; stat_kv sq_lines_b "$HB"
[ "$MA" = "$MB" ] || die "@SQ lists of A and B differ ($HA vs $HB lines)"
say "step 0: headers match, $HA @SQ lines"
samtools view -H "$B" > "$S/hdr.sam" || die "cannot read header of B"
NA=$(samtools view -c -@ "$T" "$A") || die "cannot count records in A"
stat_kv input_committed_records "$NA"
say "step 0: records in A: $NA"

say "step 1: deferred read names"
samtools view -@ "$T" -F 0xF00 "$B" \
| mawk -F'\t' -v cnt="$S/nprim.cnt" '
    { f=$2+0;
      if (int(f/64)%2)       seg=1;
      else if (int(f/128)%2) seg=2;
      else { print "[collx] record " $1 " is not paired (flag=" f ")" > "/dev/stderr"; exit 3 }
      print $1 "\t" seg }
    END { print NR > cnt }' \
| LC_ALL=C sort --parallel="$T" -S "$SORTMEM" -T "$S" -u > "$S/def.ns"
rc=$?
[ $rc -eq 0 ] || die "step 1 failed (exit $rc)"
NPRIM=$(cat "$S/nprim.cnt")
NUNIQ=$(wc -l < "$S/def.ns")
stat_kv deferred_primary_records "$NPRIM"
stat_kv deferred_uniq_name_seg "$NUNIQ"
[ "$NPRIM" = "$NUNIQ" ] || die "duplicate (name, segment) in deferred BAM: $NPRIM records, $NUNIQ unique; leviosam2 collate semantics would differ"
say "step 1: primary records in B: $NPRIM"

say "step 2: split names into pairs within B and singletons"
mawk -F'\t' -v both="$S/both.n" -v unp="$S/unp.n" '
  { if (NR>1 && $1!=p) { if (c==2) print p > both; else print p > unp }
    if ($1!=p) { p=$1; c=1 } else c++ }
  END { if (NR>0) { if (c==2) print p > both; else print p > unp } }' "$S/def.ns" \
  || die "step 2 failed"
touch "$S/both.n" "$S/unp.n"
NBOTH=$(wc -l < "$S/both.n"); NUNP=$(wc -l < "$S/unp.n")
stat_kv names_paired_in_deferred "$NBOTH"
stat_kv names_unpaired_in_deferred "$NUNP"
[ $(( 2*NBOTH + NUNP )) -eq "$NPRIM" ] || die "step 2: 2*$NBOTH+$NUNP != $NPRIM"
say "step 2: pairs within B: $NBOTH, singletons (map U): $NUNP"

say "step 3: select committed records with names in U"
samtools view -@ "$T" -b -N "$S/unp.n" -o "$S/csel.bam" "$A" || die "step 3 (samtools view -N) failed"
NCSEL=$(samtools view -c -@ "$T" "$S/csel.bam") || die "cannot count csel"
stat_kv committed_records_in_U "$NCSEL"
say "step 3: committed records with names in U: $NCSEL (remaining $(( NA - NCSEL )))"

say "step 4: primary mates and trailing non-primary records"
samtools view -@ "$T" -b -F 0xF00 -o "$S/cmatch.bam" "$S/csel.bam" || die "step 4a failed"
samtools view -@ "$T" "$S/cmatch.bam" | cut -f1 | LC_ALL=C sort --parallel="$T" -S "$SORTMEM" -T "$S" > "$S/matched.n" || die "step 4b failed"
NMATCH=$(wc -l < "$S/matched.n")
NMATCHU=$(LC_ALL=C uniq "$S/matched.n" | wc -l)
stat_kv matched_from_committed "$NMATCH"
[ "$NMATCH" = "$NMATCHU" ] || die "step 4: more than one primary committed record per name in U ($NMATCH vs $NMATCHU)"
stat_kv remaining_in_map $(( NUNP - NMATCH ))
[ "$NUNP" -eq "$NMATCH" ] || die "Incomplete deferred pairs; refusing dropped primary ends"
say "step 4: mates taken from committed $NMATCH, left in map $(( NUNP - NMATCH ))"

samtools view -@ "$T" "$S/csel.bam" \
| mawk -F'\t' '{ f=$2+0; if (int(f/256)%16) print $1 }' \
| LC_ALL=C sort --parallel="$T" -S "$SORTMEM" -T "$S" -u > "$S/npnames.n" || die "step 4c failed"
NNP=$(wc -l < "$S/npnames.n")
stat_kv names_with_nonprimary_in_U "$NNP"

if [ "$NNP" -gt 0 ]; then
  samtools view -@ "$T" -h "$S/csel.bam" \
  | mawk -F'\t' -v want="$S/npnames.n" -v cls="$S/npclass.txt" '
      BEGIN { while ((getline l < want) > 0) w[l]=1 }
      /^@/ { print; next }
      { if (!($1 in w)) next
        f=$2+0
        if (int(f/256)%16) {
          if ($1 in seen) { after++; print }
          else { before++; pend[$1]++ }
        } else { seen[$1]=1 }
      }
      END {
        for (k in pend) if (!(k in seen)) { noprim += pend[k]; before -= pend[k] }
        printf "after=%d\nbefore=%d\nnoprim=%d\n", after+0, before+0, noprim+0 > cls }' \
  | samtools view -@ "$T" -b -o "$S/cextra.bam" - || die "step 4d failed"
  NEXTRA=$(samtools view -c -@ "$T" "$S/cextra.bam") || die "cannot count cextra"
  NBEFORE=$(sed -n 's/^before=//p' "$S/npclass.txt")
  NNOPRIM=$(sed -n 's/^noprim=//p' "$S/npclass.txt")
else
  NEXTRA=0; NBEFORE=0; NNOPRIM=0
fi
stat_kv committed_nonprimary_before_primary "$NBEFORE"
stat_kv committed_nonprimary_without_primary "$NNOPRIM"
stat_kv committed_tail_nonprimary "$NEXTRA"
say "step 4: non-primary records with names in U: after primary $NEXTRA, before $NBEFORE, without primary $NNOPRIM"
if [ "$NEXTRA" -gt 0 ]; then
  say "WARNING: step 7 uses samtools merge, which removes and re-appends the RG tag,"
  say "         moving RG:Z to the end of each record. Content is unchanged,"
  say "         but records are no longer byte-identical to leviosam2 output."
fi

say "step 5: building $P-deferred.bam"
LC_ALL=C sort --parallel="$T" -S "$SORTMEM" -T "$S" -m "$S/both.n" "$S/matched.n" > "$S/defnames.n" || die "step 5a failed"
samtools view -@ "$T" -b -N "$S/defnames.n" -o "$S/dsel.bam" "$B" || die "step 5b failed"
NDSEL=$(samtools view -c -@ "$T" "$S/dsel.bam") || die "cannot count dsel"
[ "$NDSEL" -eq $(( 2*NBOTH + NMATCH )) ] || die "step 5: selected $NDSEL from B, expected $(( 2*NBOTH + NMATCH ))"

samtools cat -h "$S/hdr.sam" -o "$S/dpre.bam" "$S/dsel.bam" "$S/cmatch.bam" || die "step 5c (cat) failed"
samtools sort -n -@ "$T" -m 2G -T "$S/nsort" -o "$P-deferred.bam" "$S/dpre.bam" || die "step 5d (sort -n) failed"
NDEF=$(samtools view -c -@ "$T" "$P-deferred.bam") || die "cannot count deferred"
stat_kv out_deferred_records "$NDEF"
[ "$NDEF" -eq $(( 2*NBOTH + 2*NMATCH )) ] || die "step 5: deferred has $NDEF records, expected $(( 2*NBOTH + 2*NMATCH ))"
say "step 5: $P-deferred.bam has $NDEF records"

say "step 6: FASTQ from the deferred BAM"
samtools fastq -@ "$T" -n -F 0xF00 \
  -1 "$P-deferred-R1.fq.gz" -2 "$P-deferred-R2.fq.gz" \
  -0 "$S/other.fq" -s "$S/single.fq" "$P-deferred.bam" > "$S/fastq.log" 2>&1 \
  || { cat "$S/fastq.log" >&2; die "step 6 (samtools fastq) failed"; }
[ -s "$S/single.fq" ] && die "step 6: samtools fastq produced singletons; deferred BAM is not fully paired"
[ -s "$S/other.fq" ] && die "step 6: samtools fastq found records without the R1/R2 bit"
say "step 6: $(grep -c . "$S/fastq.log" || true) lines in samtools fastq log"

say "step 7: building $P-committed.bam"
# Names in U are routed away; everything else streams through the FIFO into the writer.
mkfifo "$S/cbase.pipe" || die "cannot create FIFO"
if [ "$NEXTRA" -gt 0 ]; then
  samtools merge -c -p -@ "$T" -f -o "$P-committed.bam" "$S/cbase.pipe" "$S/cextra.bam" &
else
  samtools view -@ "$T" -b -o "$P-committed.bam" "$S/cbase.pipe" &
fi
CPID=$!
samtools view -@ "$T" -u -N "$S/unp.n" -o /dev/null -U "$S/cbase.pipe" "$A"
rc7=$?
wait $CPID; rc7b=$?

[ $rc7 -eq 0 ] || die "step 7: splitting samtools view failed (exit $rc7)"
[ $rc7b -eq 0 ] || die "step 7: compressing samtools failed (exit $rc7b)"
NCOM=$(samtools view -c -@ "$T" "$P-committed.bam") || die "cannot count committed"
stat_kv out_committed_records "$NCOM"
[ "$NCOM" -eq $(( NA - NCSEL + NEXTRA )) ] || die "step 7: committed has $NCOM records, expected $(( NA - NCSEL + NEXTRA )) = $NA - $NCSEL + $NEXTRA"
say "step 7: $P-committed.bam has $NCOM records"

cat >&2 <<EOF

[collx] Num. unpaired reads (map U): $NUNP
[collx] Extract $NMATCH reads from BAM
[collx] Num. remaining records in the map = $(( NUNP - NMATCH ))
[collx] Outputs:
  $P-committed.bam        $NCOM records
  $P-deferred.bam         $NDEF records
  $P-deferred-R1.fq.gz    $(( NDEF/2 )) reads
  $P-deferred-R2.fq.gz    $(( NDEF/2 )) reads
EOF
say "done"
exit 0
