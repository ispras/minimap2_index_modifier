"""Write the context-contig liftmap: liftmap.py OUT GDIR [FLANK=100].

One row per *_contig_* in GDIR/aug.fasta.fai: name, inserted length, deleted length and the
1-based junction offset; haplotype contigs add ";"-joined "junction:ins:del" for later indels.
Exits non-zero if contig lengths or geometry disagree with the source VCF/TSV."""
import sys, collections
out = sys.argv[1]
G = sys.argv[2]
FLANK = int(sys.argv[3]) if len(sys.argv) > 3 else 100
chrom_lengths = {}
for line in open(G + '/base.fasta.fai'):
    fields = line.split('\t')
    chrom_lengths[fields[0]] = int(fields[1])

def right_flank(chrom, pos, ref):
    remaining = chrom_lengths[chrom] - (pos - 1 + len(ref))
    if remaining < 0:
        raise ValueError('variant REF extends past chromosome end: %s:%d' % (chrom, pos))
    return min(FLANK, remaining)
fai = collections.OrderedDict()
dupname = 0
for l in open(G + '/aug.fasta.fai'):
    t = l.split('\t')
    if '_contig_' in t[0]:
        if t[0] in fai:
            dupname += 1
        fai[t[0]] = int(t[1])
recs = {}
dup = 0
for l in open(G + '/tmp-insertions.vcf'):
    if l[0] == '#':
        continue
    t = l.split('\t', 5)
    name = '%s_contig_%s' % (t[0], t[1])
    if name in recs:
        dup += 1
    recs[name] = (int(t[1]), t[3], t[4])
hap = collections.OrderedDict()
try:
    fh = open(G + '/hap-contigs.tsv')
except FileNotFoundError:
    fh = []
for l in fh:
    if l[0] == '#':
        continue
    t = l.rstrip('\n').split('\t')
    name, ch, p1, w, L = (t[0], t[1], int(t[2]), int(t[3]), int(t[4]))
    ind = []
    for x in t[7].split(','):
        p, ref, alt = x.split(':')
        ind.append((int(p), ref, alt))
    if name in hap or name in recs:
        dup += 1
    hap[name] = (ch, p1, w, L, ind)
missing = [n for n in fai if n not in recs and n not in hap]
extra = [n for n in list(recs) + list(hap) if n not in fai]
bad_len = []
kinds = collections.Counter()
anchor_bad = 0
bad_geom = 0
with open(out, 'w') as o:
    for n, L in fai.items():
        if n in recs:
            pos, ref, alt = recs[n]
            ins, dl = (len(alt) - 1, len(ref) - 1)
            jct = min(FLANK, pos - 1) + 1
            chrom = n.split('_contig_', 1)[0]
            exp = jct - 1 + len(alt) + right_flank(chrom, pos, ref)
            if exp != L:
                bad_len.append((n, L, exp))
            if alt[0] != ref[0]:
                anchor_bad += 1
            kinds['complex' if ins > 0 and dl > 0 else 'ins' if ins > 0 else 'del' if dl > 0 else 'neither'] += 1
            o.write('%s\t%d\t%d\t%d\n' % (n, ins, dl, jct))
        elif n in hap:
            ch, p1, w, Lh, ind = hap[n]
            if w != min(FLANK, p1 - 1) or ind[0][0] != p1 or Lh != L:
                bad_geom += 1
            gaps = []
            off = w
            for t, (p, ref, alt) in enumerate(ind):
                if t > 0:
                    pp, pref, palt = ind[t - 1]
                    between = p - pp - len(pref)
                    if between < 0:
                        bad_geom += 1
                    off += len(palt) + between
                gaps.append((off + 1, len(alt) - 1, len(ref) - 1))
                if alt[0] != ref[0]:
                    anchor_bad += 1
                if len(alt) - 1 == 0 and len(ref) - 1 == 0:
                    kinds['neither'] += 1
            last_pos, last_ref, last_alt = ind[-1]
            exp = off + len(last_alt) + right_flank(ch, last_pos, last_ref)
            if exp != L:
                bad_len.append((n, L, exp))
            kinds['hap%d' % len(ind)] += 1
            j1, i1, d1 = gaps[0]
            o.write('%s\t%d\t%d\t%d\t%s\n' % (n, i1, d1, j1, ';'.join(('%d:%d:%d' % g for g in gaps[1:]))))
print('contigs in .fai: %d (single %d, haplotype %d); panel records: %d (duplicate names: %d, in .fai: %d); without description: %d; descriptions without contig: %d' % (len(fai), sum((1 for n in fai if n in recs)), sum((1 for n in fai if n in hap)), len(recs), dup, dupname, len(missing), len(extra)))
print('kinds:', dict(sorted(kinds.items())), '; anchor ALT[0]!=REF[0]:', anchor_bad, '; length mismatches:', len(bad_len), '; geometry violations:', bad_geom)
for x in bad_len[:5]:
    print('  length mismatch:', x)
if missing or extra or bad_len or kinds['neither'] or dup or dupname or bad_geom:
    print('LIFTMAP CHECK FAILED')
    sys.exit(1)
print('OK', out, 'flank', FLANK)
