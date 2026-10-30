"""Usage:
  hapdp.py VCF FAI OUT [--min-af X] [--stats FILE]
      SNV ALT table for base sequences in FAI: chrom, pos, merged ALT bases.
  hapdp.py --hapdp VCF FAI OUT [--no-contigs] [--stats FILE]
      Per-ALT haplotype carrier masks (128 bits as two hex words, high then low);
      the SYN sample is ignored and at most 128 haplotypes are allowed.
"""
import gzip
import shutil
import subprocess
import sys

def records(vcf, with_gt):
    if shutil.which('bcftools'):
        fmt = '%CHROM\t%POS\t%REF\t%ALT' + ('[\t%GT]' if with_gt else '') + '\n'
        cmd = ['bcftools', 'query', '-f', fmt]
        if with_gt:
            samples = subprocess.run(['bcftools', 'query', '-l', vcf], capture_output=True, text=True).stdout.split()
            if 'SYN' in samples:
                cmd += ['-s', '^SYN']
        p = subprocess.Popen(cmd + [vcf], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            t = line.rstrip('\n').split('\t')
            yield (t[0], t[1], t[2], t[3], t[4:])
        if p.wait() != 0:
            sys.exit(f'ERROR: bcftools query {vcf} exited with {p.returncode}')
    else:
        op = gzip.open if vcf.endswith('.gz') else open
        keep = None
        with op(vcf, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    keep = [i for i, s in enumerate(line.rstrip('\n').split('\t')) if i >= 9 and s != 'SYN']
                if line.startswith('#'):
                    continue
                t = line.rstrip('\n').split('\t')
                yield (t[0], t[1], t[3], t[4], [t[i].split(':')[0] for i in keep] if with_gt else [])

def gt_freq(gts, k):
    n = a = 0
    for g in gts:
        for x in g.replace('|', '/').split('/'):
            if x == '.' or x == '':
                continue
            n += 1
            a += x == k
    return a / n if n else 0.0

def samples_of(vcf):
    if shutil.which('bcftools'):
        return [x for x in subprocess.run(['bcftools', 'query', '-l', vcf], capture_output=True, text=True, check=True).stdout.split() if x != 'SYN']
    op = gzip.open if vcf.endswith('.gz') else open
    with op(vcf, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                return [x for x in line.rstrip('\n').split('\t')[9:] if x != 'SYN']
    sys.exit('ERROR: no #CHROM line')

def main_hapdp(vcf, fai, out, stats_fn, contigs):
    base = set()
    with open(fai) as f:
        for line in f:
            base.add(line.split('\t', 1)[0])
    names = samples_of(vcf)
    n = dict(rec=0, base=0, contig=0, contig_skipped=0, not_snv_ref=0, alt_not_snv=0, no_carrier=0, lines_out=0, alleles_carried=0)
    ploidy = None
    op = gzip.open if out.endswith('.gz') else open
    body = []
    for c, pos, ref, alt, gts in records(vcf, True):
        n['rec'] += 1
        is_ctg = c not in base
        if is_ctg:
            n['contig'] += 1
            if not contigs:
                n['contig_skipped'] += 1
                continue
        else:
            n['base'] += 1
        if len(gts) != len(names):
            sys.exit(f'ERROR: {c}:{pos}: {len(gts)} GT fields for {len(names)} samples')
        gl = [g.replace('/', '|').split('|') for g in gts]
        pl = [len(x) for x in gl]
        if ploidy is None:
            ploidy = pl
        elif pl != ploidy:
            sys.exit(f'ERROR: {c}:{pos}: sample ploidy differs from the first record')
        if len(ref) != 1 or ref.upper() not in 'ACGT':
            n['not_snv_ref'] += 1
            continue
        flat = [x for g in gl for x in g]
        for j, x in enumerate(alt.split(',')):
            x = x.upper()
            if not (len(x) == 1 and x in 'ACGT' and (x != ref.upper())):
                n['alt_not_snv'] += 1
                continue
            k = str(j + 1)
            m = 0
            for h, al in enumerate(flat):
                if al == k:
                    m |= 1 << h
            if m == 0:
                n['no_carrier'] += 1
                continue
            n['lines_out'] += 1
            n['alleles_carried'] += bin(m).count('1')
            body.append(f'{c}\t{pos}\t{x}\t{m >> 64:016x}\t{m & (1 << 64) - 1:016x}\n')
    hap = [f'{s}:{p}' for s, k in zip(names, ploidy or []) for p in range(k)]
    if len(hap) > 128:
        sys.exit(f'ERROR: {len(hap)} haplotypes > 128')
    with op(out, 'wt') as fo:
        fo.write(f"#hapdp v1: mkvaralt.py --hapdp {vcf} (contigs {('yes' if contigs else 'no')})\n#n_hap\t{len(hap)}\n")
        for h, nm in enumerate(hap):
            fo.write(f'#hap\t{h}\t{nm}\n')
        fo.writelines(body)
    msg = 'mkvaralt --hapdp: records {rec}, on base {base}, on contigs {contig} (skipped {contig_skipped}); non-SNV REF {not_snv_ref}, non-SNV ALT {alt_not_snv}, ALT without carrier (SYN) {no_carrier}; lines {lines_out}, carried alleles {alleles_carried}'.format(**n) + f'; haplotypes {len(hap)}'
    print(msg, file=sys.stderr)
    if stats_fn:
        with open(stats_fn, 'w') as f:
            f.write(msg + '\n')
            for k, v in n.items():
                f.write(f'{k}\t{v}\n')

def main():
    a = sys.argv[1:]
    stats_fn, min_af = (None, None)
    hapdp = '--hapdp' in a
    contigs = '--no-contigs' not in a
    a = [x for x in a if x not in ('--hapdp', '--no-contigs')]
    if '--stats' in a:
        i = a.index('--stats')
        stats_fn = a[i + 1]
        del a[i:i + 2]
    if '--min-af' in a:
        i = a.index('--min-af')
        min_af = float(a[i + 1])
        del a[i:i + 2]
    if len(a) != 3:
        sys.exit(__doc__)
    vcf, fai, out = a
    if hapdp:
        if min_af is not None:
            sys.exit('ERROR: --min-af cannot be combined with --hapdp')
        return main_hapdp(vcf, fai, out, stats_fn, contigs)
    base = set()
    with open(fai) as f:
        for line in f:
            base.add(line.split('\t', 1)[0])
    n = dict(rec=0, base=0, other=0, not_snv_ref=0, alt_not_snv=0, alt_low_af=0, pos_out=0, alleles_out=0, merged=0)
    op = gzip.open if out.endswith('.gz') else open
    pend_key, pend = (None, set())

    def emit(fo):
        if pend_key is None:
            return
        alts = ''.join((x for x in 'ACGT' if x in pend))
        fo.write(f'{pend_key[0]}\t{pend_key[1]}\t{alts}\n')
        n['pos_out'] += 1
        n['alleles_out'] += len(alts)
    with op(out, 'wt') as fo:
        for c, pos, ref, alt, gts in records(vcf, min_af is not None):
            n['rec'] += 1
            if c not in base:
                n['other'] += 1
                continue
            n['base'] += 1
            if len(ref) != 1 or ref.upper() not in 'ACGT':
                n['not_snv_ref'] += 1
                continue
            got = set()
            for j, x in enumerate(alt.split(',')):
                x = x.upper()
                if len(x) == 1 and x in 'ACGT' and (x != ref.upper()):
                    if min_af is not None and gt_freq(gts, str(j + 1)) <= min_af:
                        n['alt_low_af'] += 1
                        continue
                    got.add(x)
                else:
                    n['alt_not_snv'] += 1
            if not got:
                continue
            key = (c, pos)
            if key == pend_key:
                n['merged'] += 1
                pend |= got
            else:
                emit(fo)
                pend_key, pend = (key, got)
        emit(fo)
    msg = 'mkvaralt: records {rec}, on base {base}, on other sequences (contigs) {other}; on base: non-SNV REF {not_snv_ref}, non-SNV ALT {alt_not_snv}, ALT with frequency <= threshold {alt_low_af}; positions {pos_out}, alleles {alleles_out}, merged repeated positions {merged}'.format(**n)
    if min_af is not None:
        msg += f'; --min-af {min_af}'
    print(msg, file=sys.stderr)
    if stats_fn:
        with open(stats_fn, 'w') as f:
            f.write(msg + '\n')
            for k, v in n.items():
                f.write(f'{k}\t{v}\n')
if __name__ == '__main__':
    main()
