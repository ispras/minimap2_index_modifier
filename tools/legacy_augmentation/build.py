#!/usr/bin/env python3
"""Build a hap1 augmented reference (hap2-long: indel contexts, chain, modifier VCF with SYN,
liftmap, hap-DP) or a panel-consensus SNV reference (panel-reference)."""
import argparse
import gzip
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
PARAMETERS={'AF_INDEL':'0.1','FIXDEDUP':'1','HAP':'1','HAP_GAP':'100','HAP_MAXSPAN':'400','HAP_MAXT':'16','HAP_MAXGAPS':'8','MIN_INDEL_ALLELE':'2'}


def rows(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as stream:yield from stream


def inspect_inputs(reference,panel,profile):
    import pysam
    if not reference.is_file() or not Path(str(reference)+'.fai').is_file() or not panel.is_file():
        raise ValueError('Reference, existing FAI and panel are required')
    names=None;ploidy=None;last={};records=0
    with pysam.FastaFile(str(reference)) as fasta:
        for line in rows(panel):
            if line.startswith('#CHROM'):
                names=line.rstrip('\n').split('\t')[9:]
                if not names or len(set(names))!=len(names) or 'SYN' in names:raise ValueError('Invalid source sample names')
                if profile=='base1-44-all-alt' and len(names)!=45:raise ValueError('Positional selector requires exactly 45 source samples')
            if line.startswith('#'):continue
            f=line.rstrip('\n').split('\t')
            if names is None or len(f)!=9+len(names):raise ValueError('Malformed VCF sample columns')
            chrom,pos,ref,alt=f[0],int(f[1]),f[3],f[4]
            if profile=='base1-44-all-alt' and f[8].split(':')[0]!='GT':
                raise ValueError('Positional selector requires GT as the first sample subfield')
            if chrom not in fasta.references or pos<1 or pos+len(ref)-1>fasta.get_reference_length(chrom):
                raise ValueError('Variant coordinate outside reference')
            if pos<last.get(chrom,0):raise ValueError('Unsorted source positions')
            last[chrom]=pos;records+=1
            if profile.startswith('h2l-'):
                if f[8]!='GT' or not re.fullmatch(r'[A-Za-z0-9]+',chrom) or ',' in alt:
                    raise ValueError('Supported legacy path requires biallelic GT-only records and plain chromosome names')
                if set(ref.upper()+alt.upper())-set('ACGT') or fasta.fetch(chrom,pos-1,pos-1+len(ref)).upper()!=ref.upper():
                    raise ValueError('Hap1 REF mismatch or unsupported allele')
                gts=f[9:]
                if any(not re.fullmatch(r'[01](?:\|[01])?',gt) for gt in gts):
                    raise ValueError('Missing, unphased or unsupported carrier')
                widths=[len(gt.split('|')) for gt in gts]
                if ploidy is None:ploidy=widths
                if widths!=ploidy:raise ValueError('Variable ploidy')
                if sum(widths)>127:raise ValueError('No room for the synthetic SYN carrier in the 128-haplotype index')
                if len(ref)!=len(alt) and pos<=101:
                    raise ValueError('Indel too close to the chromosome start for context geometry')
    if names is None or not records:raise ValueError('Empty or incomplete source panel')
    return {'samples':names,'records':records,'ploidy':ploidy}


def synthesize(modifier,destination):
    count=0
    with open(modifier) as source,open(destination,'x') as out:
        for line in source:
            if line.startswith('##'):out.write(line);continue
            if line.startswith('#CHROM'):out.write(line.rstrip('\n')+'\tSYN\n');continue
            f=line.rstrip('\n').split('\t')
            value='0'
            if '_contig_' not in f[0] and '1' not in ''.join(f[9:]):value='1';count+=1
            out.write('\t'.join(f+[value])+'\n')
    return count


REQUIRED={'hap2-long':{'aug.fasta','aug.fasta.fai','aug.liftmap','hapdp.tsv','mod.syn.vcf.gz','mod.syn.vcf.gz.tbi','hap-contigs.tsv','tmp-insertions.vcf','aug.chain'},
          'panel-reference':{'base.fasta','base.fasta.fai','selected.vcf.gz','selected.vcf.gz.tbi','leftover.vcf.gz','leftover.vcf.gz.tbi'}}


def main():
    parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest='mode',required=True)
    for mode in ('hap2-long','panel-reference'):
        p=sub.add_parser(mode)
        for key in ('ref','panel','out'):p.add_argument('--'+key,required=True)
        p.add_argument('--profile',required=True,choices=('h2l-5k','h2l-10k') if mode=='hap2-long' else ('base1-44-all-alt','later-af075'))
        p.add_argument('--bcftools',default='bcftools')
    args=parser.parse_args()
    import pysam
    out=Path(args.out).absolute();reference=Path(args.ref).resolve();panel=Path(args.panel).resolve()
    if out.exists():raise ValueError('Output directory already exists: '+str(out))
    inspect_inputs(reference,panel,args.profile)
    out.mkdir(parents=True)
    environment=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    def run(label,argv,stdout=None):
        stdout=out/(label+'.stdout') if stdout is None else stdout
        with stdout.open('xb') as so,(out/(label+'.stderr')).open('xb') as se:
            result=subprocess.run(list(map(str,argv)),cwd=out,env=environment,stdout=so,stderr=se)
        if result.returncode:raise RuntimeError(label+' failed')
    def indexed(source,destination):
        pysam.tabix_compress(str(source),str(destination),force=False)
        pysam.tabix_index(str(destination),preset='vcf',force=False)
    if args.mode=='hap2-long':
        environment.update(PARAMETERS,MAX_INDEL_ALLELE='5001' if args.profile=='h2l-5k' else '10001',PYTHONOPTIMIZE='0')
        shutil.copyfile(str(reference)+'.fai',out/'base.fasta.fai')
        run('contexts',[sys.executable,'-B',HERE/'hap_contexts.py',reference,panel,'aug.fasta','aug.chain',str(reference)+'.fai','aug.bed'])
        pysam.faidx(str(out/'aug.fasta'))
        run('reheader',[args.bcftools,'reheader','-f','aug.fasta.fai','-o','reheadered.vcf','vcf-for-indexing.vcf'])
        run('sort',[args.bcftools,'sort','-m','512M','-T','sort-tmp','-Ov','-o','modifier.vcf','reheadered.vcf'])
        synthesize(out/'modifier.vcf',out/'mod.syn.vcf')
        indexed(out/'mod.syn.vcf',out/'mod.syn.vcf.gz')
        run('liftmap',[sys.executable,'-B',HERE/'liftmap.py','aug.liftmap',out,'100'])
        run('hapdp',[sys.executable,'-B',HERE/'hapdp.py','mod.syn.vcf.gz','base.fasta.fai','hapdp.tsv','--hapdp','--stats','hapdp.stats.txt'])
        from integrity import check_hap
        check_hap(out,reference)
    else:
        with (out/'input.vcf').open('x') as stream:
            for line in rows(panel):stream.write(line)
        indexed(out/'input.vcf',out/'panel.vcf.gz')
        run('select',['awk','-F','\t','-f',HERE/(args.profile+'.awk'),'input.vcf'])
        run('sort',[args.bcftools,'sort','-m','512M','-T','sort-tmp','-Oz','-o','sorted.vcf.gz','test.vcf'])
        run('clean',[args.bcftools,'view','-i','REF~"^[ACGT]$"','sorted.vcf.gz','-Oz','-o','selected.vcf.gz'])
        pysam.tabix_index(str(out/'selected.vcf.gz'),preset='vcf',force=False)
        run('consensus',[args.bcftools,'consensus','-f',reference,'selected.vcf.gz'],out/'base.fasta')
        run('isec',[args.bcftools,'isec','-C','panel.vcf.gz','selected.vcf.gz','-w1','-Oz','-o','leftover.vcf.gz'])
        pysam.tabix_index(str(out/'leftover.vcf.gz'),preset='vcf',force=False)
        pysam.faidx(str(out/'base.fasta'))
        with pysam.FastaFile(str(reference)) as original,pysam.FastaFile(str(out/'base.fasta')) as new:
            if original.references!=new.references or original.lengths!=new.lengths:raise ValueError('SNV reference construction changed dictionary/length')
    missing=REQUIRED[args.mode]-{p.name for p in out.iterdir()}
    if missing:raise RuntimeError('Missing outputs: '+', '.join(sorted(missing)))
    print('LEGACY_AUGMENTATION_COMPLETE',args.mode,args.profile,out)


if __name__=='__main__':main()
