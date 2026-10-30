#!/usr/bin/env python3
"""Build the augmented reference: indel contigs, liftmap, hap-DP table and the SNV VCF for the index."""
import argparse
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def check_reference(reference):
    import pysam
    path=Path(reference)
    if not path.is_file() or not Path(str(path)+'.fai').is_file():
        raise ValueError('Reference and existing FAI companion are required')
    with pysam.FastaFile(str(path)) as fasta:
        if not fasta.references or len(set(fasta.references))!=len(fasta.references):
            raise ValueError('Invalid reference dictionary')
        return dict(zip(fasta.references,fasta.lengths))


REQUIRED={'aug.fasta','aug.fasta.fai','aug.liftmap','contigs.jsonl','modifier.vcf','modifier.vcf.gz',
          'modifier.vcf.gz.tbi','panel.stats.json'}
REQUIRED_BY_MODE={'direct-grch38':{'hapdp.tsv','coverage.json'},'base1-plain':{'aug.chain','aug.bed','selected-indels.vcf'}}


def main():
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest='mode',required=True)
    for mode in ('direct-grch38','base1-plain'):
        command=sub.add_parser(mode)
        for key in ('ref','out'):command.add_argument('--'+key,required=True)
        command.add_argument('--bcftools',default='bcftools')
        if mode=='direct-grch38':
            for key in ('h1','h2','paf1','paf2'):command.add_argument('--'+key,required=True)
        else:
            command.add_argument('--panel',required=True)
    args=parser.parse_args()
    import pysam
    out=Path(args.out).absolute()
    if out.exists():
        raise ValueError('Output directory already exists: '+str(out))
    reference=Path(args.ref).resolve()
    chroms=check_reference(reference)
    if args.mode=='direct-grch38':
        from paf_selection import read_paf
        for label in ('paf1','paf2'):
            for chrom,length in chroms.items():
                for record in read_paf(getattr(args,label).replace('{chrom}',chrom)):
                    if record['tn']!=chrom or record['tl']!=length or not 0<=record['ts']<record['te']<=length or not 0<=record['qs']<record['qe']<=record['ql']:
                        raise ValueError('PAF dictionary/coordinate mismatch: '+label+':'+chrom)
    out.parent.mkdir(parents=True,exist_ok=True)
    if args.mode=='direct-grch38':
        command=[sys.executable,'-B',str(HERE/'direct_core.py'),'--ref',str(reference),
                 '--h1',str(Path(args.h1).resolve()),'--h2',str(Path(args.h2).resolve()),
                 '--paf1',str(Path(args.paf1).absolute()),'--paf2',str(Path(args.paf2).absolute()),
                 '--pafsel',str(HERE/'paf_selection.py'),'--out',str(out)]
        result=subprocess.run(command,capture_output=True,text=True)
        if out.is_dir():
            (out/'construction.stdout').write_text(result.stdout)
            (out/'construction.stderr').write_text(result.stderr)
        if result.returncode:raise RuntimeError('Direct construction failed: '+result.stderr[-2000:])
    else:
        from simple_panel import build
        build(reference,Path(args.panel).resolve(),out)
    pysam.faidx(str(out/'aug.fasta'))
    for label,command in [
        ('reheader',[args.bcftools,'reheader','-f',str(out/'aug.fasta.fai'),'-o',str(out/'modifier.headers.vcf'),str(out/'modifier.vcf')]),
        ('sort',[args.bcftools,'sort','-m','512M','-T',str(out/'sort-tmp'),'-Ov','-o',str(out/'modifier.sorted.vcf'),str(out/'modifier.headers.vcf')])]:
        result=subprocess.run(command,capture_output=True,text=True)
        (out/(label+'.stdout')).write_text(result.stdout);(out/(label+'.stderr')).write_text(result.stderr)
        if result.returncode:raise RuntimeError(label+' failed')
    pysam.tabix_compress(str(out/'modifier.sorted.vcf'),str(out/'modifier.vcf.gz'),force=False)
    pysam.tabix_index(str(out/'modifier.vcf.gz'),preset='vcf',force=False)
    missing=(REQUIRED|REQUIRED_BY_MODE[args.mode])-{p.name for p in out.iterdir()}
    if missing:raise RuntimeError('Missing outputs: '+', '.join(sorted(missing)))
    from integrity import check
    check(out,reference,args.mode)
    print('AUGMENTATION_BUILD_COMPLETE',args.mode,out)


if __name__=='__main__':main()
