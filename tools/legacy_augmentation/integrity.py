"""Independent checks of a hap2-long bundle: base sequences unchanged, context contigs and
liftmap gaps match the selected indels, modifier carriers and hapdp.tsv masks agree."""
import json
from pathlib import Path
import pysam


def check_hap(out,reference):
    out=Path(out);selected={}
    for line in (out/'tmp-insertions.vcf').read_text().splitlines():
        if line.startswith('#'):continue
        f=line.split('\t');selected[f[0]+'_contig_'+f[1]]=(f[0],int(f[1]),f[3],f[4],f[9:])
    templates={name:[row] for name,row in selected.items()}
    for line in (out/'hap-contigs.tsv').read_text().splitlines():
        if line.startswith('#'):continue
        fields=line.split('\t');events=[]
        for value in fields[7].split(','):
            pos,ref,alt=value.split(':');key=fields[1]+'_contig_'+pos
            row=selected[key]
            if row[2:4]!=(ref,alt):raise ValueError('Haplotype metadata disagrees with selected alleles')
            events.append(row)
        carriers=[gt.split('|') for gt in events[0][4]]
        common=sum(all(row[4][sample].split('|')[copy]=='1' for row in events) for sample,gts in enumerate(carriers) for copy in range(len(gts)))
        if common!=int(fields[6]) or not common:raise ValueError('Template has incorrect or absent carriers')
        templates[fields[0]]=events
    side={}
    for line in (out/'aug.liftmap').read_text().splitlines():
        f=line.split('\t')
        if f[0] in side or len(f) not in (4,5):raise ValueError('Invalid sidecar row')
        side[f[0]]=[(int(f[3]),int(f[1]),int(f[2]))]
        if len(f)==5:side[f[0]] += [tuple(map(int,g.split(':'))) for g in f[4].split(';')]
    with pysam.FastaFile(str(reference)) as base,pysam.FastaFile(str(out/'aug.fasta')) as aug:
        for name,length in zip(base.references,base.lengths):
            if aug.get_reference_length(name)!=length:raise ValueError('Hap1 dictionary changed')
            for start in range(0,length,1<<20):
                if base.fetch(name,start,min(length,start+(1<<20)))!=aug.fetch(name,start,min(length,start+(1<<20))):
                    raise ValueError('Hap1 base was modified')
        if set(templates)!=set(side) or set(templates)!=set(aug.references)-set(base.references):
            raise ValueError('Context companion dictionary mismatch')
        for name,events in templates.items():
            chrom=events[0][0];cursor=max(0,events[0][1]-101);sequence='';gaps=[]
            for ch,pos,ref,alt,_ in events:
                if ch!=chrom or pos-1<cursor:raise ValueError('Invalid context event ordering')
                if base.fetch(ch,pos-1,pos-1+len(ref)).upper()!=ref.upper():raise ValueError('REF mismatch')
                sequence+=base.fetch(ch,cursor,pos-1);gaps.append((len(sequence)+1,len(alt)-1,len(ref)-1))
                sequence+=alt;cursor=pos-1+len(ref)
            sequence+=base.fetch(chrom,cursor,min(base.get_reference_length(chrom),cursor+100))
            if aug.fetch(name)!=sequence or side[name]!=gaps:raise ValueError('Context sequence/gap mismatch')
        expected=[];samples=None;ploidy=None
        for line in (out/'mod.syn.vcf').read_text().splitlines():
            if line.startswith('#CHROM'):
                samples=line.split('\t')[9:-1];continue
            if line.startswith('#'):continue
            f=line.split('\t');name,pos,ref,alt=f[0],int(f[1]),f[3],f[4]
            if len(f)!=len(samples)+10 or name not in aug.references or pos<1 or pos>aug.get_reference_length(name):
                raise ValueError('Invalid modifier coordinates or carriers')
            if aug.fetch(name,pos-1,pos).upper()!=ref.upper():raise ValueError('Modifier REF mismatch')
            if any(not all(a in ('0','1') for a in gt.split('|')) for gt in f[9:-1]):raise ValueError('Invalid carrier')
            widths=[len(gt.split('|')) for gt in f[9:-1]]
            if ploidy is None:ploidy=widths
            if ploidy!=widths:raise ValueError('Carrier ploidy changed')
            alleles=[a for gt in f[9:-1] for a in gt.split('|')]
            mask=sum(1<<index for index,allele in enumerate(alleles) if allele=='1')
            synthetic='1' if '_contig_' not in name and not mask else '0'
            if f[-1]!=synthetic:raise ValueError('Incorrect SYN carrier')
            if mask and len(ref)==len(alt)==1 and ref!=alt and alt in 'ACGT':
                expected.append((name,str(pos),alt,f'{mask>>64:016x}',f'{mask&((1<<64)-1):016x}'))
        observed=[tuple(line.split('\t')) for line in (out/'hapdp.tsv').read_text().splitlines() if not line.startswith('#')]
        if expected!=observed:raise ValueError('Hap-DP differs from real source carriers')
