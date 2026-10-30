"""Independent output geometry and carrier checks for packaged profiles."""
import itertools
import json
from pathlib import Path
import pysam


def check(out,reference,mode):
    out=Path(out)
    with pysam.FastaFile(str(reference)) as base,pysam.FastaFile(str(out/'aug.fasta')) as aug:
        for name,length in zip(base.references,base.lengths):
            if aug.get_reference_length(name)!=length:raise ValueError('Base dictionary changed')
            for start in range(0,length,1<<20):
                if base.fetch(name,start,min(length,start+(1<<20)))!=aug.fetch(name,start,min(length,start+(1<<20))):
                    raise ValueError('Base reference sequence changed')
        side={}
        for line in (out/'aug.liftmap').read_text().splitlines():
            fields=line.split('\t')
            if len(fields) not in (4,5) or fields[0] in side:raise ValueError('Invalid liftmap')
            gaps=[(int(fields[3]),int(fields[1]),int(fields[2]))]
            if len(fields)==5:gaps += [tuple(map(int,item.split(':'))) for item in fields[4].split(';')]
            side[fields[0]]=gaps
        names=set()
        with (out/'contigs.jsonl').open() as stream:
            for line in stream:
                item=json.loads(line);name=item['name']
                if name in names:raise ValueError('Duplicate metadata')
                names.add(name);cursor=item['start'];sequence='';gaps=[]
                if not 0<=cursor<item['end']<=base.get_reference_length(item['chrom']):
                    raise ValueError('Metadata span outside reference')
                if mode=='direct-grch38' and (not item['origins'] or set(item['origins'])-{'hap1','hap2'}):
                    raise ValueError('Invalid source carriers')
                for event in item['events']:
                    pos=event['pos']-1
                    if pos<cursor:raise ValueError('Overlapping events')
                    if mode=='direct-grch38' and base.fetch(item['chrom'],pos,pos+len(event['ref'])).upper()!=event['ref']:
                        raise ValueError('Event REF mismatch')
                    sequence+=base.fetch(item['chrom'],cursor,pos)
                    gaps.append((len(sequence)+1,len(event['alt'])-1,len(event['ref'])-1))
                    sequence+=event['alt'];cursor=pos+len(event['ref'])
                sequence+=base.fetch(item['chrom'],cursor,item['end'])
                observed=aug.fetch(name)
                if sequence.upper()!=observed.upper() or len(sequence)!=item['length']:
                    raise ValueError('Contig sequence geometry mismatch')
                if gaps!=side.get(name) or [list(g) for g in gaps]!=item['gaps']:
                    raise ValueError('Gap companion mismatch')
        if names!=set(side) or names!=set(aug.references)-set(base.references):
            raise ValueError('Companion dictionary mismatch')
        with (out/'modifier.vcf').open() as stream:
            records=(line.rstrip('\n').split('\t') for line in stream if not line.startswith('#'))
            hap_stream=(out/'hapdp.tsv').open() if mode=='direct-grch38' else None
            hap_records=(line.rstrip('\n').split('\t') for line in hap_stream if not line.startswith('#')) if hap_stream else itertools.repeat(None)
            try:
                for fields,hap in itertools.zip_longest(records,hap_records) if hap_stream else zip(records,hap_records):
                    if fields is None or (hap_stream and hap is None):raise ValueError('Carrier table count mismatch')
                    name,pos,ref,alt=fields[0],int(fields[1]),fields[3],fields[4]
                    if name not in aug.references or pos<1 or pos+len(ref)-1>aug.get_reference_length(name):
                        raise ValueError('Modifier outside sequence')
                    actual=aug.fetch(name,pos-1,pos-1+len(ref)).upper()
                    if actual!=ref.upper() and (mode=='direct-grch38' or any(a in 'ACGT' and a!=b for a,b in zip(actual,ref.upper()))):
                        raise ValueError('Modifier REF mismatch')
                    if hap_stream:
                        if len(fields)!=11 or any(gt not in ('0','1') for gt in fields[9:]):
                            raise ValueError('Invalid direct source-path carriers')
                        mask=int(fields[9])|(int(fields[10])<<1)
                        if not mask or hap[:3]!=[name,str(pos),alt] or int(hap[3],16)!=0 or int(hap[4],16)!=mask:
                            raise ValueError('VCF/hap-DP carrier mismatch')
            finally:
                if hap_stream:hap_stream.close()
