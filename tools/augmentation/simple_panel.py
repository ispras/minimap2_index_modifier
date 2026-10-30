"""Plain BASE1 single-context subprofile, without genotype interpretation."""
import bisect
import collections
import gzip
import json
from pathlib import Path
import re
import shutil
import textwrap


def lines(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as stream:
        yield from stream


def selected_indels(path):
    """Preserve historical first-position occupation before allele-length checks."""
    seen=set()
    for line in lines(path):
        if line.startswith('#'):
            continue
        fields=line.rstrip('\n').split('\t')
        if len(fields)<8:
            raise ValueError('Malformed panel record')
        chrom,pos,ref=fields[0],int(fields[1]),fields[3]
        alts=fields[4].split(',');refs=ref.split(',')
        if any(label in chrom for label in ('random','chrUn','chrM')):
            continue
        insertion=any(len(alt)>len(ref) for alt in alts)
        deletion=any(len(alt)<len(ref) for alt in alts)
        if not insertion and not deletion:
            continue
        key=(chrom,pos)
        if key in seen:
            continue
        info=dict(item.split('=',1) for item in fields[7].split(';') if '=' in item)
        try:af=float(info['AF'].split(',')[0])
        except (KeyError,ValueError):continue
        if not af>0.1:
            continue
        seen.add(key)
        fields[3],fields[4]=refs[0],alts[0]
        if insertion and not 5<=len(fields[4])<=100:
            continue
        if deletion and not 5<=len(fields[3])<=100:
            continue
        yield fields


def build(reference,panel,out):
    import pysam
    fasta=pysam.FastaFile(str(reference))
    stats=collections.Counter();events=list(selected_indels(panel))
    geometry={};positions=collections.defaultdict(list)
    seen_order={}
    for fields in events:
        chrom,pos,ref,alt=fields[0],int(fields[1]),fields[3],fields[4]
        if not re.fullmatch(r'[A-Za-z0-9]+',chrom):
            raise ValueError('Historical chain naming supports only plain chromosome names')
        if chrom not in fasta.references or pos<=101 or pos+len(ref)-1+100>fasta.get_reference_length(chrom):
            raise ValueError('Unsupported historical BASE1 boundary geometry; no silent repair')
        if pos<seen_order.get(chrom,0):
            raise ValueError('Panel must retain chromosome-position order')
        seen_order[chrom]=pos
        if not ref or not alt or len(ref)==len(alt) or set(ref.upper()+alt.upper())-set('ACGT'):
            raise ValueError('Unsupported first-ALT context allele')
        # BASE1 uses a prebuilt substituted reference. Its indel anchor may differ.
        start=pos-101;stop=pos-1+len(ref)+100
        sequence=fasta.fetch(chrom,start,pos-1)+alt+fasta.fetch(chrom,pos-1+len(ref),stop)
        name=f'{chrom}_contig_{pos}'
        shared=fasta.fetch(chrom,pos-1,pos).upper() in 'ACGT' and fasta.fetch(chrom,pos-1,pos).upper()==alt[0].upper()
        first,ins,dele=(101,len(alt)-1,len(ref)-1) if shared else (100,len(alt),len(ref))
        geometry[name]=(chrom,pos,ref,alt,start,stop,sequence,first,ins,dele)
        positions[chrom].append(pos)
    out=Path(out);out.mkdir()
    shutil.copyfile(reference,out/'aug.fasta')
    with (out/'selected-indels.vcf').open('w') as stream:
        for line in lines(panel):
            if line.startswith('#'):stream.write(line)
        for fields in events:stream.write('\t'.join(fields)+'\n')
    with (out/'aug.fasta').open('a') as augmented,(out/'aug.chain').open('w') as chain, \
         (out/'aug.bed').open('w') as bed,(out/'aug.liftmap').open('w') as lift, \
         (out/'contigs.jsonl').open('w') as metadata:
        number=1
        for name,(chrom,pos,ref,alt,start,stop,sequence,first,ins,dele) in geometry.items():
            augmented.write('>'+name+'\n'+textwrap.fill(sequence,70)+'\n')
            chain.write(f'chain {len(sequence)} {name} {len(sequence)} + 0 {len(sequence)} {chrom} {fasta.get_reference_length(chrom)} + {start} {stop} {number}\n{first} {ins}\t{dele}\n100\n\n')
            bed.write(f'{chrom}\t{pos-10}\t{pos+10}\n')
            # Lift-map geometry, distinct from the chain anchor test.
            lift.write(f'{name}\t{len(alt)-1}\t{len(ref)-1}\t101\n')
            metadata.write(json.dumps(dict(name=name,chrom=chrom,start=start,end=stop,length=len(sequence),
                events=[dict(chrom=chrom,pos=pos,ref=ref,alt=alt)],origins=['population-panel'],
                carrier_semantics='GT opaque, not used; plain index modifier',gaps=[[101,len(alt)-1,len(ref)-1]]))+'\n')
            number+=1
        for chrom,length in zip(fasta.references,fasta.lengths):
            chain.write(f'chain {length} {chrom} {length} + 0 {length} {chrom} {length} + 0 {length} {number}\n{length}\n\n')
            number+=1
    with (out/'modifier.vcf').open('w') as modifier:
        for line in lines(panel):
            if line.startswith('#'):
                modifier.write(line);continue
            fields=line.rstrip('\n').split('\t')
            chrom,pos,ref,alt=fields[0],int(fields[1]),fields[3],fields[4]
            if chrom not in fasta.references or pos<1 or pos+len(ref)-1>fasta.get_reference_length(chrom):
                raise ValueError('Invalid panel coordinate/reference')
            snv=len(ref)==len(alt)==1
            reference_base=fasta.fetch(chrom,pos-1,pos).upper()
            ref_ok=not snv or reference_base not in 'ACGT' or reference_base==ref.upper()
            if snv:
                if ref_ok:modifier.write(line);stats['base_snvs']+=1
                else:stats['ref_mismatch_base']+=1
            loci=positions.get(chrom,[])
            if pos in loci:
                continue
            if not snv:
                continue
            for anchor in loci[bisect.bisect_left(loci,pos-100):bisect.bisect_right(loci,pos+100)]:
                name=f'{chrom}_contig_{anchor}';g=geometry[name]
                distance=pos-anchor
                if distance<0:offset=100+distance
                elif distance>=len(g[2]):offset=100+len(g[3])+distance-len(g[2])
                else:stats['inside_deleted_segment']+=1;continue
                if not 0<=offset<len(g[6]):stats['out_of_bounds']+=1;continue
                if not ref_ok:stats['ref_mismatch_contig']+=1;continue
                projected=fields.copy();projected[0]=name;projected[1]=str(offset+1)
                modifier.write('\t'.join(projected)+'\n');stats['contig_snvs']+=1
    fasta.close()
    (out/'panel.stats.json').write_text(json.dumps(dict(stats),sort_keys=True,indent=2)+'\n')
