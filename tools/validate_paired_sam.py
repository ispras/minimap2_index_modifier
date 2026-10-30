#!/usr/bin/env python3
"""Validate raw paired SAM before any parser can coerce invalid references."""
import argparse
import gzip
import itertools
import json
import re
import sys
from pathlib import Path

# Validation uses assertions; never run with those checks disabled.
if sys.flags.optimize:
    raise RuntimeError('Paired SAM validation requires Python without -O or PYTHONOPTIMIZE')

def reverse_complement(s):return s.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]
def canonical_name(s):
    s=s.split()[0].lstrip('@')
    return s[:-2] if s.endswith(('/1','/2')) else s

def fastq(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as fh:
        while True:
            rows=[fh.readline().rstrip('\n') for _ in range(4)]
            if not rows[0]:
                assert not any(rows),'Truncated FASTQ'
                return
            assert rows[0].startswith('@') and rows[2].startswith('+')
            assert rows[1] and len(rows[1])==len(rows[3])
            yield rows

def reference_end(row):
    """Return the exclusive zero-based reference end, ignoring clipping."""
    return int(row[3])-1+sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',row[5]) if op in 'MDN=X')

def native_tlen(row,mate):
    """Match minimap2 format.c: signed inclusive distance between mapped 5' ends."""
    flag,mflag=int(row[1]),int(mate[1])
    if flag&4 or mflag&4 or row[2]!=mate[2]:
        return 0
    own=reference_end(row)-1 if flag&16 else int(row[3])-1
    other=reference_end(mate)-1 if mflag&16 else int(mate[3])-1
    distance=other-own
    return distance+1 if distance>0 else distance-1 if distance<0 else 0

def validate(sam,r1,r2,fai,expected_pairs=None,output=None):
    with open(fai) as source:
        refs=dict((f[0],int(f[1])) for f in (l.split() for l in source))
    seen_refs={}; count=0; records=0; unprojectable=0; primary=0
    with open(sam) as fh:
        line=fh.readline();header=[]
        while line.startswith('@'):
            header.append(line)
            if line.startswith('@SQ\t'):
                fields=dict(x.split(':',1) for x in line.rstrip().split('\t')[1:])
                assert fields['SN'] not in seen_refs,'Duplicate SQ'
                seen_refs[fields['SN']]=int(fields['LN'])
            line=fh.readline()
        assert seen_refs==refs,'Reference dictionary mismatch'
        if output is not None:
            output.write(''.join(header))
        lines=itertools.chain([line] if line else [],fh)
        groups=itertools.groupby((x.rstrip('\n').split('\t') for x in lines),key=lambda f:f[0])
        left,right=fastq(r1),fastq(r2)
        for name,group in groups:
            a=next(left,None); b=next(right,None)
            assert a is not None and b is not None,'Extra SAM group'
            assert name==canonical_name(a[0])==canonical_name(b[0]),('QNAME order/content',name)
            rows=list(group); chosen={}
            for f in rows:
                assert len(f)>=11
                flag=int(f[1]); records+=1
                assert not flag&1024,'Mapper emitted duplicate flag before optional marking'
                assert f[2]=='*' or f[2] in refs,('Unknown RNAME',name,f[2])
                assert f[6] in ('*','=') or f[6] in refs,('Unknown RNEXT',name,f[6])
                tags={x[:2]:x[5:] for x in f[11:]}
                if 'SA' in tags:assert all(x.split(',')[0] in refs for x in tags['SA'].strip(';').split(';'))
                if not flag&4:
                    assert f[2]!='*' and int(f[3])>0 and f[5]!='*'
                    ops=re.findall(r'(\d+)([MIDNSHP=X])',f[5])
                    assert ''.join(n+op for n,op in ops)==f[5]
                    qspan=sum(int(n) for n,op in ops if op in 'MIS=X')
                    tspan=sum(int(n) for n,op in ops if op in 'MDN=X')
                    assert int(f[3])-1+tspan<=refs[f[2]]
                    if f[9]!='*':assert qspan==len(f[9])
                if not flag&2304:
                    assert flag&1 and bool(flag&64)!=bool(flag&128)
                    end=1 if flag&64 else 2
                    assert end not in chosen,'Duplicated primary end'
                    chosen[end]=(f,flag,tags)
                    inp=a if end==1 else b; seq,qual=inp[1],inp[3]
                    if flag&16:seq,qual=reverse_complement(seq),qual[::-1]
                    assert (f[9],f[10])==(seq,qual),('SEQ/QUAL',name,end)
                    primary+=1
                if 'lu' in tags:
                    unprojectable+=1
                    assert flag&4 and not flag&18
                    assert f[2:6]==['*','0','0','*'] and f[8]=='0'
                    assert tags['lu']=='unprojectable_contig' and '_contig_' in tags.get('lo','')
            assert set(chosen)=={1,2}
            for end,(f,flag,tags) in chosen.items():
                mate,mflag,mtags=chosen[3-end]
                assert bool(flag&8)==bool(mflag&4),('Mate unmapped flag',name,end)
                assert not flag&2 or (not flag&4 and not mflag&4 and f[2]==mate[2])
                if not mflag&4:
                    assert bool(flag&32)==bool(mflag&16)
                    actual=f[2] if f[6]=='=' else f[6]
                    assert actual==mate[2] and int(f[7])==int(mate[3]),('Mapped mate fields',name,end)
                elif 'lu' in mtags:
                    assert f[6:8]==['*','0']
                expected_tlen=native_tlen(f,mate)
                assert int(f[8])==expected_tlen,('Native 5prime TLEN',name,end,int(f[8]),expected_tlen)
                assert int(f[8])==-int(mate[8]),('TLEN asymmetry',name)
                if flag&2:
                    assert bool(flag&16)!=bool(mflag&16),'Improper orientation labelled proper pair'
            count+=1
            if output is not None:
                output.write(''.join('\t'.join(f)+'\n' for f in rows))
        assert next(left,None) is None and next(right,None) is None,'Missing SAM pairs'
    assert expected_pairs is None or count==expected_pairs,(count,expected_pairs)
    return dict(pairs=count,primary_records=primary,all_records=records,unprojectable_records=unprojectable,reference_sequences=len(refs))

def main():
    p=argparse.ArgumentParser()
    for key in ('sam','r1','r2','fai','out'):p.add_argument('--'+key,required=True)
    p.add_argument('--expected-pairs',type=int)
    p.add_argument('--forward',action='store_true')
    a=p.parse_args(); result=validate(a.sam,a.r1,a.r2,a.fai,a.expected_pairs,sys.stdout if a.forward else None)
    with open(a.out,'x') as out:json.dump(result,out,indent=2)
    print('RAW_SAM_PAIRS_VALID',json.dumps(result,sort_keys=True),file=sys.stderr)

if __name__=='__main__':main()
