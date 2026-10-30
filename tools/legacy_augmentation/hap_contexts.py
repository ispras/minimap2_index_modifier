"""Build indel context contigs for an augmented reference.

Args: REF_FASTA PANEL_VCF OUT_FASTA OUT_CHAIN REF_FAI OUT_BED. Writes the reference plus one
contig per selected indel (and, with HAP=1, haplotype contigs joining co-carried indels),
a chain to the reference, hap-contigs.tsv, tmp-insertions.vcf and vcf-for-indexing.vcf
(panel SNVs projected onto the contigs). Parameters come from environment variables."""
import sys
import string
import collections
DROPPED = collections.Counter()
WRITTEN = [0]
HAPDROPPED = collections.Counter()
HAPWRITTEN = [0]
HAPC = collections.Counter()
import os
import subprocess
import re
import textwrap
import time
import argparse
import gzip
import bisect
import pysam
parser = argparse.ArgumentParser()
parser.add_argument('files', metavar='file', type=str, nargs=6, help='REF_FASTA PANEL_VCF OUT_FASTA OUT_CHAIN REF_FAI OUT_BED')
args = parser.parse_args()
file1, file2, file3, file4, file5, file6 = args.files
file7 = 'tmp-contigs.chain'
input_vcf = file2
res = ''
count = 1
fai_list = {}
contig_half = 100
contig_half_len = contig_half
haplotype_gt = True
_hap_carr = {}

def hap_mask(parts, carriers):
    out = parts[:9]
    any1 = False
    for j, g in enumerate(parts[9:]):
        new = []
        for k, a in enumerate(g.split('|')):
            if a == '.':
                a = '0'
            if carriers is not None and (not (j < len(carriers) and k < len(carriers[j]) and (carriers[j][k] == '1'))):
                a = '0'
            if a == '1':
                any1 = True
            new.append(a)
        out.append('|'.join(new))
    return (out, any1)
start_time = time.time()
import gzip
AF_INDEL = float(os.environ.get('AF_INDEL', '0.1'))
MIN_INDEL_ALLELE = int(os.environ.get('MIN_INDEL_ALLELE', '5'))
MAX_INDEL_ALLELE = int(os.environ.get('MAX_INDEL_ALLELE', '100'))
if not 2 <= MIN_INDEL_ALLELE <= MAX_INDEL_ALLELE:
    raise ValueError('need 2 <= MIN_INDEL_ALLELE <= MAX_INDEL_ALLELE')
FIXDEDUP = os.environ.get('FIXDEDUP', '1') == '1'
HAP = os.environ.get('HAP', '1') == '1'
HAP_GAP = int(os.environ.get('HAP_GAP', '100'))
HAP_MAXSPAN = int(os.environ.get('HAP_MAXSPAN', '400'))
HAP_MAXT = int(os.environ.get('HAP_MAXT', '16'))
HAP_MAXGAPS = int(os.environ.get('HAP_MAXGAPS', '8'))
print('gen_hap: AF_INDEL=%g FIXDEDUP=%d HAP=%d HAP_GAP=%d HAP_MAXSPAN=%d HAP_MAXT=%d HAP_MAXGAPS=%d' % (AF_INDEL, FIXDEDUP, HAP, HAP_GAP, HAP_MAXSPAN, HAP_MAXT, HAP_MAXGAPS))
print('gen_hap: INDEL_ALLELE_LENGTH=%d..%d' % (MIN_INDEL_ALLELE, MAX_INDEL_ALLELE))

def parse_info_field(info_str):
    info = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            k, v = entry.split('=', 1)
            info[k] = v
        else:
            info[entry] = True
    return info

def filter_vcf(input_vcf, output_vcf, af_threshold=0.1):
    seen_positions = set()
    opener = gzip.open if input_vcf.endswith('.gz') else open
    out_opener = gzip.open if output_vcf.endswith('.gz') else open
    with opener(input_vcf, 'rt') as infile, out_opener(output_vcf, 'wt') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, ref = (fields[0], fields[1], fields[3])
            alts = fields[4].split(',')
            refs = fields[3].split(',')
            if 'random' in chrom:
                continue
            if 'chrUn' in chrom:
                continue
            if 'chrM' in chrom:
                continue
            has_insertion = any((len(alt) > len(ref) for alt in alts))
            has_deletion = any((len(alt) < len(ref) for alt in alts))
            if not has_insertion and (not has_deletion):
                continue
            position_key = f'{chrom}_{pos}'
            if position_key in seen_positions:
                continue
            info = parse_info_field(fields[7])
            af_value = None
            if 'AF' in info:
                af_str = info['AF']
                af_parts = af_str.split(',') if af_str else []
                if len(af_parts) >= 1 and af_parts[0] != '':
                    try:
                        af_value = float(af_parts[0])
                    except ValueError:
                        af_value = None
            if af_value is None or af_value <= af_threshold:
                continue
            if not FIXDEDUP:
                seen_positions.add(position_key)
            fields[4] = alts[0]
            fields[3] = refs[0]
            if has_insertion:
                if len(fields[4]) > MAX_INDEL_ALLELE:
                    continue
                if len(fields[4]) < MIN_INDEL_ALLELE:
                    continue
            if has_deletion:
                if len(fields[3]) > MAX_INDEL_ALLELE:
                    continue
                if len(fields[3]) < MIN_INDEL_ALLELE:
                    continue
            if FIXDEDUP:
                seen_positions.add(position_key)
            outfile.write('\t'.join(fields) + '\n')
_r2_chrom = [None, None]

def r2_ref_ok(chrom, pos, ref):
    if _r2_chrom[0] != chrom:
        try:
            _r2_chrom[1] = pysam_fasta.fetch(chrom).upper()
        except Exception:
            _r2_chrom[1] = None
        _r2_chrom[0] = chrom
    seq = _r2_chrom[1]
    if seq is None:
        return False
    b = seq[pos - 1:pos - 1 + len(ref)]
    if len(b) < len(ref):
        return False
    r = ref.upper()
    for i in range(len(b)):
        if b[i] not in 'ACGT':
            continue
        if b[i] != r[i]:
            return False
    return True

def filter_vcf_by_positions(vcf_file, output_file, contig_map, window_left, window_right, contig_to_vcf, contig_geom):
    seen_exact_matches = set()
    alt_length_by_contig = {}
    _op = gzip.open if vcf_file.endswith('.gz') else open
    with _op(vcf_file, 'rt') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            line2 = line
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ok_ref = True
            if len(parts[3]) == len(parts[4]) and len(parts[3]) == 1:
                ok_ref = r2_ref_ok(chrom, pos, parts[3])
                if ok_ref:
                    outfile.write('\t'.join(hap_mask(parts, None)[0]) + '\n' if haplotype_gt else line2)
                else:
                    DROPPED['ref_mismatch_base'] += 1
            if chrom not in contig_map:
                continue
            positions = contig_map[chrom]
            key = (chrom, pos)
            idx = bisect.bisect_left(positions, pos)
            if idx < len(positions) and positions[idx] == pos:
                if key not in seen_exact_matches:
                    modified_parts = parts.copy()
                    modified_parts[0] = f'{chrom}_contig_{positions[idx]}'
                    modified_parts[1] = str(contig_half_len + 1)
                    alt = parts[4]
                    alt_length_by_contig[modified_parts[0]] = len(alt)
                    seen_exact_matches.add(key)
                continue
            left = bisect.bisect_left(positions, pos - window_right)
            right = bisect.bisect_right(positions, pos + window_left)
            if left < right:
                ref = parts[3]
                alt = parts[4]
                if ',' not in alt and len(ref) == 1 and (len(alt) == 1):
                    for i in range(left, right):
                        contig_pos = positions[i]
                        name = f'{chrom}_contig_{contig_pos}'
                        g = contig_geom.get(name)
                        if g is None:
                            DROPPED['no_geometry'] += 1
                            continue
                        w, lr, la, Lc = g
                        d = pos - contig_pos
                        if d < 0:
                            off0 = w + d
                        elif d >= lr:
                            off0 = w + la + (d - lr)
                        else:
                            DROPPED['inside_deleted_segment'] += 1
                            continue
                        if not (0 <= off0 and off0 + len(ref) <= Lc):
                            DROPPED['out_of_bounds'] += 1
                            continue
                        if not ok_ref:
                            DROPPED['ref_mismatch_contig'] += 1
                            continue
                        modified_parts = parts.copy()
                        modified_parts[0] = name
                        modified_parts[1] = str(off0 + 1)
                        if haplotype_gt:
                            if name not in _hap_carr:
                                if len(_hap_carr) > 200000:
                                    _hap_carr.clear()
                                _hap_carr[name] = [g.split('|') for g in contig_to_vcf[name].rstrip('\n').split('\t')[9:]]
                            modified_parts, any1 = hap_mask(modified_parts, _hap_carr[name])
                            if not any1:
                                DROPPED['hap_no_alt_carrier'] += 1
                                continue
                        WRITTEN[0] += 1
                        outfile.write('\t'.join(modified_parts) + '\n')
            if HAP and chrom in hap_win and (',' not in parts[4]) and (len(parts[3]) == 1) and (len(parts[4]) == 1):
                hst, hl = hap_win[chrom]
                hi = bisect.bisect_right(hst, pos) - 1
                while hi >= 0 and hl[hi][0] >= pos - HAP_MAXWIN[0]:
                    s0, e0, name = hl[hi]
                    hi -= 1
                    if not s0 <= pos <= e0:
                        continue
                    w, p1, geo, Lc = hap_geom[name]
                    off0 = w + (pos - p1)
                    inside = False
                    for p, lr, la in geo:
                        if p < pos < p + lr:
                            inside = True
                            break
                        if pos >= p + lr:
                            off0 += la - lr
                    if inside:
                        HAPDROPPED['inside_deleted_segment'] += 1
                        continue
                    if not (0 <= off0 and off0 + len(parts[3]) <= Lc):
                        HAPDROPPED['out_of_bounds'] += 1
                        continue
                    if not ok_ref:
                        HAPDROPPED['ref_mismatch_contig'] += 1
                        continue
                    modified_parts = parts.copy()
                    modified_parts[0] = name
                    modified_parts[1] = str(off0 + 1)
                    if haplotype_gt:
                        modified_parts, any1 = hap_mask(modified_parts, hap_carr[name])
                        if not any1:
                            HAPDROPPED['hap_no_alt_carrier'] += 1
                            continue
                    HAPWRITTEN[0] += 1
                    outfile.write('\t'.join(modified_parts) + '\n')
    print('contig records written: %d' % WRITTEN[0])
    print('contig records dropped: %d' % sum(DROPPED.values()))
    for _k in sorted(DROPPED):
        print('   dropped, reason %-26s %d' % (_k, DROPPED[_k]))
    if HAP:
        print('haplotype contig records written: %d' % HAPWRITTEN[0])
        for _k in sorted(HAPDROPPED):
            print('   hap dropped, reason %-26s %d' % (_k, HAPDROPPED[_k]))
filter_vcf(file2, 'tmp-insertions.vcf', af_threshold=AF_INDEL)
print('Insertions from VCF filtered.')
file2 = 'tmp-insertions.vcf'

def format_header(header, pos):
    match = re.match('>(\\w+):(\\d+)-(\\d+)', header)
    if match:
        chrom = match.groups()[0]
        return f'>{chrom}_contig_{pos}'
    return header
pysam_fasta = pysam.FastaFile(file1)
contig_to_vcf = {}
contig_geom = {}
hap_geom = {}
hap_carr = {}
hap_win = {}
HAP_MAXWIN = [0]

def hap_templates(cl):
    tm = collections.OrderedDict()
    nsamp = len(cl[0][3])
    for j in range(nsamp):
        nph = max((len(r[3][j]) for r in cl))
        for k in range(nph):
            key = tuple((i for i, r in enumerate(cl) if k < len(r[3][j]) and r[3][j][k] == '1'))
            if len(key) >= 2:
                tm[key] = tm.get(key, 0) + 1
    return list(tm.items())

def build_hap_contigs(new):
    t0 = time.time()
    by_chr = collections.OrderedDict()
    for nm, ln in contig_to_vcf.items():
        f = ln.rstrip('\n').split('\t')
        by_chr.setdefault(f[0], []).append((int(f[1]), f[3], f[4], [g.split('|') for g in f[9:]], nm))
    names = set(contig_to_vcf)
    tsv = open('hap-contigs.tsv', 'w')
    tsv.write('#name\tchrom\tp1\tw\tlen\tncar_exact\tncar_all\tindels\n')
    for ch, recs in by_chr.items():
        recs.sort(key=lambda r: r[0])
        clusters, cur, cur_end = ([], [recs[0]], recs[0][0] + len(recs[0][1]) - 1)
        for r in recs[1:]:
            if r[0] - cur_end <= HAP_GAP:
                cur.append(r)
                cur_end = max(cur_end, r[0] + len(r[1]) - 1)
            else:
                clusters.append(cur)
                cur, cur_end = ([r], r[0] + len(r[1]) - 1)
        clusters.append(cur)
        wins = []
        for cl in clusters:
            HAPC['clusters total'] += 1
            if len(cl) < 2:
                HAPC['clusters with 1 indel'] += 1
                continue
            HAPC['clusters with >=2 indels'] += 1
            HAPC['indels in clusters with >=2'] += len(cl)
            span = max((r[0] + len(r[1]) - 1 for r in cl)) - cl[0][0] + 1
            if span > HAP_MAXSPAN:
                HAPC['skipped: cluster longer than HAP_MAXSPAN'] += 1
                HAPC['skipped: indels in such clusters'] += len(cl)
                continue
            tms = hap_templates(cl)
            HAPC['templates with >=2 indels (with carrier)'] += len(tms)
            good = []
            for key, ncar in tms:
                bad = any((cl[key[t + 1]][0] < cl[key[t]][0] + len(cl[key[t]][1]) for t in range(len(key) - 1)))
                if bad:
                    HAPC['skipped: template with overlapping REF'] += 1
                    continue
                if len(key) > HAP_MAXGAPS:
                    HAPC['skipped: template longer than HAP_MAXGAPS indels'] += 1
                    continue
                good.append((key, ncar))
            good.sort(key=lambda kv: (-kv[1], [cl[i][0] for i in kv[0]]))
            if len(good) > HAP_MAXT:
                HAPC['skipped: templates beyond HAP_MAXT'] += len(good) - HAP_MAXT
                good = good[:HAP_MAXT]
            if good:
                HAPC['clusters with haplotype contigs'] += 1
            for K, (key, ncar) in enumerate(good, 1):
                ind = [cl[i] for i in key]
                p1 = ind[0][0]
                w = min(contig_half, p1 - 1)
                parts = [pysam_fasta.fetch(ch, p1 - 1 - w, p1 - 1)]
                for t, (p, ref, alt, gt, _) in enumerate(ind):
                    parts.append(alt)
                    if t + 1 < len(ind):
                        parts.append(pysam_fasta.fetch(ch, p - 1 + len(ref), ind[t + 1][0] - 1))
                    else:
                        parts.append(pysam_fasta.fetch(ch, p - 1 + len(ref), p - 1 + len(ref) + contig_half))
                seq = ''.join(parts)
                name = '%s_contig_%d_h%d' % (ch, p1, K)
                if name in names:
                    raise SystemExit('ERROR: duplicate haplotype contig name: ' + name)
                names.add(name)
                new.write('>%s\n' % name)
                new.write(textwrap.fill(seq, width=70))
                new.write('\n')
                geo = [(p, len(ref), len(alt)) for p, ref, alt, gt, _ in ind]
                hap_geom[name] = (w, p1, geo, len(seq))
                nsamp = len(ind[0][3])
                carr, nall = ([], 0)
                for j in range(nsamp):
                    nph = max((len(r[3][j]) for r in ind))
                    row = []
                    for k in range(nph):
                        c = all((k < len(r[3][j]) and r[3][j][k] == '1' for r in ind))
                        row.append('1' if c else '0')
                        nall += c
                    carr.append(row)
                hap_carr[name] = carr
                eK = ind[-1][0] + len(ind[-1][1])
                wins.append((p1 - w, eK + contig_half - 1, name))
                HAPC['haplotype contigs'] += 1
                HAPC['haplotype contigs with %d indels' % len(ind)] += 1
                tsv.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n' % (name, ch, p1, w, len(seq), ncar, nall, ','.join(('%d:%s:%s' % (p, ref, alt) for p, ref, alt, gt, _ in ind))))
        if wins:
            wins.sort()
            hap_win[ch] = ([x[0] for x in wins], wins)
            HAP_MAXWIN[0] = max(HAP_MAXWIN[0], max((e - s + 1 for s, e, _ in wins)))
    tsv.close()
    print('haplotype contigs: %.1f s' % (time.time() - t0))
    for k in sorted(HAPC):
        print('   hap %-48s %d' % (k, HAPC[k]))
with open(file1, 'r') as fasta, open(file2, 'r') as vcf, open(file3, 'w') as new, open(file4, 'w') as chainfile, open(file5, 'r') as fai, open(file6, 'w') as bedfile, open(file7, 'w') as tmpcontigs:
    while string != '':
        string = fasta.readline()
        if string == '':
            break
        new.write(string)
    string = 'start'
    for line in fai:
        columns = line.strip().split()
        if len(columns) >= 2:
            key = columns[0]
            value = columns[1]
            fai_list[key] = value
    while string != '':
        string = vcf.readline()
        if string == '':
            break
        tmp = string
        full_sequence = []
        if tmp[0] != '#':
            tmp = tmp.split()
            chr = tmp[0]
            alts = tmp[4].split(',')
            refs = tmp[3].split(',')
            ref = refs[0]
            alt = alts[0]
            has_insertion = len(alt) > len(ref)
            has_deletion = len(alt) < len(ref)
            pos = int(tmp[1])
            pos_start = pos - contig_half_len
            contig_half_len = contig_half
            pos_end = pos + contig_half_len
            if pos_start <= 0:
                pos_start = 1
                contig_half_len = pos - pos_start
            bedfile.write(f'{chr}\t{pos - 10}\t{pos + 10}\n')
            left_seq = ''
            if pos_start <= pos - 1:
                try:
                    left_seq = pysam_fasta.fetch(chr, pos_start - 1, pos - 1)
                except Exception:
                    left_seq = ''
            header_line = f'>{chr}:{pos_start}-{pos - 1}'
            contig_header = format_header(header_line, pos)
            new.write(f'{contig_header}\n')
            if left_seq:
                full_sequence.append(left_seq.strip())
            full_sequence.append(tmp[4].strip())
            skip = len(ref) - 1
            right_seq = ''
            if pos + 1 <= pos_end:
                try:
                    right_seq = pysam_fasta.fetch(chr, pos + skip, pos_end + skip)
                except Exception:
                    right_seq = ''
            if right_seq:
                full_sequence.append(right_seq.strip())
            combined_sequence = ''.join(full_sequence)
            wrapped_sequence = textwrap.fill(combined_sequence, width=70)
            new.write(wrapped_sequence)
            new.write('\n')
            contig_header2 = contig_header.split('_')
            contig_header2 = contig_header2[0]
            half = int(contig_half_len)
            alt_s, ref_s = (tmp[4], tmp[3])
            try:
                ref_base = pysam_fasta.fetch(chr, pos_start - 1 + half, pos_start + half)
            except Exception:
                ref_base = ''
            anchor_shared = len(ref_base) == 1 and ref_base.upper() in 'ACGT' and (ref_base.upper() == alt_s[0].upper())
            if anchor_shared:
                first_blk, dt_gap, dq_gap = (half + 1, len(alt_s) - 1, len(ref_s) - 1)
            else:
                first_blk, dt_gap, dq_gap = (half, len(alt_s), len(ref_s))
            last_blk = len(combined_sequence) - half - len(alt_s)
            q_end = pos_start - 1 + first_blk + dq_gap + last_blk
            chainfile.write(f'chain {len(combined_sequence)} {contig_header[1:]} {len(combined_sequence)} + 0 {len(combined_sequence)} {contig_header2[1:]} {fai_list.get(contig_header2[1:])} + {pos_start - 1} {q_end} {count}\n')
            tmpcontigs.write(f'chain {len(combined_sequence)} {contig_header[1:]} {len(combined_sequence)} + 0 {len(combined_sequence)} {contig_header2[1:]} {fai_list.get(contig_header2[1:])} + {pos_start - 1} {q_end} {count}\n')
            count += 1
            chainfile.write(f'{first_blk} {dt_gap}\t{dq_gap}\n{last_blk}\n\n')
            new_header = contig_header[1:]
            contig_to_vcf[new_header] = string
            contig_geom[new_header] = (pos - pos_start, len(tmp[3]), len(tmp[4]), len(combined_sequence))
    if HAP:
        build_hap_contigs(new)
    for chrom, chrom_len in fai_list.items():
        chainfile.write(f'chain {chrom_len} {chrom} {chrom_len} + 0 {chrom_len} {chrom} {chrom_len} + 0 {chrom_len} {count}\n')
        count += 1
        chainfile.write(f'{chrom_len}\n\n')
print('New Fasta and ChainFile generated.')
filename = file4
input_file = input_vcf
output_file = 'vcf-for-indexing.vcf'
contig_map = {}
with open(filename, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('chain'):
            parts = line.split()
            if len(parts) >= 3:
                if parts[-3] != '0':
                    contig_name = parts[2]
                    chromosome = parts[7]
                    match = re.search('_(\\d+)$', contig_name)
                    if match:
                        contig_suffix = int(match.group(1))
                        contig_map.setdefault(chromosome, []).append(contig_suffix)
contig_half_len = contig_half
filter_vcf_by_positions(input_file, output_file, contig_map, contig_half_len, contig_half_len, contig_to_vcf, contig_geom)
print('VCF created.')
