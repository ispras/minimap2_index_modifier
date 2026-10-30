#!/usr/bin/env python3
"""GRCh38 augmentation adapter; mapping algorithms are not changed."""
import argparse
import bisect
import collections
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import pysam
from panel_core import Event, construct, common_carriers, callable_blocks


def covered(intervals, pos, end):
    i = bisect.bisect_right(intervals, (pos, 10**12)) - 1
    return i >= 0 and intervals[i][0] <= pos and end <= intervals[i][1]


def main():
    ap = argparse.ArgumentParser()
    for key in ("ref", "h1", "h2", "paf1", "paf2", "pafsel", "out"):
        ap.add_argument("--" + key, required=True)
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        raise SystemExit("Refusing existing output directory")
    out.mkdir(parents=True)
    spec = importlib.util.spec_from_file_location("pafsel", a.pafsel)
    ps = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ps)
    ref = pysam.FastaFile(a.ref)
    coverage, stats = {}, collections.Counter()
    for h, pat in enumerate((a.paf1, a.paf2)):
        for chrom in ref.references:
            p = pat.replace("{chrom}", chrom)
            if not Path(p).exists():
                coverage[h, chrom] = []
                continue
            kept, st = ps.select(ps.read_paf(p))
            blocks = []
            for r in kept:
                part, end = callable_blocks(r["ts"], r["cs"])
                if end != r["te"]:
                    raise ValueError("PAF/cs target span mismatch")
                blocks.extend(part)
            coverage[h, chrom] = sorted(blocks)
            for k, v in st.items():
                stats[f"hap{h+1}_paf_{k}"] += v
    indels, snvs = collections.defaultdict(dict), collections.defaultdict(dict)
    cache_chr, cache_seq = None, None
    for h, path in enumerate((a.h1, a.h2)):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip().split("\t")
                chrom, pos, r, alt = f[0], int(f[1]), f[3].upper(), f[4].upper()
                if chrom != cache_chr:
                    cache_chr, cache_seq = chrom, ref.fetch(chrom).upper()
                e = Event(chrom, pos, r, alt)
                e.validate(cache_seq)
                if r == alt:
                    raise ValueError("Reference-identical input allele")
                if len(r) == len(alt) == 1:
                    stats[f"hap{h+1}_snv_input"] += 1
                    if not all(covered(coverage[k, chrom], pos-1, pos) for k in (0, 1)):
                        stats[f"hap{h+1}_snv_omitted_coverage"] += 1
                        continue
                    snvs[chrom][e] = snvs[chrom].get(e, 0) | (1 << h)
                elif len(r) != len(alt):
                    stats[f"hap{h+1}_indel_input"] += 1
                    if not (2 <= max(len(r), len(alt)) <= 5001):
                        stats[f"hap{h+1}_indel_outside_cap"] += 1
                        continue
                    if min(len(r), len(alt)) != 1 or r[0] != alt[0]:
                        raise ValueError("Unexpected non-simple anchored indel")
                    indels[chrom].setdefault(e, set()).add("hap" + str(h+1))
                else:
                    raise ValueError("Unexpected multibase substitution")
    with open(out / "coverage.json", "w") as fh:
        json.dump({f"hap{h+1}:{c}": v for (h, c), v in coverage.items()}, fh)
    shutil.copyfile(a.ref, out / "aug.fasta")
    shutil.copyfile(a.ref + ".fai", out / "base.fasta.fai")
    vcf = open(out / "modifier.vcf", "w")
    hd = open(out / "hapdp.tsv", "w")
    meta = open(out / "contigs.jsonl", "w")
    lm = open(out / "aug.liftmap", "w")
    fa = open(out / "aug.fasta", "a")
    vcf.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=Source-path carrier>\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHAP1\tHAP2\n")
    hd.write("#hapdp v1\n#n_hap\t2\n#hap\t0\tHAP1:0\n#hap\t1\tHAP2:0\n")

    def write_snv(name, pos, r, alt, mask):
        if not mask:
            return
        vcf.write(f"{name}\t{pos}\t.\t{r}\t{alt}\t.\tPASS\t.\tGT\t{mask & 1}\t{(mask >> 1) & 1}\n")
        hd.write(f"{name}\t{pos}\t{alt}\t0000000000000000\t{mask:016x}\n")
        stats["modifier_records"] += 1

    for chrom in ref.references:
        seq = ref.fetch(chrom).upper()
        variants = sorted(snvs[chrom])
        snvpos = [e.pos for e in variants]
        seen = {}
        for e in variants:
            mask = snvs[chrom][e]
            for h in (0, 1):
                if mask & (1 << h):
                    key = (e.pos, h)
                    if key in seen and seen[key] != e.alt:
                        raise ValueError("Conflicting source SNV alleles")
                    seen[key] = e.alt
            write_snv(chrom, e.pos, e.ref, e.alt, mask)
        carrier_map = indels[chrom]
        events = sorted(carrier_map)
        templates = {tuple([e]): carrier_map[e] for e in events}
        clusters, cur, end = [], [], -1
        for e in events:
            if cur and e.pos - end > 100:
                clusters.append(cur)
                cur = []
            cur.append(e)
            end = max(end, e.pos + len(e.ref) - 1)
        if cur:
            clusters.append(cur)
        for cl in clusters:
            if max(e.pos + len(e.ref)-1 for e in cl) - cl[0].pos + 1 > 400:
                stats["clusters_over_400"] += 1
                continue
            for origin in ("hap1", "hap2"):
                t = tuple(e for e in cl if origin in carrier_map[e])
                if not 2 <= len(t) <= 8:
                    continue
                if any(y.pos < x.pos + len(x.ref) for x, y in zip(t, t[1:])):
                    stats["templates_overlap_rejected"] += 1
                    continue
                templates[t] = common_carriers(t, carrier_map)
        names = collections.Counter()
        geometry_seen = set()
        for template, origins in sorted(templates.items()):
            sequence, start, stop, gaps = construct(seq, template)
            geometry = (start, stop, sequence, tuple(gaps))
            if geometry in geometry_seen:
                raise ValueError("Duplicate template geometry needs explicit carrier merge")
            geometry_seen.add(geometry)
            pos = template[0].pos
            names[pos] += 1
            name = f"{chrom}_contig_{pos}_h{names[pos]}"
            fa.write(f">{name}\n")
            fa.write("\n".join(sequence[i:i+70] for i in range(0, len(sequence), 70)) + "\n")
            j, ins, dele = gaps[0]
            fields = [name, str(ins), str(dele), str(j)]
            if len(gaps) > 1:
                fields.append(";".join(f"{j}:{i}:{d}" for j, i, d in gaps[1:]))
            lm.write("\t".join(fields) + "\n")
            meta.write(json.dumps(dict(name=name, chrom=chrom, start=start, end=stop,
                length=len(sequence), origins=sorted(origins), events=[e.__dict__ for e in template], gaps=gaps)) + "\n")
            stats["contigs"] += 1
            stats["multi_event_contigs" if len(template)>1 else "single_event_contigs"] += 1
            mask = sum(1 << h for h in (0, 1) if f"hap{h+1}" in origins)
            for k in range(bisect.bisect_left(snvpos, start+1), bisect.bisect_right(snvpos, stop)):
                e = variants[k]
                offset = e.pos-1-start
                inside = False
                for t in template:
                    if t.pos <= e.pos < t.pos + len(t.ref):
                        inside = True
                        break
                    if e.pos >= t.pos + len(t.ref):
                        offset += len(t.alt)-len(t.ref)
                if inside:
                    continue
                if sequence[offset] != e.ref:
                    raise ValueError("Projected SNV REF mismatch")
                write_snv(name, offset+1, e.ref, e.alt, snvs[chrom][e] & mask)
    for fh in (fa, lm, meta, vcf, hd):
        fh.close()
    with open(out / "panel.stats.json", "w") as fh:
        json.dump(dict(stats), fh, indent=2, sort_keys=True)
    print("PANEL_GENERATED", json.dumps(dict(stats), sort_keys=True))


if __name__ == "__main__":
    main()
