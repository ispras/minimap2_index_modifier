"""Pure representation helpers for GRCh38-backed haplotype augmentation."""
from dataclasses import dataclass
import re


def callable_blocks(start, cs):
    """Target bases actually aligned to source bases, excluding deletions."""
    cursor, blocks = start, []
    tokens = re.findall(r':[0-9]+|\*[A-Za-z][A-Za-z]|[=+\-][A-Za-z]+', cs)
    if ''.join(tokens) != cs:
        raise ValueError("Unsupported cs token")
    for tok in tokens:
        op = tok[0]
        n = int(tok[1:]) if op == ':' else (1 if op == '*' else len(tok)-1)
        if op in ':=*':
            if blocks and blocks[-1][1] == cursor:
                blocks[-1] = (blocks[-1][0], cursor+n)
            else:
                blocks.append((cursor, cursor+n))
            cursor += n
        elif op == '-':
            cursor += n
    return blocks, cursor


@dataclass(frozen=True, order=True)
class Event:
    chrom: str
    pos: int
    ref: str
    alt: str

    def validate(self, reference):
        if self.pos < 1 or not self.ref or not self.alt:
            raise ValueError("Invalid VCF coordinates or empty allele")
        if set(self.ref + self.alt) - set("ACGT"):
            raise ValueError("Non-ACGT allele")
        if reference[self.pos - 1:self.pos - 1 + len(self.ref)] != self.ref:
            raise ValueError("REF mismatch")


def merge_carriers(labelled):
    """Deduplicate events, never positions; retain explicit source identity."""
    merged = {}
    for origin, event in labelled:
        if origin not in ("hap1", "hap2"):
            raise ValueError("Unknown carrier")
        merged.setdefault(event, set()).add(origin)
    return merged


def construct(reference, events, flank=100):
    """Build one compatible template and its exact projection geometry."""
    if not events or len(events) > 8:
        raise ValueError("Template gap count outside supported range")
    events = sorted(events)
    if len({e.chrom for e in events}) != 1:
        raise ValueError("Mixed chromosomes")
    start = max(0, events[0].pos - 1 - flank)
    cursor = start
    pieces, gaps = [], []
    offset = 0
    for event in events:
        event.validate(reference)
        site = event.pos - 1
        if site < cursor:
            raise ValueError("Overlapping reference alleles")
        if len(event.ref) == len(event.alt) or event.ref[0] != event.alt[0]:
            raise ValueError("Expected anchored indel")
        between = reference[cursor:site]
        pieces.extend((between, event.alt))
        offset += len(between)
        gaps.append((offset + 1, len(event.alt) - 1, len(event.ref) - 1))
        offset += len(event.alt)
        cursor = site + len(event.ref)
    end = min(len(reference), cursor + flank)
    pieces.append(reference[cursor:end])
    return "".join(pieces), start, end, gaps


def common_carriers(events, carrier_map):
    """Reject artificial templates assembled from incompatible source paths."""
    carriers = {"hap1", "hap2"}
    for event in events:
        carriers.intersection_update(carrier_map[event])
    return carriers
