#!/usr/bin/env python
"""
wobble.py -- wobble-tolerant isoform clustering for the per-cohort 2/3 consensus
filter (consensus_from_gffcompare.py). NOT used across cohorts: the cross-cohort
union uses gffcompare's exact intron-chain match (see cross_cohort_reference.py),
because there the goal is maximal retention -- only true structural duplicates are
collapsed. Wobble is for WITHIN a cohort, where the 3 assemblers disagree by a few
bp on identical alignments.

gffcompare matches transcripts by EXACT intron chain, so two tools that call the
same isoform with a few-bp shift at a splice site land in separate transfrags and
never reach consensus. Here we group transcripts whose intron chains agree within
+/- wobble bp (Union-Find over a bucketed, sorted sweep -- same idea IsoQuant/isomatch
use). Single-exon transcripts have no chain and are grouped by terminal overlap
within wobble instead.

Adjacency sweep after sorting by the full junction tuple unions near-identical
chains; transitive chains fall out of Union-Find. For small wobble this matches the
intended clustering; very large wobble can over-chain (don't set it large).

Wobble is standard for long-read isoform matching: IsoQuant uses delta 0/4/6/12
(exact/precise/PacBio / default/ONT / loose; Prjibelski 2023), Cupcake
max_fuzzy_junction=5, TAMA calls it "wobble". These are INTERNAL splice-site
tolerances (small, ~3-6bp); transcript-end fuzziness (TSS/TES) is much larger
(~50bp, e.g. SQANTI3 reference-match) -- so this is applied to junctions only for
multi-exon transcripts. Single-exon transcripts (no junctions) fall back to terminal
overlap within the SAME small tolerance, which is deliberately conservative (won't
over-merge distinct single-exon genes; may leave a few single-exon near-dups).
"""
import re
from collections import defaultdict


def parse_attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def intron_chains(gtf_path):
    """transcript_id -> (chrom, strand, introns_tuple, (start, end)).

    introns_tuple is the flat (donor, acceptor, donor, acceptor, ...) of internal
    junctions; empty for single-exon. (start, end) is the transcript span (for
    single-exon overlap matching).
    """
    exons = defaultdict(list)
    meta = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "exon":
                continue
            tid = parse_attr(c[8], "transcript_id")
            if not tid:
                continue
            exons[tid].append((int(c[3]), int(c[4])))
            meta[tid] = (c[0], c[6])
    chains = {}
    for tid, ex in exons.items():
        ex.sort()
        introns = tuple(x for i in range(len(ex) - 1) for x in (ex[i][1], ex[i + 1][0]))
        chrom, strand = meta[tid]
        chains[tid] = (chrom, strand, introns, (ex[0][0], ex[-1][1]))
    return chains


class _UF:
    def __init__(self, ids):
        self.p = {i: i for i in ids}

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def cluster_by_wobble(chains, wobble):
    """Return a list of clusters (each a list of transcript_ids).

    Multi-exon: same chrom/strand/intron-count and every junction within wobble.
    Single-exon: same chrom/strand and both ends within wobble. wobble<=0 puts
    every transcript in its own singleton cluster (exact behaviour upstream).
    """
    uf = _UF(list(chains.keys()))
    if wobble > 0:
        multi = defaultdict(list)
        single = defaultdict(list)
        for tid, (chrom, strand, introns, span) in chains.items():
            if introns:
                multi[(chrom, strand, len(introns))].append(tid)
            else:
                single[(chrom, strand)].append(tid)
        for tids in multi.values():
            tids.sort(key=lambda t: chains[t][2])
            for i in range(1, len(tids)):
                a, b = tids[i - 1], tids[i]
                if all(abs(x - y) <= wobble for x, y in zip(chains[a][2], chains[b][2])):
                    uf.union(a, b)
        for tids in single.values():
            tids.sort(key=lambda t: chains[t][3])
            for i in range(1, len(tids)):
                a, b = tids[i - 1], tids[i]
                (sa, ea), (sb, eb) = chains[a][3], chains[b][3]
                if abs(sa - sb) <= wobble and abs(ea - eb) <= wobble:
                    uf.union(a, b)
    groups = defaultdict(list)
    for tid in chains:
        groups[uf.find(tid)].append(tid)
    return list(groups.values())
