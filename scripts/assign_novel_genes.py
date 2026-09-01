#!/usr/bin/env python
"""
assign_novel_genes.py

Assign a gene_id to novel transcripts by same-strand exonic overlap with the
reference, instead of inheriting it from gffcompare's cmp_ref.

cmp_ref is a nearest-TRANSCRIPT annotation, not a host-GENE annotation. gffcompare
emits one for class codes that share no exonic sequence with the reference gene at all
(i = inside a reference intron, y = a reference gene inside our intron) and for codes
that are on the opposite strand (x, s); for the rest it picks a single best transcript
by class-code priority, which is not the gene sharing the most sequence. Inheriting
gene_id from it puts antisense and intronic transcripts into the isoform pool of a gene
they are not isoforms of, which corrupts every usage denominator computed from it.

the rule here:
  1. exonic-base overlap against every same-strand reference gene; most bases wins
  2. no same-strand exonic overlap -> a fresh locus, not the nearest gene
  3. transcripts that reach step 2 and overlap each other are pooled into one locus,
     so an antisense locus with several isoforms stays testable rather than becoming
     a set of single-isoform genes that dtu drops

bridging is recorded, not resolved: when a second same-strand gene whose span is
disjoint from the assigned gene also holds exonic overlap, its id goes in a `bridging`
attribute. read-through transcription is real, and picking one of the two loci is what
hides the other.

as a library, call assign(); as a script, rewrite gene_id on a built reference GTF.
transcript ids and sequences are untouched, so a rewritten reference is
quantification-neutral: existing salmon output stays valid and only the tx2gene join
downstream changes.
"""
import argparse
import bisect
import gzip
import re
from collections import defaultdict

ATTR = re.compile(r'(\S+) "([^"]*)"')


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def merge_intervals(iv):
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


class ExonIndex:
    """same-strand exonic overlap lookup over a set of genes"""

    def __init__(self):
        self._raw = defaultdict(list)
        self.gene_span = {}
        self._built = False

    def add_gene(self, gene_id, chrom, strand, exons):
        merged = merge_intervals(exons)
        for s, e in merged:
            self._raw[(chrom, strand)].append((s, e, gene_id))
        lo = min(s for s, _ in merged)
        hi = max(e for _, e in merged)
        prev = self.gene_span.get(gene_id)
        if prev:
            lo, hi = min(lo, prev[1]), max(hi, prev[2])
        self.gene_span[gene_id] = (chrom, lo, hi)
        self._built = False

    def build(self):
        self.idx, self.starts, self.maxend = {}, {}, {}
        for key, lst in self._raw.items():
            lst.sort()
            self.idx[key] = lst
            self.starts[key] = [x[0] for x in lst]
            run, cummax = 0, []
            for _, e, _g in lst:
                run = max(run, e)
                cummax.append(run)
            self.maxend[key] = cummax
        self._built = True

    def overlap(self, chrom, strand, exons):
        """gene_id -> exonic bases shared, same strand only"""
        if not self._built:
            self.build()
        key = (chrom, strand)
        if key not in self.idx:
            return {}
        lst, st, me = self.idx[key], self.starts[key], self.maxend[key]
        tot = defaultdict(int)
        for s, e in exons:
            j = bisect.bisect_right(st, e) - 1
            while j >= 0 and me[j] >= s:
                rs, re_, g = lst[j]
                ov = min(e, re_) - max(s, rs) + 1
                if ov > 0:
                    tot[g] += ov
                j -= 1
        return dict(tot)


def _spans_disjoint(a, b):
    if a[0] != b[0]:
        return True
    return a[2] < b[1] or b[2] < a[1]


def _pool_unassigned(unassigned, prefix, start_at=0):
    """single-linkage pooling of leftover transcripts by same-strand exonic overlap"""
    by_key = defaultdict(list)
    for tid, chrom, strand, exons in unassigned:
        by_key[(chrom, strand)].append((min(s for s, _ in exons),
                                        max(e for _, e in exons), tid, exons))
    out, n = {}, start_at
    for key in sorted(by_key):
        items = sorted(by_key[key])
        cluster, hi = [], None
        groups = []
        for lo, end, tid, exons in items:
            if cluster and lo > hi:
                groups.append(cluster)
                cluster = []
            cluster.append((tid, exons))
            hi = end if hi is None or lo > hi else max(hi, end)
        if cluster:
            groups.append(cluster)
        for grp in groups:
            n += 1
            gid = "%sG_%06d" % (prefix, n)
            for tid, _ in grp:
                out[tid] = gid
    return out, n


def assign(novel, index, prefix="NOVEL"):
    """novel: list of (tid, chrom, strand, [(s,e)...]).
    returns tid -> (gene_id, how, bridging_list)"""
    result, unassigned = {}, []
    for tid, chrom, strand, exons in novel:
        ov = index.overlap(chrom, strand, exons)
        if not ov:
            unassigned.append((tid, chrom, strand, exons))
            continue
        best = max(sorted(ov), key=lambda g: ov[g])
        span = index.gene_span[best]
        bridging = sorted(g for g in ov
                          if g != best and _spans_disjoint(span, index.gene_span[g]))
        result[tid] = (best, "exonic_overlap", bridging)
    pooled, _ = _pool_unassigned(unassigned, prefix)
    for tid, gid in pooled.items():
        result[tid] = (gid, "new_locus", [])
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True,
                    help="built reference GTF: reference records plus novel records")
    ap.add_argument("--novel-source", default="ONT",
                    help="GTF source column marking novel records (default ONT)")
    ap.add_argument("--prefix", default="XNOVEL", help="prefix for fresh loci")
    ap.add_argument("--mapping-out", help="TSV of old vs new gene_id per novel transcript")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    index = ExonIndex()
    ref_exons, ref_meta = defaultdict(list), {}
    nov_exons, nov_meta = defaultdict(list), {}

    with opener(args.gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "exon":
                continue
            gid = attr(c[8], "gene_id")
            s, e = int(c[3]), int(c[4])
            if c[1] == args.novel_source:
                tid = attr(c[8], "transcript_id")
                nov_exons[tid].append((s, e))
                nov_meta.setdefault(tid, (c[0], c[6], gid))
            else:
                ref_exons[gid].append((s, e))
                ref_meta[gid] = (c[0], c[6])

    for gid, exons in ref_exons.items():
        chrom, strand = ref_meta[gid]
        index.add_gene(gid, chrom, strand, exons)
    index.build()

    novel = [(tid, nov_meta[tid][0], nov_meta[tid][1], sorted(ex))
             for tid, ex in nov_exons.items()]
    assigned = assign(novel, index, prefix=args.prefix)

    changed = sum(1 for tid, (gid, _, _) in assigned.items() if gid != nov_meta[tid][2])
    fresh = sum(1 for v in assigned.values() if v[1] == "new_locus")
    bridged = sum(1 for v in assigned.values() if v[2])

    if args.mapping_out:
        with open(args.mapping_out, "w") as out:
            out.write("transcript_id\told_gene_id\tnew_gene_id\tassignment\tbridging\n")
            for tid in sorted(assigned):
                gid, how, br = assigned[tid]
                out.write("%s\t%s\t%s\t%s\t%s\n"
                          % (tid, nov_meta[tid][2], gid, how, ",".join(br)))

    op = gzip.open if args.output.endswith(".gz") else open
    with opener(args.gtf) as fh, op(args.output, "wt") as out:
        for line in fh:
            if line[0] == "#":
                out.write(line)
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[1] != args.novel_source:
                out.write(line)
                continue
            tid = attr(c[8], "transcript_id")
            if tid not in assigned:
                out.write(line)
                continue
            gid, how, br = assigned[tid]
            kept = dict(ATTR.findall(c[8]))
            old = kept.get("gene_id", "")
            fields = ['gene_id "%s"' % gid, 'transcript_id "%s"' % tid]
            if c[2] == "transcript":
                fields += ['cmp_ref_gene "%s"' % old,
                           'gene_assignment "%s"' % how,
                           'bridging "%s"' % ",".join(br)]
                for k in ("ref_gene", "class_code", "cohorts"):
                    if k in kept:
                        fields.append('%s "%s"' % (k, kept[k]))
            c[8] = "; ".join(fields) + ";"
            out.write("\t".join(c) + "\n")

    print("* %d novel transcripts: %d gene_id changed, %d on fresh loci, %d bridging -> %s"
          % (len(assigned), changed, fresh, bridged, args.output))


if __name__ == "__main__":
    main()
