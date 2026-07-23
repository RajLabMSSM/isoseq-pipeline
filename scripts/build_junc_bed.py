#!/usr/bin/env python3
"""
build_junc_bed.py  —  minimap2 --junc-bed from annotation + short-read junctions

Concatenates the reference annotation (BED12, passed through verbatim) with
empirical STAR junctions (SJ.out.tab -> BED6 introns) so minimap2 snaps noisy
long reads to splice sites validated by both. minimap2 2.31 reads mixed column
counts per line: >=12 cols -> BED12 (introns from blocks); <12 cols -> the line
is one intron, strand taken from col6 (hence BED6, not a bare 5-col file).

Usage:
    build_junc_bed.py --annotation anno.bed --sj_folder dir/ -o combined.bed
                      [--min_uniq 2] [--canonical_only]
"""

import argparse
import glob
import os
import sys
from collections import defaultdict

# STAR SJ.out.tab columns: 1 chrom, 2 intron_start(1-based), 3 intron_end(1-based),
# 4 strand(0/1/2), 5 motif(0=noncanonical), 6 annotated, 7 unique reads, 8 multi, 9 overhang
STRAND = {"1": "+", "2": "-", "0": "."}


def collect_sj(sj_folder, min_uniq, canonical_only):
    """Aggregate introns across all SJ.out.tab; return BED6 lines."""
    junc = defaultdict(lambda: [0, "."])   # (chrom, start0, end, strand) -> [sum_uniq, motif_ok]
    files = glob.glob(os.path.join(sj_folder, "*.SJ.out.tab"))
    if not files:
        sys.exit("No *.SJ.out.tab under " + sj_folder)
    for f in files:
        with open(f) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 7:
                    continue
                strand = STRAND.get(c[3], ".")
                key = (c[0], int(c[1]) - 1, int(c[2]), strand)   # 0-based intron
                junc[key][0] += int(c[6])
                if c[4] != "0":
                    junc[key][1] = "canonical"

    rows = []
    for (chrom, st, en, strand), (uniq, motif) in junc.items():
        if uniq < min_uniq:
            continue
        if canonical_only and motif != "canonical":
            continue
        rows.append("\t".join([chrom, str(st), str(en), ".",
                               str(min(uniq, 1000)), strand]))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotation", required=True, help="reference junction BED12")
    p.add_argument("--sj_folder",  required=True, help="folder of STAR *.SJ.out.tab")
    p.add_argument("--min_uniq",   type=int, default=2,
                   help="min summed unique reads to keep a short-read junction")
    p.add_argument("--canonical_only", action="store_true",
                   help="drop short-read junctions with non-canonical motif")
    p.add_argument("-o", "--out", required=True)
    args = p.parse_args()

    sr = collect_sj(args.sj_folder, args.min_uniq, args.canonical_only)
    with open(args.out, "w") as out:
        with open(args.annotation) as anno:
            for line in anno:
                out.write(line)
        out.write("\n".join(sr) + "\n")

    sys.stderr.write("short-read junctions kept: %d -> %s\n" % (len(sr), args.out))


if __name__ == "__main__":
    main()
