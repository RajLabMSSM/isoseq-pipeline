#!/usr/bin/env python3
"""minimap2 --junc-bed: the annotation bed12 verbatim plus the cohort's STAR introns as bed6."""
import argparse
import glob
import os
import sys
from collections import defaultdict

STRAND = {"1": "+", "2": "-", "0": "."}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--sj_folder", required=True)
    ap.add_argument("--min_uniq", type=int, default=2)
    ap.add_argument("--canonical_only", action="store_true")
    ap.add_argument("--exclude_samples", default="")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    exclude = {x for x in args.exclude_samples.split(",") if x}
    files = [f for f in glob.glob(os.path.join(args.sj_folder, "*.SJ.out.tab"))
             if os.path.basename(f).split(".")[0] not in exclude]
    if not files:
        sys.exit("No *.SJ.out.tab under " + args.sj_folder)

    # (chrom, 0-based start, end, strand) -> [summed unique reads, motif]
    junc = defaultdict(lambda: [0, "."])
    for f in files:
        with open(f) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 7:
                    continue
                key = (c[0], int(c[1]) - 1, int(c[2]), STRAND.get(c[3], "."))
                junc[key][0] += int(c[6])
                if c[4] != "0":
                    junc[key][1] = "canonical"

    sr = ["\t".join([chrom, str(st), str(en), ".", str(min(uniq, 1000)), strand])
          for (chrom, st, en, strand), (uniq, motif) in junc.items()
          if uniq >= args.min_uniq and (motif == "canonical" or not args.canonical_only)]
    with open(args.out, "w") as out, open(args.annotation) as anno:
        for line in anno:
            out.write(line)
        out.write("\n".join(sr) + "\n")
    sys.stderr.write("short-read junctions kept: %d -> %s\n" % (len(sr), args.out))


if __name__ == "__main__":
    main()
