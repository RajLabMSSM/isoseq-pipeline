#!/usr/bin/env python3
"""minimap2 --junc-bed: the annotation bed12 verbatim plus the cohort's STAR introns as bed6."""
import argparse
import glob
import os
import sys
from collections import defaultdict

# STAR SJ.out.tab strand code -> bed strand
STRAND = {"1": "+", "2": "-", "0": "."}
# bed score column must be 0-1000, so the summed read count is capped there
BED_MAX_SCORE = 1000

# STAR SJ.out.tab columns used here (0-based; STAR manual, section 4.4)
SJ_CHROM = 0
SJ_INTRON_START = 1  # first intron base, 1-based
SJ_INTRON_END = 2  # last intron base, 1-based
SJ_STRAND = 3  # 0 undefined, 1 +, 2 -
SJ_MOTIF = 4  # 0 non-canonical, 1-6 canonical / semi-canonical motifs
SJ_UNIQUE_READS = 6  # uniquely mapping reads crossing the junction
SJ_MIN_COLUMNS = 7


def main():
    # read the options
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--sj_folder", required=True)
    ap.add_argument("--min_uniq", type=int, default=2)
    ap.add_argument("--canonical_only", action="store_true")
    ap.add_argument("--exclude_samples", default="")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    # collect the STAR junction files, leaving out excluded samples (sample = file name before the first dot)
    exclude = {x for x in args.exclude_samples.split(",") if x}
    files = [f for f in glob.glob(os.path.join(args.sj_folder, "*.SJ.out.tab"))
             if os.path.basename(f).split(".")[0] not in exclude]
    if not files:
        sys.exit("No *.SJ.out.tab under " + args.sj_folder)

    # sum unique reads per intron across samples and remember whether any sample called its motif canonical
    # (chrom, 0-based start, end, strand) -> [summed unique reads, motif]
    junc = defaultdict(lambda: [0, "."])
    for f in files:
        with open(f) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < SJ_MIN_COLUMNS:
                    continue
                # bed starts are 0-based, STAR starts are 1-based
                key = (cols[SJ_CHROM], int(cols[SJ_INTRON_START]) - 1, int(cols[SJ_INTRON_END]),
                       STRAND.get(cols[SJ_STRAND], "."))
                junc[key][0] += int(cols[SJ_UNIQUE_READS])
                if cols[SJ_MOTIF] != "0":
                    junc[key][1] = "canonical"

    # keep introns with enough unique reads (and a canonical motif when asked) as bed6 lines
    short_read_lines = []
    for (chrom, start, end, strand), (uniq, motif) in junc.items():
        if uniq < args.min_uniq:
            continue
        if args.canonical_only and motif != "canonical":
            continue
        score = min(uniq, BED_MAX_SCORE)
        short_read_lines.append("\t".join([chrom, str(start), str(end), ".", str(score), strand]))

    # write the annotation bed as is, then the short-read introns after it
    with open(args.out, "w") as out, open(args.annotation) as anno:
        for line in anno:
            out.write(line)
        out.write("\n".join(short_read_lines) + "\n")
    sys.stderr.write("short-read junctions kept: %d -> %s\n" % (len(short_read_lines), args.out))


if __name__ == "__main__":
    main()
