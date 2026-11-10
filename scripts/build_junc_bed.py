#!/usr/bin/env python3
"""minimap2 --junc-bed: the annotation bed12 verbatim plus the cohort's STAR introns as bed6."""
import argparse
import glob
import os
import sys
from collections import defaultdict

# STAR SJ.out.tab strand code -> bed strand (aka - STAR value of 1 corresponds to positive strand)
STRAND = {"1": "+", "2": "-", "0": "."}
# the score column for a BED file must be 0-1000, so the summed read count is capped there
BED_MAX_SCORE = 1000

# Columns from the SJ.out.tab file (that STAR produced) used in this script (0-based; STAR manual, section 4.4)
# Python uses zero-based indexing, so column 1 is accessed using index 0, column 2 using index 1, and so on
SJ_CHROM = 0
SJ_INTRON_START = 1  # first base of the intron, 1-based
SJ_INTRON_END = 2  # last base of the intron, 1-based
SJ_STRAND = 3  # 0 means undefined, 1 means +, 2 means -
SJ_MOTIF = 4  # 0 means non-canonical, 1-6 are codes for different type of canonical motifs
SJ_UNIQUE_READS = 6  # uniquely mapping reads crossing the junction
SJ_MIN_COLUMNS = 7


def main():
    # read the options; only those with "required=TRUE" are necessary arguments for the command line
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotation", required=True) # give an existing BED12 annotation file (from GENCODE)
    ap.add_argument("--sj_folder", required=True) # directory containing STAR junction files
    ap.add_argument("--min_uniq", type=int, default=2) # minimum total uniquely mapped read support
    ap.add_argument("--canonical_only", action="store_true") # keep only junctions with canonical/semi canonical motif evidence
    ap.add_argument("--exclude_samples", default="") # comma-separated sample IDs to exclude
    ap.add_argument("-o", "--out", required=True) # path to output BED file
    args = ap.parse_args()

    # collect the STAR junction files, leaving out excluded samples (sample = file name before the first dot)
    # ex: if i want to exclude sample3.SJ.out.tab, i would use --exclude_samples sample3
    exclude = {x for x in args.exclude_samples.split(",") if x}
    files = [f for f in glob.glob(os.path.join(args.sj_folder, "*.SJ.out.tab"))
             if os.path.basename(f).split(".")[0] not in exclude]
    if not files:
        sys.exit("No *.SJ.out.tab under " + args.sj_folder)

    # goal = sum unique reads per intron across samples and remember whether any sample called its motif canonical
    # each splice junction is identified by 4 characteristics: chrom, start, end, strand) 
    # then, for each junction, the script stores total_unique_reads and motif_status 
    junc = defaultdict(lambda: [0, "."])     # makes a dictionary where each splice junction is identified by chrom/start/end/strand
    for f in files:                          # loops through every STAR junction file
        with open(f) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < SJ_MIN_COLUMNS:
                    continue
                key = (cols[SJ_CHROM], int(cols[SJ_INTRON_START]) - 1, int(cols[SJ_INTRON_END]),    # bed starts are 0-based, STAR starts are 1-based
                       STRAND.get(cols[SJ_STRAND], "."))                                            # ex: STAR says start is 1001, but BED says 1000
                junc[key][0] += int(cols[SJ_UNIQUE_READS])      # summing all samples to give a total # of unique supporting reads per junction
                if cols[SJ_MOTIF] != "0":                       # any nonzero STAR motif is fine, because 0 means non-canoincal, other numbers are canonical with different levels
                    junc[key][1] = "canonical"

    # filters to keep only introns with enough unique reads (and a canonical motif when asked) as bed6 lines
    # the >= 2 unique reads + canonical filter is our own choice, not taken from a paper
    # each BED6 entry represents a single intron with 6 columns (chromosome, start, end, name, supporting read count, strand)
    short_read_lines = []
    for (chrom, start, end, strand), (uniq, motif) in junc.items():
        if uniq < args.min_uniq:
            continue
        if args.canonical_only and motif != "canonical":
            continue
        score = min(uniq, BED_MAX_SCORE)
        short_read_lines.append("\t".join([chrom, str(start), str(end), ".", str(score), strand]))  # this converts STAR junctions to BED6

    # Write the final output file which has two kinds of records
    # part 1 is existing BED12 annotations (12 columns), copied directly as is, preserving original record
    # part 2 is the STAR-derived BED6 junctions (that we filtered), see above for the 6 columns used 
    # this combined junction BED file is used as input to minimap2 for --junc-bed
    with open(args.out, "w") as out, open(args.annotation) as anno:
        for line in anno:
            out.write(line)
        out.write("\n".join(short_read_lines) + "\n")
    sys.stderr.write("short-read junctions kept: %d -> %s\n" % (len(short_read_lines), args.out))


if __name__ == "__main__":
    main()
