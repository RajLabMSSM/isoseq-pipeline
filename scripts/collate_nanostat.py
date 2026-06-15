#!/usr/bin/env python
"""
collate_nanostat.py

Collate per-sample NanoStat --tsv outputs into one tidy table for plotting
(e.g. Glinos et al. GTEx long-read Supp. Fig. 1A/B: raw-vs-aligned read number
and median read length, one point per sample).

NanoStat --tsv writes a 2-column key->value table per file:
    Metrics             dataset
    number_of_reads     12345
    median_read_length  2499.0
    ...
This script globs the per-sample raw + aligned files written by the qc pipeline
    <sample>.nanostat.raw.tsv      (NanoStat --fastq on the raw input)
    <sample>.nanostat.aligned.tsv  (NanoStat --bam  on the aligned primaries)
and emits one row per (sample, stage) with the summary metrics as columns.

Mirrors collate_read_lengths.R's interface: --inFolder to search, -o for output.
Brooke Friedman
"""
import argparse
import glob
import os

# NanoStat --tsv keys we keep (others, e.g. per-read top-5 lists, are dropped)
KEEP = [
    "number_of_reads",
    "number_of_bases",
    "median_read_length",
    "mean_read_length",
    "read_length_stdev",
    "n50",
    "mean_qual",
    "median_qual",
]


def parse_nanostat(path):
    """Return {metric: value} for the KEEP metrics found in a NanoStat --tsv file."""
    vals = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            key, val = parts
            if key in KEEP:
                vals[key] = val
    return vals


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inFolder", default=".",
                    help="root dir to search recursively for *.nanostat.{raw,aligned}.tsv")
    ap.add_argument("-o", "--output", default="nanostat_collated.tsv",
                    help="output .tsv path")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.inFolder, "**", "*.nanostat.raw.tsv"),
                             recursive=True)
                   + glob.glob(os.path.join(args.inFolder, "**", "*.nanostat.aligned.tsv"),
                               recursive=True))
    print("* %d NanoStat files found under %s" % (len(files), args.inFolder))

    rows = []
    for f in files:
        base = os.path.basename(f)
        # <sample>.nanostat.<stage>.tsv
        stage = "aligned" if base.endswith(".nanostat.aligned.tsv") else "raw"
        sample = base.replace(".nanostat.%s.tsv" % stage, "")
        vals = parse_nanostat(f)
        rows.append((sample, stage, vals))

    header = ["sample", "stage"] + KEEP
    with open(args.output, "w") as out:
        out.write("\t".join(header) + "\n")
        for sample, stage, vals in rows:
            out.write("\t".join([sample, stage] + [vals.get(k, "NA") for k in KEEP]) + "\n")

    print("* wrote %d rows -> %s" % (len(rows), args.output))


if __name__ == "__main__":
    main()
