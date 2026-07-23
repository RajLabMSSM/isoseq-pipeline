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

# EVERY metric NanoStat --tsv emits is captured -- nothing is dropped. The compound fields
# are split into their components so each is a real column:
#   "Reads >Q10:"                  "7625392 (93.3%) 6118.6Mb" -> n_ / pct_ / mb_reads_above_q10
#   "longest_read_(with_Q):1"      "390790 (20.6)"            -> longest_read_1_len / _qual
#   "highest_Q_read_(with_length):1" "24.3 (5186)"            -> highest_q_read_1_qual / _len
SCALARS = [
    "number_of_reads",
    "number_of_bases",
    "median_read_length",
    "mean_read_length",
    "read_length_stdev",
    "n50",
    "mean_qual",
    "median_qual",
    # aligned stage only (NanoStat --bam); "NA" for raw
    "number_of_bases_aligned",
    "fraction_bases_aligned",
    "average_identity",
    "median_identity",
]

QCUTS = [5, 7, 10, 12, 15]
TOPN = [1, 2, 3, 4, 5]

QCOLS = []
for _q in QCUTS:
    QCOLS += ["n_reads_above_q%d" % _q,
              "pct_reads_above_q%d" % _q,
              "mb_reads_above_q%d" % _q]

TOPCOLS = []
for _i in TOPN:
    TOPCOLS += ["longest_read_%d_len" % _i, "longest_read_%d_qual" % _i]
for _i in TOPN:
    TOPCOLS += ["highest_q_read_%d_qual" % _i, "highest_q_read_%d_len" % _i]

ALLCOLS = SCALARS + QCOLS + TOPCOLS


def _pair(val):
    """'390790 (20.6)' -> ('390790', '20.6'); tolerant of missing parens."""
    val = val.strip()
    if "(" not in val:
        return val, "NA"
    a, b = val.split("(", 1)
    return a.strip(), b.split(")")[0].strip()


def parse_nanostat(path):
    """Return {column: value} for EVERY metric in a NanoStat --tsv file."""
    vals = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            key, val = parts
            if key in SCALARS:
                vals[key] = val
                continue

            flat = key.replace(" ", "")
            matched = False
            for q in QCUTS:
                if flat == "Reads>Q%d:" % q:
                    # "7625392 (93.3%) 6118.6Mb"
                    n = val.split("(")[0].strip() if "(" in val else "NA"
                    pct = val.split("(")[1].split("%")[0] if "(" in val else "NA"
                    mb = val.split(")")[-1].replace("Mb", "").strip() if ")" in val else "NA"
                    vals["n_reads_above_q%d" % q] = n
                    vals["pct_reads_above_q%d" % q] = pct
                    vals["mb_reads_above_q%d" % q] = mb or "NA"
                    matched = True
                    break
            if matched:
                continue

            for i in TOPN:
                if flat == "longest_read_(with_Q):%d" % i:
                    ln, q = _pair(val)
                    vals["longest_read_%d_len" % i] = ln
                    vals["longest_read_%d_qual" % i] = q
                elif flat == "highest_Q_read_(with_length):%d" % i:
                    qv, ln = _pair(val)
                    vals["highest_q_read_%d_qual" % i] = qv
                    vals["highest_q_read_%d_len" % i] = ln
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

    header = ["sample", "stage"] + ALLCOLS
    with open(args.output, "w") as out:
        out.write("\t".join(header) + "\n")
        for sample, stage, vals in rows:
            out.write("\t".join([sample, stage]
                                + [vals.get(k, "NA") for k in ALLCOLS]) + "\n")

    print("* wrote %d rows -> %s" % (len(rows), args.output))


if __name__ == "__main__":
    main()
