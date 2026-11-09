#!/usr/bin/env python
"""Per-library long-read QC at three stages (raw reads, after pychopper, after minimap2): a log2 read-length histogram of
the primary mapped reads (bins of 0.1), a one-row-per-library summary of read counts, median lengths and yields, and
every NanoStat metric of every stage kept in its own table."""

# In simpler terms, the script answers four questions:
# 1. How long are the aligned Nanopore reads?
# 2. How many reads did each sample have before pychopper, after pychopper and after minimap2?
# 3. How much was lost at pychopper, and how much at the alignment?
# 4. What fraction of raw reads did Pychopper find primers on at both ends?

# outputs:
# read_length_hist.tsv.gz =	Histogram data showing the distribution of primary mapped read lengths for every sample
# longread_qc_summary.tsv =	One row per sample containing the read counts, median lengths and yields that figures draw from
# nanostat_all_metrics.tsv = One row per sample per stage holding every NanoStat metric, kept to look back at, not plotted


import argparse
import os

import numpy as np
import pandas as pd

# histogram on the log2 read length: bins of 0.1 from 0 to 18 (2^18 = 262 kb, past any cDNA read)
LOG2_BIN_WIDTH = 0.1
LOG2_MAX = 18.0
# longer reads are clipped just under the top edge so they land in the last bin
LOG2_CLIP_TOP = 17.999


def sample_of(path, suffix):
    # sample id = file name minus its known suffix
    return os.path.basename(path)[:-len(suffix)]


def nanostat(paths, suffix):
    # every metric NanoStat reports, one row per sample. values stay as written, because some are text like "189507 (18.2)"
    rows = {}
    for p in paths:
        rows[sample_of(p, suffix)] = pd.read_csv(p, sep="\t", index_col=0).iloc[:, 0]
    return pd.DataFrame(rows).T


def pychopper_counts(paths):
    # the counts each pychopper stats table holds, one row per sample. a read with several usable segments is "rescued":
    # pychopper's Rescue value counts its segments, and the RescueSegmentNr rows count the reads themselves
    rows = {}
    for p in paths:
        stats = pd.read_csv(p, sep="\t")
        classes = stats[stats["Category"] == "Classification"].set_index("Name")["Value"]
        read_stats = stats[stats["Category"] == "ReadStats"].set_index("Name")["Value"]
        rows[sample_of(p, ".stats.tsv")] = {
            "pychopper_primers_found": int(classes["Primers_found"]),
            "pychopper_rescued_reads": int(stats.loc[stats["Category"] == "RescueSegmentNr", "Value"].sum()),
            "pychopper_rescued_segments": int(classes["Rescue"]),
            "pychopper_unusable": int(classes["Unusable"]),
            "pychopper_length_fail": int(read_stats["LenFail"]),
            "pychopper_quality_fail": int(read_stats["QcFail"])}
    return pd.DataFrame.from_dict(rows, orient="index")


def main():
    # read the options
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True)                        # sample metadata TSV
    ap.add_argument("--read-lengths", nargs="+", required=True)         # individual primary mapped-read lengths
    ap.add_argument("--nanostat-raw", nargs="+", required=True)         # nanostat reports of the raw reads
    ap.add_argument("--nanostat-pychopper", nargs="+", required=True)   # nanostat reports of the reads pychopper kept
    ap.add_argument("--nanostat-aligned", nargs="+", required=True)     # nanostat reports of the primary mapped reads
    ap.add_argument("--pychopper", nargs="+", required=True)            # pychopper classification stats
    ap.add_argument("--hist-out", required=True)                        # output is a histogram TSV
    ap.add_argument("--summary-out", required=True)                     # output is a QC summary TSV
    ap.add_argument("--nanostat-out", required=True)                    # output is the every-metric NanoStat TSV
    args = ap.parse_args()

    # bin edges and bin midpoints of the log2 length histogram
    edges = np.round(np.arange(0, LOG2_MAX + LOG2_BIN_WIDTH / 2, LOG2_BIN_WIDTH), 1)
    mids = np.round(edges[:-1] + LOG2_BIN_WIDTH / 2, 2)
    # one histogram per sample, as counts and as the share of that sample's reads
    hist = []
    for p in sorted(args.read_lengths):
        lengths = pd.read_csv(p, header=None, dtype=np.int64)[0].to_numpy()
        counts, _ = np.histogram(np.clip(np.log2(lengths), 0, LOG2_CLIP_TOP), bins=edges)
        hist.append(pd.DataFrame({"sample": sample_of(p, ".mapped.readlengths.txt"),
                                  "log2_length_mid": mids, "n_reads": counts, "share": counts / counts.sum()}))
    pd.concat(hist).to_csv(args.hist_out, sep="\t", index=False)

    # read every NanoStat metric at the three stages
    raw = nanostat(args.nanostat_raw, ".nanostat.raw.tsv")
    chopped = nanostat(args.nanostat_pychopper, ".nanostat.pychopper.tsv")
    aln = nanostat(args.nanostat_aligned, ".nanostat.aligned.tsv")

    # keep all of them in one table, one row per sample per stage, so nothing NanoStat measured is thrown away
    all_metrics = pd.concat([table.rename_axis("sample").reset_index().assign(stage=stage)
                             for stage, table in [("raw", raw), ("after_pychopper", chopped), ("after_minimap2", aln)]])
    first_columns = ["sample", "stage"]
    all_metrics = all_metrics[first_columns + [c for c in all_metrics.columns if c not in first_columns]]
    all_metrics.sort_values(first_columns).to_csv(args.nanostat_out, sep="\t", index=False)

    # one summary row per sample: read count and median read length at each stage, the columns figures draw from
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str).set_index("sample")
    out = pd.DataFrame({"sample": raw.index, "condition": meta.loc[raw.index, "condition"].to_numpy()})
    for stage, table in [("raw", raw), ("pychopper", chopped), ("aligned", aln)]:
        out["number_of_reads_" + stage] = pd.to_numeric(table.loc[raw.index, "number_of_reads"]).astype(int).to_numpy()
    for stage, table in [("raw", raw), ("pychopper", chopped), ("aligned", aln)]:
        out["median_read_length_" + stage] = pd.to_numeric(table.loc[raw.index, "median_read_length"]).to_numpy()

    # pychopper_retention = sequences pychopper wrote out over raw reads, the yield of that step. a fused read can yield
    # more than one sequence, so this counts sequences per raw read and could in principle pass 1
    out["pychopper_retention"] = out["number_of_reads_pychopper"] / out["number_of_reads_raw"]
    # mapping_rate = primary mapped reads over the sequences given to minimap2, the yield of the alignment alone
    out["mapping_rate"] = out["number_of_reads_aligned"] / out["number_of_reads_pychopper"]
    # aligned_frac = primary mapped reads (one record per read, no supplementary) over the raw reads
    # before pychopper, so it folds the pychopper loss into the alignment rate
    out["aligned_frac"] = out["number_of_reads_aligned"] / out["number_of_reads_raw"]

    # the pychopper classification counts, kept so each share below can be checked by hand
    counts = pychopper_counts(args.pychopper)
    out = out.join(counts, on="sample")
    # full_length_share = reads with primers at both ends over all raw reads. read counts on both sides: pychopper's own
    # percentage divides by a total that counts rescued reads once per segment and leaves out reads failing its quality
    # filter. it is measured before the 50 nt length filter, so it can exceed pychopper_retention
    out["full_length_share"] = out["pychopper_primers_found"] / out["number_of_reads_raw"]
    out.sort_values("sample").to_csv(args.summary_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
