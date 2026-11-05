#!/usr/bin/env python
"""Per-library long-read QC: a log2 read-length histogram of the aligned reads (bins of 0.1), and
read counts, aligned fraction and median lengths from NanoStat plus the pychopper full-length share
(Primers_found over all classified reads, as MultiQC defines it)."""
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
    # read count and median read length per sample from the NanoStat tsv (one value column)
    rows = {}
    for p in paths:
        stats = pd.read_csv(p, sep="\t", index_col=0).iloc[:, 0]
        rows[sample_of(p, suffix)] = {"number_of_reads": int(float(stats["number_of_reads"])),
                                      "median_read_length": float(stats["median_read_length"])}
    return pd.DataFrame.from_dict(rows, orient="index")


def main():
    # read the options
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--read-lengths", nargs="+", required=True)
    ap.add_argument("--nanostat-raw", nargs="+", required=True)
    ap.add_argument("--nanostat-aligned", nargs="+", required=True)
    ap.add_argument("--pychopper", nargs="+", required=True)
    ap.add_argument("--hist-out", required=True)
    ap.add_argument("--summary-out", required=True)
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

    # one summary row per sample: raw and aligned read counts and median lengths from NanoStat
    raw = nanostat(args.nanostat_raw, ".nanostat.raw.tsv")
    aln = nanostat(args.nanostat_aligned, ".nanostat.aligned.tsv")
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str).set_index("sample")
    out = pd.DataFrame({"sample": raw.index, "condition": meta.loc[raw.index, "condition"].to_numpy(),
                        "number_of_reads_raw": raw["number_of_reads"].to_numpy(),
                        "number_of_reads_aligned": aln.loc[raw.index, "number_of_reads"].to_numpy(),
                        "median_read_length_raw": raw["median_read_length"].to_numpy(),
                        "median_read_length_aligned": aln.loc[raw.index, "median_read_length"].to_numpy()})
    # aligned_frac = mapped records of the alignment (supplementary records included) over the raw reads
    # before pychopper, so it folds the pychopper loss into the alignment rate
    out.insert(4, "aligned_frac", out["number_of_reads_aligned"] / out["number_of_reads_raw"])
    # pychopper full-length share = reads with both primers found over all classified reads
    full_length_share = {}
    for p in args.pychopper:
        pychopper_stats = pd.read_csv(p, sep="\t")
        classes = pychopper_stats[pychopper_stats["Category"] == "Classification"].set_index("Name")["Value"]
        classified = classes[["Primers_found", "Rescue", "Unusable"]].sum()
        full_length_share[sample_of(p, ".stats.tsv")] = classes["Primers_found"] / classified
    out["full_length_share"] = out["sample"].map(full_length_share)
    out.sort_values("sample").to_csv(args.summary_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
