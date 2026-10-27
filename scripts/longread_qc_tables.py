#!/usr/bin/env python
"""Per-library long-read QC: a log2 read-length histogram of the aligned reads (bins of 0.1), and
read counts, aligned fraction and median lengths from NanoStat plus the pychopper full-length share
(Primers_found over all classified reads, as MultiQC defines it)."""
import argparse
import os

import numpy as np
import pandas as pd


def sample_of(path, suffix):
    return os.path.basename(path)[:-len(suffix)]


def nanostat(paths, suffix):
    rows = {}
    for p in paths:
        d = pd.read_csv(p, sep="\t", index_col=0).iloc[:, 0]
        rows[sample_of(p, suffix)] = {"number_of_reads": int(float(d["number_of_reads"])),
                                      "median_read_length": float(d["median_read_length"])}
    return pd.DataFrame.from_dict(rows, orient="index")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--read-lengths", nargs="+", required=True)
    ap.add_argument("--nanostat-raw", nargs="+", required=True)
    ap.add_argument("--nanostat-aligned", nargs="+", required=True)
    ap.add_argument("--pychopper", nargs="+", required=True)
    ap.add_argument("--hist-out", required=True)
    ap.add_argument("--summary-out", required=True)
    args = ap.parse_args()

    edges = np.round(np.arange(0, 18.05, 0.1), 1)
    mids = np.round(edges[:-1] + 0.05, 2)
    hist = []
    for p in sorted(args.read_lengths):
        x = pd.read_csv(p, header=None, dtype=np.int64)[0].to_numpy()
        n, _ = np.histogram(np.clip(np.log2(x), 0, 17.999), bins=edges)
        hist.append(pd.DataFrame({"sample": sample_of(p, ".mapped.readlengths.txt"),
                                  "log2_length_mid": mids, "n_reads": n, "share": n / n.sum()}))
    pd.concat(hist).to_csv(args.hist_out, sep="\t", index=False)

    raw = nanostat(args.nanostat_raw, ".nanostat.raw.tsv")
    aln = nanostat(args.nanostat_aligned, ".nanostat.aligned.tsv")
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str).set_index("sample")
    out = pd.DataFrame({"sample": raw.index, "condition": meta.loc[raw.index, "condition"].to_numpy(),
                        "number_of_reads_raw": raw["number_of_reads"].to_numpy(),
                        "number_of_reads_aligned": aln.loc[raw.index, "number_of_reads"].to_numpy(),
                        "median_read_length_raw": raw["median_read_length"].to_numpy(),
                        "median_read_length_aligned": aln.loc[raw.index, "median_read_length"].to_numpy()})
    out.insert(4, "aligned_frac", out["number_of_reads_aligned"] / out["number_of_reads_raw"])
    fl = {}
    for p in args.pychopper:
        d = pd.read_csv(p, sep="\t")
        c = d[d["Category"] == "Classification"].set_index("Name")["Value"]
        fl[sample_of(p, ".stats.tsv")] = c["Primers_found"] / c[["Primers_found", "Rescue", "Unusable"]].sum()
    out["full_length_share"] = out["sample"].map(fl)
    out.sort_values("sample").to_csv(args.summary_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
