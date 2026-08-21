#!/usr/bin/env python
"""
lapa_diff.py

Two-group differential alternative polyadenylation (APA) from a LAPA result dir.

LAPA 0.0.5's built-in `beta_binomial_test()` scores each poly(A) site per-sample
against a pooled baseline -- it is NOT a group-A-vs-group-B contrast (LAPA's
two-condition differential is marked "coming soon" in its docs). So we run the
standard APA two-group test on LAPA's OWN per-sample count matrix
(`LapaResult.counts()`): for each gene with >=2 poly(A) sites, compare each site's
usage (PAU = site reads / gene reads) between the two groups with Fisher's exact test
on pooled counts, report delta-PAU, and BH-correct.

Must run in the lapa conda env (imports lapa). Reads the sample->group mapping from
the LAPA sample-config CSV (sample,dataset,path); `dataset` == group.

Output columns: polya_site, Chromosome, Start, End, Strand, gene_id,
  count_<A>, count_<B>, PAU_<A>, PAU_<B>, delta_PAU (B-A), pval, padj, n_sites_gene.
"""
import argparse
import csv

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from lapa.result import LapaResult


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def group_samples(config_csv, group_a, group_b):
    a, b = [], []
    with open(config_csv) as fh:
        for row in csv.DictReader(fh):
            if row["dataset"] == group_a:
                a.append(row["sample"])
            elif row["dataset"] == group_b:
                b.append(row["sample"])
    return a, b


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lapa-dir", required=True, help="LAPA output dir (has polyA_clusters.bed etc.)")
    ap.add_argument("--config", required=True, help="LAPA sample-config CSV (sample,dataset,path)")
    ap.add_argument("--group-a", required=True)
    ap.add_argument("--group-b", required=True)
    ap.add_argument("--min-gene-count", type=int, default=10,
                    help="min pooled gene read count in EACH group to test (default 10)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    samples_a, samples_b = group_samples(args.config, args.group_a, args.group_b)
    if not samples_a or not samples_b:
        raise SystemExit("no samples for one of the groups (%s / %s)" % (args.group_a, args.group_b))

    res = LapaResult(args.lapa_dir)
    df = res.counts()                      # per poly(A) site: coords + gene_id + per-sample counts
    df = df.reset_index()
    sa = [s for s in samples_a if s in df.columns]
    sb = [s for s in samples_b if s in df.columns]
    df["count_a"] = df[sa].sum(axis=1)
    df["count_b"] = df[sb].sum(axis=1)

    gene = df.groupby("gene_id")
    df["gene_a"] = gene["count_a"].transform("sum")
    df["gene_b"] = gene["count_b"].transform("sum")
    df["n_sites_gene"] = gene["gene_id"].transform("size")

    test = df[(df["n_sites_gene"] >= 2)
              & (df["gene_a"] >= args.min_gene_count)
              & (df["gene_b"] >= args.min_gene_count)].copy()

    pvals = []
    for _, r in test.iterrows():
        a_in, b_in = r["count_a"], r["count_b"]
        a_out, b_out = r["gene_a"] - a_in, r["gene_b"] - b_in
        _, p = fisher_exact([[a_in, a_out], [b_in, b_out]])
        pvals.append(p)
    test["PAU_a"] = test["count_a"] / test["gene_a"]
    test["PAU_b"] = test["count_b"] / test["gene_b"]
    test["delta_PAU"] = test["PAU_b"] - test["PAU_a"]
    test["pval"] = pvals
    test["padj"] = bh_adjust(pvals) if pvals else []

    keep = ["Chromosome", "Start", "End", "Strand", "gene_id",
            "count_a", "count_b", "PAU_a", "PAU_b", "delta_PAU", "pval", "padj", "n_sites_gene"]
    keep = [c for c in keep if c in test.columns]
    out = test[keep].rename(columns={
        "count_a": "count_%s" % args.group_a, "count_b": "count_%s" % args.group_b,
        "PAU_a": "PAU_%s" % args.group_a, "PAU_b": "PAU_%s" % args.group_b,
        "delta_PAU": "delta_PAU_%s_minus_%s" % (args.group_b, args.group_a)})
    out = out.sort_values("padj")
    out.to_csv(args.output, sep="\t", index=False)
    n_sig = int((test["padj"] < 0.05).sum()) if len(test) else 0
    print("* differential APA: %d sites tested (%s n=%d vs %s n=%d), %d sig (padj<0.05) -> %s"
          % (len(test), args.group_a, len(sa), args.group_b, len(sb), n_sig, args.output))


if __name__ == "__main__":
    main()
