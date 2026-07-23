#!/usr/bin/env python
"""
cross_cohort_global_provenance.py

Global per-isoform provenance for the cross-cohort reference: for each cross-cohort novel
isoform (XNOVEL_*), report which COHORTS contributed AND the POOLED tool support (union of
bambu/isoquant/stringtie across every contributing cohort's isoform). Joins the
cross-cohort union's provenance (src_ids column) to each cohort's per-isoform
_novel_provenance.tsv.
"""
import argparse
from collections import defaultdict

TOOLS = ["bambu", "isoquant", "stringtie"]


def load_cohort_prov(path):
    """NOVEL id -> {tool: 0/1}."""
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[0]] = {t: f[idx[t]] for t in TOOLS if t in idx}
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cross-provenance", required=True,
                    help="cross_cohort_novel_provenance.tsv (has src_ids)")
    ap.add_argument("--cohort-provenance", action="append", required=True, metavar="NAME=PATH",
                    help="per-cohort _novel_provenance.tsv; repeat per cohort")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    cohort_prov = {}
    for spec in args.cohort_provenance:
        name, path = spec.split("=", 1)
        cohort_prov[name] = load_cohort_prov(path)

    with open(args.cross_provenance) as fh, open(args.output, "w") as out:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(hdr)}
        out.write("new_id\tgene_id\tclass_code\tn_cohorts\tcohorts\tn_tools_pooled\t%s\n"
                  % "\t".join(TOOLS))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            nid = f[idx["new_id"]]
            gid = f[idx["gene_id"]]
            cc = f[idx["class_code"]]
            cohorts = f[idx["cohorts"]]
            src_ids = f[idx["src_ids"]]
            pooled = {t: 0 for t in TOOLS}
            for pair in src_ids.split("|"):
                if ":" not in pair:
                    continue
                coh, orig = pair.split(":", 1)
                row = cohort_prov.get(coh, {}).get(orig)
                if not row:
                    continue
                for t in TOOLS:
                    if row.get(t) == "1":
                        pooled[t] = 1
            flags = "\t".join(str(pooled[t]) for t in TOOLS)
            out.write("%s\t%s\t%s\t%s\t%s\t%d\t%s\n"
                      % (nid, gid, cc, f[idx["n_cohorts"]], cohorts, sum(pooled.values()), flags))
    print("* wrote global provenance -> %s" % args.output)


if __name__ == "__main__":
    main()
