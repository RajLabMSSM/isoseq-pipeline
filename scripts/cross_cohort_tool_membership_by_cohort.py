#!/usr/bin/env python3
"""Per cross-cohort novel isoform, one 0/1 column per cohort x tool: 1 when any of its members from
that cohort was built by that tool (from the cohort's novel_provenance.tsv)."""
import argparse

TOOLS = ["bambu", "isoquant", "stringtie"]


def load_cohort(path):
    d = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[idx["transcript_id"]]] = tuple(
                int(f[idx[t]]) if f[idx[t]] not in ("NA", "") else 0 for t in TOOLS)
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cross-provenance", required=True)
    ap.add_argument("--cohort", action="append", required=True, metavar="NAME=PATH")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    order = []
    per_cohort = {}
    for spec in args.cohort:
        name, path = spec.split("=", 1)
        order.append(name)
        per_cohort[name] = load_cohort(path)

    cols = ["new_id", "gene_id", "class_code", "n_cohorts", "cohorts"]
    for c in order:
        cols += ["%s_%s" % (c, t) for t in TOOLS]

    with open(args.cross_provenance) as fh, open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        h = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            flags = {c: [0, 0, 0] for c in order}
            for src in f[idx["src_ids"]].split("|"):
                if ":" not in src:
                    continue
                coh, tid = src.split(":", 1)
                if coh in per_cohort and tid in per_cohort[coh]:
                    tv = per_cohort[coh][tid]
                    flags[coh] = [a | b for a, b in zip(flags[coh], tv)]
            row = [f[idx["new_id"]], f[idx["gene_id"]], f[idx["class_code"]],
                   f[idx["n_cohorts"]], f[idx["cohorts"]]]
            for c in order:
                row += [str(v) for v in flags[c]]
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
