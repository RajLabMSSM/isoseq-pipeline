#!/usr/bin/env python3
"""Record which assembly tools support each cross-cohort novel isoform in each cohort.

The output has one 0/1 column for every combination of cohort and tool.
A column is 1 when any of the isoform's members from that cohort was built by that tool.
The tool calls come from each cohort's novel_provenance.tsv."""
import argparse

TOOLS = ["bambu", "isoquant", "stringtie"]


def load_cohort(path):
    # Read one cohort's provenance table into a dictionary of transcript id to 0/1 flags for bambu,
    # isoquant and stringtie. NA or an empty value counts as 0.
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

    # Load every cohort's tool flags, keeping the order in which the cohorts were given.
    order = []
    per_cohort = {}
    for spec in args.cohort:
        name, path = spec.split("=", 1)
        order.append(name)
        per_cohort[name] = load_cohort(path)

    # The output starts with the cross-cohort columns, then one <cohort>_<tool> column for every cohort and tool.
    cols = ["new_id", "gene_id", "class_code", "n_cohorts", "cohorts"]
    for c in order:
        cols += ["%s_%s" % (c, t) for t in TOOLS]

    with open(args.cross_provenance) as fh, open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # Set a cohort's tool flag to 1 when any member from that cohort was built by that tool.
            flags = {c: [0] * len(TOOLS) for c in order}
            for src in f[idx["src_ids"]].split("|"):
                if ":" not in src:
                    continue
                cohort, tid = src.split(":", 1)
                if cohort in per_cohort and tid in per_cohort[cohort]:
                    member_flags = per_cohort[cohort][tid]
                    flags[cohort] = [a | b for a, b in zip(flags[cohort], member_flags)]
            row = [f[idx["new_id"]], f[idx["gene_id"]], f[idx["class_code"]],
                   f[idx["n_cohorts"]], f[idx["cohorts"]]]
            for c in order:
                row += [str(v) for v in flags[c]]
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
