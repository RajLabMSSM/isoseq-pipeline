#!/usr/bin/env python3
"""Per cross-cohort novel isoform, one 0/1 column per cohort x tool: 1 when any of its members from
that cohort was built by that tool (from the cohort's novel_provenance.tsv)."""
import argparse

TOOLS = ["bambu", "isoquant", "stringtie"]


def load_cohort(path):
    # cohort provenance -> {transcript id: (bambu, isoquant, stringtie) 0/1}; NA or empty counts as 0
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

    # every cohort's tool flags, in the order the cohorts were given
    order = []
    per_cohort = {}
    for spec in args.cohort:
        name, path = spec.split("=", 1)
        order.append(name)
        per_cohort[name] = load_cohort(path)

    # output columns: the cross-cohort ids, then <cohort>_<tool> for every cohort and tool
    cols = ["new_id", "gene_id", "class_code", "n_cohorts", "cohorts"]
    for c in order:
        cols += ["%s_%s" % (c, t) for t in TOOLS]

    with open(args.cross_provenance) as fh, open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # a tool flag of a cohort is 1 when any member from that cohort was built by the tool
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
