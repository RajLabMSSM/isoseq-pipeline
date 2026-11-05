#!/usr/bin/env python
"""Per-cohort funnel as metric<TAB>value: raw per tool -> novel candidates -> single-exon and
unstranded dropped -> tama_merge -> tool overlap -> SQANTI filter -> final novel.

raw_*_all includes the reference passthrough of bambu and isoquant. collapsed_removed = spliced
candidates - consensus models. lr_specific / sr_supported split the consensus models by SQANTI's
min_cov (0 = a junction no short read covers); a model built by two tools counts under both."""
import argparse
from collections import Counter
from itertools import combinations

TOOLS = ["bambu", "isoquant", "stringtie"]
# every tool combination, the three-tool one first, then pairs, then single tools
COMBOS = [c for k in range(len(TOOLS), 0, -1) for c in combinations(TOOLS, k)]


def read_table(path):
    # yield each row of a tab-separated table with a header as a dict
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(hdr, line.rstrip("\n").split("\t")))


def overlap(tools_of, keep_ids):
    # among the models in keep_ids: count per tool combination (venn) and per tool (a model built by two
    # tools counts under both)
    venn = Counter(present for tid, present in tools_of.items() if tid in keep_ids)
    per_tool = Counter()
    for present, n in venn.items():
        for t in present:
            per_tool[t] += n
    return venn, per_tool


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--membership", required=True)
    ap.add_argument("--sqanti-classification", required=True)
    ap.add_argument("--pass-ids", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    # per-tool counts from the consensus stats table
    stats = {r["tool"]: r for r in read_table(args.stats)}
    raw = {t: int(stats[t]["raw_all"]) for t in TOOLS}
    novelcand = {t: int(stats[t]["novel_candidates"]) for t in TOOLS}
    mono = {t: int(stats[t]["monoexon_dropped"]) for t in TOOLS}
    unstr = {t: int(stats[t]["unstranded_dropped"]) for t in TOOLS}
    # consensus model -> tools that built it, from the membership table
    tools = {r["transcript_id"]: tuple(t for t in TOOLS if r[t] == "1") for r in read_table(args.membership)}
    # SQANTI: models passing its filter, and the lowest short-read coverage of any junction (min_cov)
    passed = set()
    min_cov = {}
    for r in read_table(args.sqanti_classification):
        if r["filter_result"] == "Isoform":
            passed.add(r["isoform"])
        if r["min_cov"] != "NA":
            min_cov[r["isoform"]] = int(r["min_cov"])
    # final novel ids, and the number of reference transcripts (transcript lines of the reference gtf)
    with open(args.pass_ids) as fh:
        final_ids = {line.strip() for line in fh if line.strip()}
    n_ref = 0
    with open(args.reference) as fh:
        for line in fh:
            fields = line.split("\t")
            if not line.startswith("#") and len(fields) > 2 and fields[2] == "transcript":
                n_ref += 1

    # tool overlap of all consensus models and of the SQANTI-passed ones
    n_cons = len(tools)
    sqanti_pass = len(passed)
    venn, cons_per_tool = overlap(tools, set(tools))
    venn_sq, sq_per_tool = overlap(tools, final_ids)
    # long-read specific (a junction with no short-read coverage) versus short-read supported models
    lr_ids = {t for t, v in min_cov.items() if v == 0}
    sr_ids = {t for t, v in min_cov.items() if v > 0}
    _, lr = overlap(tools, lr_ids)
    _, sr = overlap(tools, sr_ids)
    # spliced, stranded candidates that went into tama_merge
    spliced = {t: novelcand[t] - mono[t] - unstr[t] for t in TOOLS}

    # metric rows in funnel order; each per-tool block ends with its total
    rows = [("cohort", args.cohort)]
    for name, d in (("raw_%s_all", raw), ("novel_candidates_%s", novelcand),
                    ("novel_monoexon_dropped_%s", mono), ("novel_unstranded_dropped_%s", unstr),
                    ("novel_spliced_candidates_%s", spliced)):
        rows += [(name % t, d[t]) for t in TOOLS]
        rows.append((name.replace("%s", "total"), sum(d.values())))
    rows += [("consensus_novel_total", n_cons),
             ("collapsed_removed", sum(spliced.values()) - n_cons)]
    rows += [("overlap_%s" % "+".join(c), venn[c]) for c in COMBOS]
    rows += [("consensus_supported_by_%s" % t, cons_per_tool[t]) for t in TOOLS]
    rows += [("sqanti_overlap_%s" % "+".join(c), venn_sq[c]) for c in COMBOS]
    rows += [("sqanti_pass_%s" % t, sq_per_tool[t]) for t in TOOLS]
    rows += [("sqanti_pass_%dtool" % k, sum(v for c, v in venn_sq.items() if len(c) == k)) for k in (1, 2, 3)]
    rows += [("lr_specific_%s" % t, lr[t]) for t in TOOLS]
    rows.append(("lr_specific_total", len(lr_ids & tools.keys())))
    rows += [("sr_supported_%s" % t, sr[t]) for t in TOOLS]
    rows.append(("sr_supported_total", len(sr_ids & tools.keys())))
    # novel_pre_sqanti and sqanti_removed count every merged model, including the few that the merged
    # gffcompare re-classified = or c, which are left out of the SQANTI input
    rows += [("novel_pre_sqanti", n_cons), ("sqanti_pass", sqanti_pass),
             ("sqanti_removed", n_cons - sqanti_pass), ("final_novel", len(final_ids)),
             ("reference_transcripts", n_ref), ("final_total", n_ref + len(final_ids))]

    with open(args.output, "w") as out:
        out.write("metric\tvalue\n")
        out.writelines("%s\t%s\n" % r for r in rows)
    print("* tracking_numbers -> %s (%d metrics)" % (args.output, len(rows)))


if __name__ == "__main__":
    main()
