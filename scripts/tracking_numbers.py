#!/usr/bin/env python
"""Funnel of the pooled isoform build as metric<TAB>value: raw per tool -> novel candidates -> single-exon and
unstranded dropped -> tama_merge -> tool overlap -> SQANTI filter -> final novel.

Without --only the funnel covers every cohort. With --only <cohort> it covers that cohort's share: its own
assembler counts, and the merged models that at least one of its assemblies built, with the tools of that cohort
only. A model built in both cohorts therefore counts in both cohort funnels.

raw_*_all counts the union after the transcripts that carry a GENCODE id were removed. collapsed_removed = spliced
candidates - consensus models. lr_specific / sr_supported split the consensus models by SQANTI's min_cov (0 = a
junction no short read covers). A model built by two tools counts under both."""
import argparse
from collections import Counter
from itertools import combinations

TOOLS = ["bambu", "isoquant", "stringtie"]
# Every tool combination, the three-tool one first, then pairs, then single tools.
COMBOS = [c for k in range(len(TOOLS), 0, -1) for c in combinations(TOOLS, k)]


def read_table(path):
    # Yield each row of a tab-separated table with a header as a dict.
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(hdr, line.rstrip("\n").split("\t")))


def overlap(tools_of, keep_ids):
    # Among the models in keep_ids, count the models per tool combination (venn) and per tool. A model built by two
    # tools counts under both tools.
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
    ap.add_argument("--cohort", required=True, help="label written in the first row")
    ap.add_argument("--only", help="restrict the funnel to the models built by this cohort")
    ap.add_argument("--metric-prefix", default="", help="put in front of every metric name")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    # Sum the per-cohort, per-tool counts of the consensus stats table, over all cohorts or over the one given.
    raw, novelcand, mono, unstr = Counter(), Counter(), Counter(), Counter()
    for r in read_table(args.stats):
        if args.only and r["cohort"] != args.only:
            continue
        raw[r["tool"]] += int(r["raw_all"])
        novelcand[r["tool"]] += int(r["novel_candidates"])
        mono[r["tool"]] += int(r["monoexon_dropped"])
        unstr[r["tool"]] += int(r["unstranded_dropped"])

    # Map each consensus model in scope to the tools that built it. Under --only the tools are read from that
    # cohort's <cohort>_<tool> columns, so the other cohort's tools do not count. The cohorts of every model are
    # kept for the cohort combination counts.
    tools = {}
    cohorts_of = {}
    for r in read_table(args.membership):
        cohorts = r["cohorts"].split(",")
        if args.only and args.only not in cohorts:
            continue
        if args.only:
            tools[r["transcript_id"]] = tuple(t for t in TOOLS if r["%s_%s" % (args.only, t)] == "1")
        else:
            tools[r["transcript_id"]] = tuple(t for t in TOOLS if r[t] == "1")
        cohorts_of[r["transcript_id"]] = r["cohorts"]

    # Read the SQANTI result: the models in scope that passed its filter, and the lowest short-read coverage of any
    # of their junctions (min_cov).
    passed = set()
    min_cov = {}
    for r in read_table(args.sqanti_classification):
        if r["isoform"] not in tools:
            continue
        if r["filter_result"] == "Isoform":
            passed.add(r["isoform"])
        if r["min_cov"] != "NA":
            min_cov[r["isoform"]] = int(r["min_cov"])

    # Read the final novel ids in scope, and count the reference transcripts (transcript lines of the reference).
    with open(args.pass_ids) as fh:
        final_ids = {line.strip() for line in fh if line.strip()} & tools.keys()
    n_ref = 0
    with open(args.reference) as fh:
        for line in fh:
            fields = line.split("\t")
            if not line.startswith("#") and len(fields) > 2 and fields[2] == "transcript":
                n_ref += 1

    # Count the tool overlap of all consensus models in scope and of the SQANTI-passed ones.
    n_cons = len(tools)
    sqanti_pass = len(passed)
    venn, cons_per_tool = overlap(tools, set(tools))
    venn_sq, sq_per_tool = overlap(tools, final_ids)
    # Split the models into long-read specific ones (a junction with no short-read coverage) and short-read
    # supported ones.
    lr_ids = {t for t, v in min_cov.items() if v == 0}
    sr_ids = {t for t, v in min_cov.items() if v > 0}
    _, lr = overlap(tools, lr_ids)
    _, sr = overlap(tools, sr_ids)
    # The spliced, stranded candidates that went into tama_merge.
    spliced = {t: novelcand[t] - mono[t] - unstr[t] for t in TOOLS}

    # Build the metric rows in funnel order. Each per-tool block ends with its total.
    rows = []
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
    rows.append(("lr_specific_total", len(lr_ids)))
    rows += [("sr_supported_%s" % t, sr[t]) for t in TOOLS]
    rows.append(("sr_supported_total", len(sr_ids)))
    # novel_pre_sqanti and sqanti_removed count every merged model, including the few that the merged gffcompare
    # re-classified = or c, which are left out of the SQANTI input.
    rows += [("novel_pre_sqanti", n_cons), ("sqanti_pass", sqanti_pass),
             ("sqanti_removed", n_cons - sqanti_pass), ("final_novel", len(final_ids))]
    # The final novel models per cohort combination, as in the cohorts attribute of the final GTF. These rows are
    # written for the whole build only.
    if not args.only:
        by_combo = Counter(cohorts_of[t] for t in final_ids)
        rows += [("by_combo:%s" % combo, n) for combo, n in sorted(by_combo.items(), key=lambda kv: (-kv[1], kv[0]))]
        shared = sum(n for combo, n in by_combo.items() if "," in combo)
        rows += [("single_cohort_only", len(final_ids) - shared), ("shared_by_2plus_cohorts", shared)]
    rows += [("reference_transcripts", n_ref), ("final_total", n_ref + len(final_ids))]

    with open(args.output, "w") as out:
        out.write("metric\tvalue\n")
        out.write("cohort\t%s\n" % args.cohort)
        out.writelines("%s%s\t%s\n" % (args.metric_prefix, m, v) for m, v in rows)
    print("* tracking_numbers -> %s (%d metrics)" % (args.output, len(rows) + 1))


if __name__ == "__main__":
    main()
