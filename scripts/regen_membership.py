#!/usr/bin/env python
"""
regen_membership.py

Regenerate the per-TCONS tool-membership TSV (bambu/isoquant/stringtie support per novel
consensus isoform) keyed on the REAL merged_consensus TCONS ids, WITHOUT re-running the
assemblers. Used to backfill provenance for a run whose union/tool_membership were temp()
and deleted. Going forward tool_membership is no longer temp() (consensus_pipeline.smk),
so this is only needed for the pre-fix v5 references.

Method: the surviving gffcmp.annotated.gtf IS the gffcompare-annotated union -- every
assembler transcript, with class_code inline and tool encoded in the transcript_id prefix
(<tool>__<group>__<id>). Rebuild the union in memory, re-run the SAME end-aware clustering
(imports consensus_endaware.cluster), take each cluster's tool set, then assign it to the
real TCONS whose representative matches the cluster representative within tolerance. Any
real TCONS with no matching cluster (a handful of class-code edge cases that can't be
reproduced from the annotated gtf alone) is written NA and reported. Validate the result
against the PRESERVED <dc>_tool_overlap.tsv Venn (aggregate is exact).
"""
import argparse
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consensus_endaware import cluster, introns_match, choose_rep, TOOLS  # noqa: E402


def attr(f, k):
    m = re.search(k + r' "([^"]*)"', f)
    return m.group(1) if m else None


def load(path, tool_from_prefix):
    exons = defaultdict(list)
    tmeta = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("transcript", "exon"):
                continue
            tid = attr(c[8], "transcript_id")
            if tid is None:
                continue
            if c[2] == "transcript":
                tmeta[tid] = (c[0], c[6], attr(c[8], "class_code"))
            else:
                exons[tid].append((int(c[3]), int(c[4])))
    tx = {}
    for tid, (chrom, strand, cc) in tmeta.items():
        ex = sorted(exons[tid])
        intr = tuple(x for i in range(len(ex) - 1) for x in (ex[i][1], ex[i + 1][0]))
        lo, hi = ex[0][0], ex[-1][1]
        tss, tes = (hi, lo) if strand == "-" else (lo, hi)
        tx[tid] = {"chrom": chrom, "strand": strand, "introns": intr, "tss": tss, "tes": tes,
                   "span": (lo, hi), "len": sum(e - s + 1 for s, e in ex), "cc": cc,
                   "tool": tid.split("__")[0] if tool_from_prefix else None}
    return tx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", required=True, help="<dc>_gffcmp.annotated.gtf (the union)")
    ap.add_argument("--merged", required=True, help="<dc>_merged_consensus.gtf (real TCONS)")
    ap.add_argument("--exclude-codes", default="=,c")
    ap.add_argument("--wobble", type=int, default=6)
    ap.add_argument("--tss-tol", type=int, default=100)
    ap.add_argument("--tes-tol", type=int, default=100)
    ap.add_argument("--ignore-tss-multiexon", action="store_true",
                    help="must MATCH the consensus run or membership will not line up")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    drop = set(c for c in args.exclude_codes.split(",") if c)

    uni = load(args.annotated, True)
    novel = [u for u, d in uni.items() if d["cc"] not in drop]
    clusters = cluster(novel, uni, args.wobble, args.tss_tol, args.tes_tol,
                       args.ignore_tss_multiexon)

    # cluster representatives (bucketed) -> tool set
    reps = defaultdict(list)   # (chrom,strand,nintr) -> [(rep_dict, tools)]
    for members in clusters:
        rep = uni[choose_rep(members, uni)]
        tools = set(uni[u]["tool"] for u in members)
        reps[(rep["chrom"], rep["strand"], len(rep["introns"]))].append((rep, tools))

    cons = load(args.merged, False)
    combo = defaultdict(int)
    na = 0
    with open(args.output, "w") as o:
        o.write("transcript_id\tclass_code\tis_novel\tn_tools\t%s\n" % "\t".join(TOOLS))
        for tid, d in cons.items():
            best, bestd = None, None
            for rep, tools in reps.get((d["chrom"], d["strand"], len(d["introns"])), []):
                skip_tss = args.ignore_tss_multiexon and len(d["introns"]) > 0
                if ((skip_tss or abs(d["tss"] - rep["tss"]) <= args.tss_tol)
                        and abs(d["tes"] - rep["tes"]) <= args.tes_tol
                        and introns_match(d["introns"], rep["introns"], args.wobble)):
                    dist = abs(d["tss"] - rep["tss"]) + abs(d["tes"] - rep["tes"])
                    if bestd is None or dist < bestd:
                        best, bestd = tools, dist
            if best is None:
                o.write("%s\t%s\t1\tNA\t%s\n" % (tid, d["cc"], "\t".join(["NA"] * len(TOOLS))))
                na += 1
                continue
            flags = "\t".join("1" if t in best else "0" for t in TOOLS)
            o.write("%s\t%s\t1\t%d\t%s\n" % (tid, d["cc"], len(best), flags))
            combo["+".join(t for t in TOOLS if t in best)] += 1

    print("* regenerated membership for %d TCONS (%d NA) -> %s" % (len(cons), na, args.output))
    print("* Venn by tool combination (compare to <dc>_tool_overlap.tsv):")
    for k in sorted(combo, key=lambda x: (x.count("+"), x)):
        print("    %-28s %d" % (k, combo[k]))
    print("    total_matched %d" % sum(combo.values()))


if __name__ == "__main__":
    main()
