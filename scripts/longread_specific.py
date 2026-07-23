#!/usr/bin/env python
"""
longread_specific.py

Per-tool count of LONG-READ-SPECIFIC novel isoforms: novel candidates whose splice
structure could NOT be reconstructed from the short reads. An isoform is long-read-
specific when at least one of its introns is absent from the pooled short-read STAR
junction set (SJ.out.tab). Reported per assembler (bambu / isoquant / stringtie), the
same tool split tracking_numbers already uses for novel_candidates_<tool>.

Source of per-tool isoforms = the consensus gffcompare-annotated GTF
(<dc>_gffcmp.annotated.gtf): every input transcript carries its tool in the id prefix
(<tool>__<group>__<origid>) and a gffcompare class_code, so "novel candidate" here is
class_code not in --exclude-codes, matching consensus_stats/tracking_numbers exactly.

Short-read junction universe = every STAR SJ.out.tab in --sj-folder (col2 = intron
start, col3 = intron end, col7 = uniquely-mapped read count), pooled across samples with
total uniq >= --sj-min-uniq, optionally excluding samples in --exclude-samples. GTF
introns are (exon_i.end+1, exon_{i+1}.start-1), i.e. the same 1-based inclusive intron
boundaries STAR reports, so they compare directly. Matching is strand-agnostic (assembler
strand can flip on ONT cDNA) with +/- --wobble bp tolerance, implemented as a dilated
query against the exact SJ set.

Buckets per tool:
  novel_candidates   class_code not in exclude
  multiexon          >= 1 intron (has a junction to test)
  single_exon        0 introns (no junction; not classifiable by this test)
  sr_supported       multiexon AND every intron is in the SR set
  lr_specific        multiexon AND >= 1 intron is NOT in the SR set   <-- the metric

Brooke Friedman
"""
import argparse
import glob
import os
import re
import sys
from collections import defaultdict

TOOLS = ["bambu", "isoquant", "stringtie"]


def attr(field9, key):
    i = field9.find(key + ' "')
    if i < 0:
        return None
    i += len(key) + 2
    return field9[i:field9.find('"', i)]


def load_sj(folder, min_uniq, exclude_samples):
    """Exact set of (chrom, start, end) short-read introns, pooled uniq >= min_uniq."""
    pooled = defaultdict(int)
    files = sorted(glob.glob(os.path.join(folder, "*.SJ.out.tab")))
    used = []
    for path in files:
        base = os.path.basename(path).replace(".SJ.out.tab", "")
        if base in exclude_samples:
            continue
        used.append(base)
        with open(path) as fh:
            for ln in fh:
                c = ln.split("\t")
                if len(c) < 7:
                    continue
                pooled[(c[0], int(c[1]), int(c[2]))] += int(c[6])
    idx = {k for k, v in pooled.items() if v >= min_uniq}
    print("* SJ: %d samples pooled (%s); %d junctions total, %d with uniq>=%d"
          % (len(used), ",".join(used), len(pooled), len(idx), min_uniq), file=sys.stderr)
    return idx


def sr_hit(chrom, s, e, sj, wobble):
    for ds in range(-wobble, wobble + 1):
        for de in range(-wobble, wobble + 1):
            if (chrom, s + ds, e + de) in sj:
                return True
    return False


def parse_annotated(path, exclude_codes):
    """Yield (tool, [introns]) for each NOVEL transcript in the gffcmp-annotated GTF.
    Exons are accumulated per transcript first (the GTF is not guaranteed tx-grouped)."""
    exclude = set(exclude_codes.split(","))
    tx_exons = {}
    tx_tool = {}
    tx_novel = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 9:
                continue
            if f[2] == "transcript":
                tid = attr(f[8], "transcript_id")
                cc = attr(f[8], "class_code")
                tx_novel[tid] = (cc not in exclude) if cc is not None else True
            elif f[2] == "exon":
                tid = attr(f[8], "transcript_id")
                tx_exons.setdefault(tid, (f[0], [])) [1].append((int(f[3]), int(f[4])))
                if tid not in tx_tool:
                    tx_tool[tid] = tid.split("__", 1)[0]
    for tid, (chrom, exons) in tx_exons.items():
        if not tx_novel.get(tid, False):
            continue
        exons.sort()
        introns = [(chrom, a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:])]
        yield tx_tool.get(tid, "NA"), introns


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotated", required=True, help="<dc>_gffcmp.annotated.gtf")
    ap.add_argument("--sj-folder", required=True, help="folder of STAR *.SJ.out.tab")
    ap.add_argument("--sj-min-uniq", type=int, default=1,
                    help="min pooled uniquely-mapped reads for a SR junction (default 1)")
    ap.add_argument("--exclude-samples", default="",
                    help="comma-separated sample basenames to drop from the SJ pool")
    ap.add_argument("--exclude-codes", default="=,c")
    ap.add_argument("--wobble", type=int, default=2)
    ap.add_argument("--cohort", default="")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    excl = {s for s in args.exclude_samples.split(",") if s}
    sj = load_sj(args.sj_folder, args.sj_min_uniq, excl)

    counts = {t: defaultdict(int) for t in TOOLS + ["NA"]}
    for tool, introns in parse_annotated(args.annotated, args.exclude_codes):
        c = counts[tool]
        c["novel_candidates"] += 1
        if not introns:
            c["single_exon"] += 1
            continue
        c["multiexon"] += 1
        all_in_sr = all(sr_hit(ch, s, e, sj, args.wobble) for ch, s, e in introns)
        if all_in_sr:
            c["sr_supported"] += 1
        else:
            c["lr_specific"] += 1

    fields = ["novel_candidates", "multiexon", "single_exon", "sr_supported", "lr_specific"]
    with open(args.output, "w") as out:
        if args.cohort:
            out.write("# cohort\t%s\tsj_min_uniq=%d\twobble=%d\n"
                      % (args.cohort, args.sj_min_uniq, args.wobble))
        hdr = ["tool"] + fields + ["pct_lr_specific_of_multiexon", "pct_lr_specific_of_novel"]
        out.write("\t".join(hdr) + "\n")
        totals = defaultdict(int)
        for t in TOOLS:
            c = counts[t]
            for k in fields:
                totals[k] += c[k]
            me, lr, nc = c["multiexon"], c["lr_specific"], c["novel_candidates"]
            row = [t] + [str(c[k]) for k in fields]
            row += ["%.1f" % (100 * lr / me) if me else "",
                    "%.1f" % (100 * lr / nc) if nc else ""]
            out.write("\t".join(row) + "\n")
        me, lr, nc = totals["multiexon"], totals["lr_specific"], totals["novel_candidates"]
        row = ["total"] + [str(totals[k]) for k in fields]
        row += ["%.1f" % (100 * lr / me) if me else "",
                "%.1f" % (100 * lr / nc) if nc else ""]
        out.write("\t".join(row) + "\n")

    print("* longread_specific -> %s" % args.output, file=sys.stderr)
    for t in TOOLS:
        c = counts[t]
        me = c["multiexon"] or 1
        print("  %-9s novel=%d  multiexon=%d  lr_specific=%d (%.1f%% of multiexon)"
              % (t, c["novel_candidates"], c["multiexon"], c["lr_specific"],
                 100 * c["lr_specific"] / me), file=sys.stderr)


if __name__ == "__main__":
    main()
