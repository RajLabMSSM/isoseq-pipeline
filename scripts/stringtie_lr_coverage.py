#!/usr/bin/env python
"""
stringtie_lr_coverage.py

Per-tool count of novel candidates with ZERO long-read coverage across the model body
-- the strictest "not long-read supported" flag. Motivation: StringTie --mix builds
models from BOTH short and long reads, so it can emit a novel isoform supported only by
short reads (zero long-read evidence); bambu / IsoQuant are long-read-only, so their
novels are LR-supported by construction and serve as a control (expect ~0 flagged).

"model body" = the transcript's EXONS (introns excluded). A model is flagged
no_lr_coverage when the summed long-read base coverage over ALL its exons, across ALL
long-read samples, is exactly 0 -- i.e. no long read overlaps any exon of the model.

Source of per-tool novels = the consensus gffcompare-annotated GTF
(<dc>_gffcmp.annotated.gtf), same as longread_specific.py / tracking_numbers: each
transcript carries its tool in the id prefix (<tool>__<group>__<origid>) and a
class_code, so "novel candidate" = class_code not in --exclude-codes.

Coverage = samtools bedcov over the per-sample long-read CRAMs (primary alignments only,
-G 0x900 drops secondary+supplementary), summed across samples. CRAMs need --reference.
"""
import argparse
import glob
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

TOOLS = ["bambu", "isoquant", "stringtie"]


def attr(field9, key):
    i = field9.find(key + ' "')
    if i < 0:
        return None
    i += len(key) + 2
    return field9[i:field9.find('"', i)]


def parse_annotated(path, exclude_codes, target_tools):
    """Return {tid: (tool, chrom, [(start,end)...])} for NOVEL transcripts of the target
    tools. Exons accumulated per transcript (GTF is not guaranteed tx-grouped)."""
    exclude = set(exclude_codes.split(","))
    tx_exons = {}
    tx_tool = {}
    tx_novel = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "transcript":
                tid = attr(f[8], "transcript_id")
                cc = attr(f[8], "class_code")
                tx_novel[tid] = (cc not in exclude) if cc is not None else True
            elif f[2] == "exon":
                tid = attr(f[8], "transcript_id")
                d = tx_exons.setdefault(tid, (f[0], []))
                d[1].append((int(f[3]), int(f[4])))
                if tid not in tx_tool:
                    tx_tool[tid] = tid.split("__", 1)[0]
    out = {}
    for tid, (chrom, exons) in tx_exons.items():
        tool = tx_tool.get(tid, "NA")
        if tool not in target_tools:
            continue
        if not tx_novel.get(tid, False):
            continue
        exons.sort()
        out[tid] = (tool, chrom, exons)
    return out


def bedcov_zero(models, alns, reference, mapq, tmpdir, limit=0):
    """Return set of tids whose summed long-read base coverage over all exons is 0.
    Writes one BED line per exon (name = tid), runs samtools bedcov, sums the per-file
    coverage columns per exon, aggregates per tid. `alns` = BAMs (fast, no --reference)
    or CRAMs (needs `reference`); pooled BAMs are strongly preferred (CRAM adds a
    multi-second reference-decode startup per file)."""
    tids = list(models)
    if limit:
        tids = tids[:limit]
    bed = os.path.join(tmpdir, "exons.bed")
    n_exons = 0
    with open(bed, "w") as b:
        for tid in tids:
            _, chrom, exons = models[tid]
            for s, e in exons:
                b.write("%s\t%d\t%d\t%s\n" % (chrom, s - 1, e, tid))  # BED 0-based start
                n_exons += 1
    # sort so bedcov random-access is coherent
    bed_sorted = os.path.join(tmpdir, "exons.sorted.bed")
    subprocess.check_call("sort -k1,1 -k2,2n %s > %s" % (bed, bed_sorted), shell=True)
    ref_args = ["--reference", reference] if reference else []
    cmd = ["samtools", "bedcov", "-G", "0x900", "-Q", str(mapq)] + ref_args + \
          [bed_sorted] + list(alns)
    print("* bedcov: %d exons x %d aln files (mapq>=%d, primary only)"
          % (n_exons, len(alns), mapq), file=sys.stderr)
    cov = defaultdict(int)  # tid -> summed base coverage over exons, all crams
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for ln in proc.stdout:
        c = ln.rstrip("\n").split("\t")
        tid = c[3]
        cov[tid] += sum(int(x) for x in c[4:])  # one column per cram
    if proc.wait() != 0:
        raise RuntimeError("samtools bedcov failed")
    return {tid for tid in tids if cov.get(tid, 0) == 0}, tids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotated", required=True, help="<dc>_gffcmp.annotated.gtf")
    ap.add_argument("--bam", nargs="+", default=[],
                    help="pooled long-read BAM(s) (PREFERRED: fast, no --reference)")
    ap.add_argument("--crams", nargs="+", default=[],
                    help="per-sample long-read CRAMs (needs --reference; slow startup)")
    ap.add_argument("--reference", default="", help="genome fasta (required for --crams)")
    ap.add_argument("--exclude-codes", default="=,c")
    ap.add_argument("--tools", default="stringtie",
                    help="comma-list of tools to score (default stringtie; 'all' = %s)"
                         % ",".join(TOOLS))
    ap.add_argument("--mapq", type=int, default=0, help="min MAPQ for bedcov (default 0)")
    ap.add_argument("--cohort", default="")
    ap.add_argument("--limit", type=int, default=0, help="score only first N models (test)")
    ap.add_argument("--list-out", default="", help="optional: write flagged tids here")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    target = TOOLS if args.tools == "all" else [t for t in args.tools.split(",") if t]
    alns = [a for a in (args.bam or args.crams) if os.path.exists(a)]
    reference = "" if args.bam else args.reference
    if not alns:
        sys.exit("no alignment files found (--bam or --crams)")
    if args.crams and not args.bam and not reference:
        sys.exit("--crams requires --reference")

    models = parse_annotated(args.annotated, args.exclude_codes, set(target))
    per_tool_total = defaultdict(int)
    for tid, (tool, _, _) in models.items():
        per_tool_total[tool] += 1
    print("* novel candidates by tool: %s"
          % {t: per_tool_total[t] for t in target}, file=sys.stderr)

    with tempfile.TemporaryDirectory(dir=os.path.dirname(args.output) or ".") as td:
        zero, scored = bedcov_zero(models, alns, reference, args.mapq, td, args.limit)

    scored_by_tool = defaultdict(int)
    zero_by_tool = defaultdict(int)
    for tid in scored:
        tool = models[tid][0]
        scored_by_tool[tool] += 1
        if tid in zero:
            zero_by_tool[tool] += 1

    with open(args.output, "w") as out:
        if args.cohort:
            out.write("# cohort\t%s\taln_files=%d\tsource=%s\tmapq=%d\tlimit=%d\n"
                      % (args.cohort, len(alns), "bam" if args.bam else "cram",
                         args.mapq, args.limit))
        out.write("tool\tnovel_candidates_scored\tno_lr_coverage\tpct_no_lr_coverage\n")
        tot_s = tot_z = 0
        for t in target:
            s, z = scored_by_tool[t], zero_by_tool[t]
            tot_s += s
            tot_z += z
            out.write("%s\t%d\t%d\t%s\n" % (t, s, z, "%.2f" % (100 * z / s) if s else ""))
        out.write("total\t%d\t%d\t%s\n"
                  % (tot_s, tot_z, "%.2f" % (100 * tot_z / tot_s) if tot_s else ""))

    if args.list_out:
        with open(args.list_out, "w") as lf:
            for tid in sorted(zero):
                lf.write("%s\t%s\n" % (models[tid][0], tid))

    print("* stringtie_lr_coverage -> %s" % args.output, file=sys.stderr)
    for t in target:
        s, z = scored_by_tool[t], zero_by_tool[t]
        print("  %-9s scored=%d  no_lr_coverage=%d (%.2f%%)"
              % (t, s, z, 100 * z / s if s else 0), file=sys.stderr)


if __name__ == "__main__":
    main()
