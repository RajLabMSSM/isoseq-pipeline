#!/usr/bin/env python3
"""
junction_qc.py  —  splice junction quality summary from a sorted BAM

Extracts intron spans (CIGAR N operations) from mapped reads, checks canonical
splice sites (GT-AG / CT-AC on reverse strand) against a reference FASTA, and
reports per-sample stats that indicate alignment preset quality.

Usage:
    python junction_qc.py --bam sample.bam --fasta genome.fa [--max-reads N]
    python junction_qc.py --bam_dir results/v3/ --fasta genome.fa --out junc_qc.tsv

Output columns:
    sample, mapped_reads, spliced_reads, pct_spliced,
    unique_junctions, junctions_5x+, canonical_junctions,
    pct_canonical, median_intron_size
"""

import argparse
import sys
import os
import re
import statistics
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.exit("pysam not found — activate the isoseq-pipeline conda env first")


def parse_junctions(read):
    """Return list of (chrom, start, end) intron spans from a read's CIGAR."""
    if read.cigartuples is None:
        return []
    junctions = []
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):   # M / = / X  — consumes reference
            pos += length
        elif op == 2:          # D — deletion in reference
            pos += length
        elif op == 3:          # N — intron skip
            junctions.append((read.reference_name, pos, pos + length))
            pos += length
        # I / S / H don't advance reference position
    return junctions


def analyze(bam_path, fasta, max_reads):
    junc_counts  = defaultdict(int)
    intron_sizes = []
    total_mapped = 0
    spliced      = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if total_mapped >= max_reads:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            total_mapped += 1
            juncs = parse_junctions(read)
            if juncs:
                spliced += 1
                for j in juncs:
                    junc_counts[j] += 1
                    intron_sizes.append(j[2] - j[1])

    # canonical splice site check (GT-AG rule; also CT-AC on minus strand reads)
    canonical = non_canonical = 0
    if fasta:
        for (chrom, start, end), count in junc_counts.items():
            try:
                donor    = fasta.fetch(chrom, start,     start + 2).upper()
                acceptor = fasta.fetch(chrom, end   - 2, end      ).upper()
                if (donor == "GT" and acceptor == "AG") or \
                   (donor == "CT" and acceptor == "AC"):
                    canonical += count
                else:
                    non_canonical += count
            except Exception:
                pass

    unique_j = len(junc_counts)
    hc_j     = sum(1 for c in junc_counts.values() if c >= 5)
    med_intron = int(statistics.median(intron_sizes)) if intron_sizes else 0
    total_checked = canonical + non_canonical
    pct_can = 100 * canonical / total_checked if total_checked else float("nan")

    return {
        "mapped_reads":       total_mapped,
        "spliced_reads":      spliced,
        "pct_spliced":        round(100 * spliced / total_mapped, 2) if total_mapped else 0,
        "unique_junctions":   unique_j,
        "junctions_5x+":      hc_j,
        "canonical_junctions":canonical,
        "pct_canonical":      round(pct_can, 2),
        "median_intron_size": med_intron,
    }


def find_bams(bam_dir):
    """Walk bam_dir for *aligned.bam files, return {sample: path}."""
    bams = {}
    for root, _, files in os.walk(bam_dir):
        for f in files:
            if f.endswith(".aligned.bam") and not f.endswith(".bai"):
                sample = f.replace(".aligned.bam", "")
                bams[sample] = os.path.join(root, f)
    return bams


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--bam",     help="single BAM file")
    g.add_argument("--bam_dir", help="directory tree to scan for *.aligned.bam")
    p.add_argument("--fasta",     required=True, help="reference genome FASTA (indexed)")
    p.add_argument("--max_reads", type=int, default=200_000,
                   help="reads per sample to examine (default 200 000)")
    p.add_argument("--out", default="-", help="output TSV path (default: stdout)")
    args = p.parse_args()

    fasta = pysam.FastaFile(args.fasta) if args.fasta else None

    if args.bam:
        bams = {os.path.basename(args.bam).replace(".aligned.bam", ""): args.bam}
    else:
        bams = find_bams(args.bam_dir)
        if not bams:
            sys.exit(f"No *.aligned.bam files found under {args.bam_dir}")

    cols = ["sample", "mapped_reads", "spliced_reads", "pct_spliced",
            "unique_junctions", "junctions_5x+", "canonical_junctions",
            "pct_canonical", "median_intron_size"]

    out = open(args.out, "w") if args.out != "-" else sys.stdout
    print("\t".join(cols), file=out)

    for sample in sorted(bams):
        print(f"  analyzing {sample} ...", file=sys.stderr)
        stats = analyze(bams[sample], fasta, args.max_reads)
        row = [sample] + [str(stats[c]) for c in cols[1:]]
        print("\t".join(row), file=out)

    if fasta:
        fasta.close()
    if args.out != "-":
        out.close()
        print(f"Written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
