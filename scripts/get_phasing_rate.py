#!/usr/bin/env python3
"""
get_phasing_rate.py
Usage: python get_phasing_rate.py --results-dir results --out phasing_rates.tsv [-t 8]
"""

import argparse
import os
import sys
import pysam
import concurrent.futures


def get_phasing_rate(args):
    sample, bam_path = args
    try:
        total = 0
        usable = 0
        phased = 0

        with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                total += 1
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < 10:
                    continue
                usable += 1
                if read.has_tag("HP"):
                    phased += 1

        usable_rate = usable / total if total > 0 else 0.0
        phasing_rate = phased / usable if usable > 0 else 0.0
        return (sample, total, usable, phased, usable_rate, phasing_rate, None)

    except Exception as e:
        return (sample, None, None, None, None, None, str(e))


def discover_bams(results_dir, sample_list=None):
    candidates = sample_list if sample_list else sorted(os.listdir(results_dir))
    found = []
    for sample in candidates:
        phased_dir = os.path.join(results_dir, sample, "phased")
        for suffix in (f"{sample}.phased.bam", f"{sample}.bam"):
            bam_path = os.path.join(phased_dir, suffix)
            if os.path.exists(bam_path):
                found.append((sample, bam_path))
                break
        else:
            print(f"  Skipping {sample}: no BAM found in {phased_dir}",
                  file=sys.stderr, flush=True)
    return found


HEADER = ["sample", "total_reads", "usable_reads", "phased_reads",
          "usable_rate", "phasing_rate"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--samples", nargs="+", default=None)
    args = parser.parse_args()

    print(f"Discovering BAMs in {args.results_dir}...", flush=True)
    bams = discover_bams(args.results_dir, args.samples)
    print(f"  Found {len(bams)} BAM(s)", flush=True)

    if not bams:
        print("Error: no BAMs found.", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(bams)} BAM(s) with {args.threads} worker(s)...", flush=True)
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        for i, result in enumerate(executor.map(get_phasing_rate, bams), 1):
            sample, total, usable, phased, usable_rate, phasing_rate, err = result
            if err:
                print(f"  [{i}/{len(bams)}] {sample}: ERROR — {err}",
                      file=sys.stderr, flush=True)
            else:
                print(f"  [{i}/{len(bams)}] {sample}: "
                      f"total={total:,} usable={usable:,} ({usable_rate*100:.1f}%) "
                      f"phased={phased:,} ({phasing_rate*100:.1f}% of usable)", flush=True)
            results.append(result)

    results.sort(key=lambda r: r[0])

    with open(args.out, "w") as f:
        f.write("\t".join(HEADER) + "\n")
        for sample, total, usable, phased, usable_rate, phasing_rate, err in results:
            if err:
                f.write(f"{sample}\t" + "\t".join(["NA"] * 5) + "\n")
            else:
                f.write(f"{sample}\t{total}\t{usable}\t{phased}\t"
                        f"{usable_rate:.4f}\t{phasing_rate:.4f}\n")

    n_ok  = sum(1 for r in results if r[6] is None)
    n_err = sum(1 for r in results if r[6] is not None)
    print(f"Done. {n_ok} succeeded, {n_err} failed.", flush=True)
