#!/usr/bin/env python
"""
phasing_rate.py

Calculate the phasing rate (fraction of reads with an HP tag) for each
longcallR phased BAM file.

Expects the same file structure as the other longcallR scripts:
    results/<sample>/phased/<sample>.phased.bam

Usage:
    python phasing_rate.py \
        --results-dir results \
        --out phasing_rates.tsv \
        [-t 8] \
        [--samples sample1 sample2 ...]
"""

import argparse
import os
import sys
import subprocess
import concurrent.futures


def get_phasing_rate(args):
    sample, bam_path = args
    try:
        # Total reads
        total = int(subprocess.check_output(
            ["samtools", "view", "-c", bam_path],
            stderr=subprocess.DEVNULL
        ).strip())

        # Reads with HP tag — use samtools view -d to filter by tag
        # -d HP selects reads that have the HP tag (any value)
        phased = int(subprocess.check_output(
            ["samtools", "view", "-c", "-d", "HP", bam_path],
            stderr=subprocess.DEVNULL
        ).strip())

        rate = phased / total if total > 0 else 0.0
        return (sample, total, phased, rate, None)

    except subprocess.CalledProcessError as e:
        return (sample, None, None, None, str(e))
    except Exception as e:
        return (sample, None, None, None, str(e))


def discover_bams(results_dir, sample_list=None):
    candidates = sample_list if sample_list else sorted(os.listdir(results_dir))
    found = []
    for sample in candidates:
        phased_dir = os.path.join(results_dir, sample, "phased")
        # Accept either <sample>.phased.bam or <sample>.bam
        for suffix in (f"{sample}.phased.bam", f"{sample}.bam"):
            bam_path = os.path.join(phased_dir, suffix)
            if os.path.exists(bam_path):
                found.append((sample, bam_path))
                break
        else:
            print(f"  Skipping {sample}: no BAM found in {phased_dir}",
                  file=sys.stderr, flush=True)
    return found


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate phasing rate for longcallR phased BAM files."
    )
    parser.add_argument("--results-dir", required=True,
                        help="Root results directory containing per-sample subdirectories")
    parser.add_argument("--out", required=True,
                        help="Output TSV file path")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Parallel workers. Default: 4")
    parser.add_argument("--samples", nargs="+", default=None,
                        help="Optional: restrict to specific sample names")
    args = parser.parse_args()

    print(f"Discovering BAMs in {args.results_dir}...", flush=True)
    bams = discover_bams(args.results_dir, args.samples)
    print(f"  Found {len(bams)} BAM(s)", flush=True)

    if not bams:
        print("Error: no BAMs found.", file=sys.stderr)
        sys.exit(1)

    print(f"Computing phasing rates with {args.threads} worker(s)...", flush=True)
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        for i, result in enumerate(executor.map(get_phasing_rate, bams), 1):
            sample, total, phased, rate, err = result
            if err:
                print(f"  [{i}/{len(bams)}] {sample}: ERROR — {err}",
                      file=sys.stderr, flush=True)
            else:
                print(f"  [{i}/{len(bams)}] {sample}: "
                      f"{phased:,}/{total:,} = {rate*100:.1f}%", flush=True)
            results.append(result)

    # Sort by sample name
    results.sort(key=lambda r: r[0])

    print(f"Writing results to {args.out}...", flush=True)
    with open(args.out, "w") as f:
        f.write("sample\ttotal_reads\tphased_reads\tphasing_rate\n")
        for sample, total, phased, rate, err in results:
            if err:
                f.write(f"{sample}\tNA\tNA\tNA\n")
            else:
                f.write(f"{sample}\t{total}\t{phased}\t{rate:.4f}\n")

    n_ok  = sum(1 for r in results if r[4] is None)
    n_err = sum(1 for r in results if r[4] is not None)
    print(f"Done. {n_ok} succeeded, {n_err} failed.", flush=True)
