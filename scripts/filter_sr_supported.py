#!/usr/bin/env python
"""
filter_sr_supported.py

Derive a short-read-supported variant of the cross-cohort reference for Salmon
quantification. Known transcripts (ENST) are always kept. A novel transcript
(XNOVEL_) is kept only when it is multi-exon AND every one of its introns is
present in the pooled short-read STAR junction set (SJ.out.tab). Mono-exon
novels have no junction to corroborate and are dropped.

The "every intron in SR" rule is the sr_supported bucket of longread_specific.py;
SJ loading, intron boundary convention (exon_i.end+1, exon_{i+1}.start-1) and the
+/- wobble strand-agnostic match are shared with that script for consistency.

Outputs, next to the input reference:
  <base>.sr_supported.gtf.gz          filtered GTF
  <base>.sr_supported.fa.gz           filtered transcript FASTA
  <base>.sr_supported.provenance.tsv  per-novel decision + reason
  <base>.sr_supported.summary.tsv     counts
"""
import argparse
import glob
import gzip
import os
import re
import sys
from collections import defaultdict

TXID_RE = re.compile(r'transcript_id "([^"]+)"')
GENEID_RE = re.compile(r'gene_id "([^"]+)"')


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_sj(folders, min_uniq):
    pooled = defaultdict(int)
    files = []
    for folder in folders:
        files.extend(sorted(glob.glob(os.path.join(folder, "*.SJ.out.tab"))))
    for path in files:
        with open(path) as fh:
            for ln in fh:
                c = ln.split("\t")
                if len(c) < 7:
                    continue
                pooled[(c[0], int(c[1]), int(c[2]))] += int(c[6])
    idx = {k for k, v in pooled.items() if v >= min_uniq}
    print("* SJ: %d files pooled; %d junctions, %d with uniq>=%d"
          % (len(files), len(pooled), len(idx), min_uniq), file=sys.stderr)
    return idx


def sr_hit(chrom, s, e, sj, wobble):
    for ds in range(-wobble, wobble + 1):
        for de in range(-wobble, wobble + 1):
            if (chrom, s + ds, e + de) in sj:
                return True
    return False


def introns_in_sr(chrom, exons, sj, wobble):
    exons.sort()
    total = 0
    hit = 0
    for a, b in zip(exons, exons[1:]):
        total += 1
        if sr_hit(chrom, a[1] + 1, b[0] - 1, sj, wobble):
            hit += 1
    return total, hit


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, help="cross-cohort reference GTF (.gtf/.gtf.gz)")
    ap.add_argument("--fasta", required=True, help="cross-cohort reference FASTA (.fa/.fa.gz)")
    ap.add_argument("--sj-folder", action="append", required=True,
                    help="folder of STAR *.SJ.out.tab (repeatable)")
    ap.add_argument("--sj-min-uniq", type=int, default=1)
    ap.add_argument("--wobble", type=int, default=2)
    ap.add_argument("--novel-prefix", default="XNOVEL_")
    ap.add_argument("--out-prefix", required=True,
                    help="output path prefix; writes .sr_supported.{gtf.gz,fa.gz,provenance.tsv,summary.tsv}")
    args = ap.parse_args()

    sj = load_sj(args.sj_folder, args.sj_min_uniq)

    tx_chrom = {}
    tx_exons = defaultdict(list)
    print("* pass 1: collecting exons", file=sys.stderr)
    with opener(args.gtf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = TXID_RE.search(f[8])
            if not m:
                continue
            tid = m.group(1)
            if not tid.startswith(args.novel_prefix):
                continue
            tx_chrom[tid] = f[0]
            tx_exons[tid].append((int(f[3]), int(f[4])))

    keep = set()
    prov = []
    n_mono = n_multi_pass = n_multi_fail = 0
    for tid, exons in tx_exons.items():
        if len(exons) < 2:
            n_mono += 1
            prov.append((tid, "novel_monoexon", 0, 0, "drop", "no_junction"))
            continue
        total, hit = introns_in_sr(tx_chrom[tid], exons, sj, args.wobble)
        if hit == total:
            keep.add(tid)
            n_multi_pass += 1
            prov.append((tid, "novel_multiexon", total, hit, "keep", "all_introns_sr"))
        else:
            n_multi_fail += 1
            prov.append((tid, "novel_multiexon", total, hit, "drop", "missing_intron_in_sr"))

    op = args.out_prefix
    prov_path = op + ".sr_supported.provenance.tsv"
    with open(prov_path, "w") as out:
        out.write("transcript_id\tcategory\tn_introns\tn_introns_in_sr\tdecision\treason\n")
        for r in prov:
            out.write("\t".join(str(x) for x in r) + "\n")

    gtf_path = op + ".sr_supported.gtf.gz"
    n_known_tx = set()
    print("* pass 2: writing filtered GTF", file=sys.stderr)
    with opener(args.gtf) as fh, gzip.open(gtf_path, "wt") as out:
        for ln in fh:
            if ln.startswith("#"):
                out.write(ln)
                continue
            f = ln.split("\t")
            if len(f) < 9:
                out.write(ln)
                continue
            m = TXID_RE.search(f[8])
            if not m:
                out.write(ln)  # gene-level feature: keep
                continue
            tid = m.group(1)
            if not tid.startswith(args.novel_prefix):
                out.write(ln)
                n_known_tx.add(tid)
            elif tid in keep:
                out.write(ln)

    fa_path = op + ".sr_supported.fa.gz"
    print("* pass 3: writing filtered FASTA", file=sys.stderr)
    kept_fa = 0
    with opener(args.fasta) as fh, gzip.open(fa_path, "wt") as out:
        emit = False
        for ln in fh:
            if ln.startswith(">"):
                tid = ln[1:].split()[0]
                emit = (not tid.startswith(args.novel_prefix)) or (tid in keep)
                if emit:
                    kept_fa += 1
            if emit:
                out.write(ln)

    summ_path = op + ".sr_supported.summary.tsv"
    n_known = len(n_known_tx)
    with open(summ_path, "w") as out:
        out.write("metric\tvalue\n")
        out.write("known_kept\t%d\n" % n_known)
        out.write("novel_monoexon_dropped\t%d\n" % n_mono)
        out.write("novel_multiexon_kept\t%d\n" % n_multi_pass)
        out.write("novel_multiexon_dropped\t%d\n" % n_multi_fail)
        out.write("novel_kept_total\t%d\n" % n_multi_pass)
        out.write("total_transcripts_kept\t%d\n" % (n_known + n_multi_pass))
        out.write("fasta_records_kept\t%d\n" % kept_fa)
        out.write("sj_min_uniq\t%d\n" % args.sj_min_uniq)
        out.write("wobble\t%d\n" % args.wobble)

    print("* done", file=sys.stderr)
    print("  known kept:            %d" % n_known, file=sys.stderr)
    print("  novel multiexon kept:  %d" % n_multi_pass, file=sys.stderr)
    print("  novel multiexon drop:  %d" % n_multi_fail, file=sys.stderr)
    print("  novel monoexon drop:   %d" % n_mono, file=sys.stderr)
    print("  fasta records kept:    %d" % kept_fa, file=sys.stderr)
    print("  -> %s" % gtf_path, file=sys.stderr)
    print("  -> %s" % fa_path, file=sys.stderr)


if __name__ == "__main__":
    main()
