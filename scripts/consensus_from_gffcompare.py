#!/usr/bin/env python
"""
consensus_from_gffcompare.py

Build a cross-tool CONSENSUS GTF from a multi-query gffcompare run.

We assemble each group's pooled BAM with three tools (bambu, isoquant, stringtie3)
and want to keep only isoforms that >= N of them agree on. gffcompare, given the
three GTFs as queries, merges transcripts with the SAME INTRON CHAIN into one
"transfrag" (TCONS_*) and records, per query, whether that query contained a match.
That presence/absence lives in the .tracking file:

    TCONS_id  XLOC_locus  ref_gene|ref_tx  class_code  q1:...  q2:...  q3:...

columns 5+ are the per-query (per-tool) hits, in the ORDER the GTFs were passed on
the gffcompare command line. A column is "-" when that tool has no matching
transcript. So the number of tools supporting a transfrag = count of non-"-" query
columns. We keep TCONS ids with support >= --min-tools and emit their records from
the gffcompare .combined.gtf.

Intron-chain note: gffcompare matches MULTI-exon transcripts by exact intron chain
(ignoring 5'/3' ends, which is what we want for ONT end-fuzziness). SINGLE-exon
transcripts have no intron chain and are matched by overlap instead, a looser
criterion; use --multi-exon-only to drop single-exon transfrags from the consensus
(handle them separately with orthogonal evidence, e.g. CAGE / short-read support).

Brooke Friedman
"""
import argparse


def kept_tcons(tracking_path, min_tools):
    """Return the set of TCONS ids supported by >= min_tools query columns."""
    keep = set()
    with open(tracking_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            tcons = fields[0]
            queries = fields[4:]               # one column per input GTF (tool)
            support = sum(1 for q in queries if q.strip() != "-")
            if support >= min_tools:
                keep.add(tcons)
    return keep


def transcript_id(attr_field):
    """Pull transcript_id "TCONS_..." out of a GTF attribute column."""
    for part in attr_field.split(";"):
        part = part.strip()
        if part.startswith("transcript_id"):
            # transcript_id "TCONS_00000001"
            return part.split('"')[1]
    return None


def multi_exon_tcons(combined_path):
    """TCONS ids that have >1 exon in the combined GTF (i.e. have an intron chain)."""
    exon_counts = {}
    with open(combined_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            tid = transcript_id(cols[8])
            if tid:
                exon_counts[tid] = exon_counts.get(tid, 0) + 1
    return {t for t, n in exon_counts.items() if n > 1}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracking", required=True, help="gffcompare <prefix>.tracking")
    ap.add_argument("--combined", required=True, help="gffcompare <prefix>.combined.gtf")
    ap.add_argument("--min-tools", type=int, default=2,
                    help="min number of tools a transcript must appear in (default 2)")
    ap.add_argument("--multi-exon-only", action="store_true",
                    help="drop single-exon transfrags (matched by overlap, not intron chain)")
    ap.add_argument("-o", "--output", required=True, help="consensus GTF out")
    args = ap.parse_args()

    keep = kept_tcons(args.tracking, args.min_tools)
    if args.multi_exon_only:
        keep &= multi_exon_tcons(args.combined)

    n_written = 0
    seen = set()
    with open(args.combined) as fh, open(args.output, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            tid = transcript_id(cols[8])
            if tid in keep:
                out.write(line)
                if tid not in seen:
                    seen.add(tid)
                    n_written += 1

    print("* %d transfrags supported by >= %d tools%s -> %d transcripts written to %s"
          % (len(keep), args.min_tools,
             " (multi-exon only)" if args.multi_exon_only else "",
             n_written, args.output))


if __name__ == "__main__":
    main()
