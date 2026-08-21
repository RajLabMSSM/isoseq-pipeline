#!/usr/bin/env python
"""
tracking_numbers.py

One per-cohort funnel table documenting how many isoforms survive each step:

  raw per tool (all transcripts, incl. reference passthrough for bambu/isoquant)
    -> novel candidates per tool (class_code not in =,c)
    -> end-aware consensus collapse (dedup: introns+-wobble, 5'+-tss_tol, 3'+-tes_tol)
    -> tool-overlap breakdown of the collapsed novel set (the Venn)
    -> SQANTI rules filter
    -> final novel in the reference.

Emitted as a tidy metric<TAB>value TSV (easy to read and to plot).

Notes on definitions (so the numbers are not misread):
  * "collapsed_removed" = novel_candidates_total - consensus_novel_total. In UNION
    mode (min_consensus_tools=1) NOTHING is dropped by voting; this number is purely
    duplicate/near-duplicate STRUCTURE merged across tools + within-tool end dups.
  * raw_*_all includes the reference transcripts bambu/isoquant pass through verbatim,
    so it is much larger than the novel candidate count -- both are reported.
"""
import argparse
import glob
import gzip
import re
from collections import defaultdict


def opentext(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

TOOLS = ["bambu", "isoquant", "stringtie"]


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def read_stats(path):
    """per-tool raw_all + novel_candidates from consensus_endaware --stats-out."""
    raw, novelcand = {}, {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                raw[f[0]] = int(f[1])
                novelcand[f[0]] = int(f[2])
    return raw, novelcand


def read_membership(path, keep=None):
    """returns (n_total, venn combo->count, per_tool present count).
    keep: optional set of transcript_ids to restrict to (e.g. SQANTI-passed)."""
    venn = defaultdict(int)
    per_tool = defaultdict(int)
    n = 0
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idi = hdr.index("transcript_id") if "transcript_id" in hdr else 0
        ti = {t: hdr.index(t) for t in TOOLS if t in hdr}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            if keep is not None and f[idi] not in keep:
                continue
            n += 1
            present = tuple(t for t in TOOLS if t in ti and f[ti[t]] == "1")
            venn[present] += 1
            for t in present:
                per_tool[t] += 1
    return n, venn, per_tool


def read_pass_ids(path):
    """isoform IDs (the 'isoform' column) from a SQANTI filter classification TSV."""
    ids = set()
    with opentext(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ii = hdr.index("isoform") if "isoform" in hdr else 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > ii:
                ids.add(f[ii])
    return ids


def read_exon_split(path):
    """mono/multi-exon split of the SQANTI-passed novels, plus the purely-intronic
    mono-exon subset. Mono-exons bypass SQANTI's ML classifier, so this split is the
    headline sensitivity axis: every novel count should be quotable with and without them."""
    mono = multi = mono_intronic = 0
    with opentext(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            ei = hdr.index("exons")
        except ValueError:
            return None
        ci = hdr.index("structural_category") if "structural_category" in hdr else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= ei:
                continue
            try:
                n = int(f[ei])
            except ValueError:
                continue
            if n == 1:
                mono += 1
                if ci is not None and len(f) > ci and f[ci].strip().replace("_", " ") == "genic intron":
                    mono_intronic += 1
            else:
                multi += 1
    return mono, multi, mono_intronic


def count_lines_minus_header(path):
    with open(path) as fh:
        return max(0, sum(1 for _ in fh) - 1)


def count_ref_tx(path):
    n = 0
    with opentext(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) > 2 and c[2] == "transcript":
                n += 1
    return n


def read_longread_specific(path):
    """per-tool + total lr_specific / sr_supported counts from a longread_specific.py TSV."""
    lr, sr = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "tool":
                hdr = {name: i for i, name in enumerate(f)}
                continue
            lr[f[0]] = int(f[hdr["lr_specific"]])
            sr[f[0]] = int(f[hdr["sr_supported"]])
    return lr, sr


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats", required=True, help="consensus_endaware --stats-out TSV")
    ap.add_argument("--membership", required=True, help="consensus membership TSV")
    ap.add_argument("--sqanti-classification", required=True,
                    help="<dc>_filter_sqanti_classification.tsv (SQANTI-passed novels)")
    ap.add_argument("--reference", required=True, help="plain reference GTF (for final total)")
    ap.add_argument("--longread-specific", default=None,
                    help="optional longread_specific.py TSV -> adds lr_specific_<tool> rows")
    ap.add_argument("--exclude-codes", default="=,c")
    ap.add_argument("--cohort", default="")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    raw, novelcand = read_stats(args.stats)
    n_consensus, venn, cons_per_tool = read_membership(args.membership)
    sqanti_ids = read_pass_ids(args.sqanti_classification)
    sqanti_pass = len(sqanti_ids)
    _, venn_sq, sq_per_tool = read_membership(args.membership, keep=sqanti_ids)
    n_ref = count_ref_tx(args.reference)

    raw_total = sum(raw.values())
    novelcand_total = sum(novelcand.values())
    from itertools import combinations
    combos = []
    for k in range(len(TOOLS), 0, -1):
        combos.extend(combinations(TOOLS, k))

    rows = []
    if args.cohort:
        rows.append(("cohort", args.cohort))
    for t in TOOLS:
        rows.append(("raw_%s_all" % t, raw.get(t, 0)))
    rows.append(("raw_total_all", raw_total))
    for t in TOOLS:
        rows.append(("novel_candidates_%s" % t, novelcand.get(t, 0)))
    rows.append(("novel_candidates_total", novelcand_total))
    rows.append(("consensus_novel_total", n_consensus))
    rows.append(("collapsed_removed", novelcand_total - n_consensus))
    for c in combos:
        rows.append(("overlap_%s" % "+".join(c), venn.get(tuple(c), 0)))
    for t in TOOLS:
        rows.append(("consensus_supported_by_%s" % t, cons_per_tool.get(t, 0)))
    for c in combos:
        rows.append(("sqanti_overlap_%s" % "+".join(c), venn_sq.get(tuple(c), 0)))
    for t in TOOLS:
        rows.append(("sqanti_pass_%s" % t, sq_per_tool.get(t, 0)))
    for k in (1, 2, 3):
        rows.append(("sqanti_pass_%dtool" % k, sum(v for cc, v in venn_sq.items() if len(cc) == k)))
    if args.longread_specific:
        lr, sr = read_longread_specific(args.longread_specific)
        for t in TOOLS:
            rows.append(("lr_specific_%s" % t, lr.get(t, 0)))
        rows.append(("lr_specific_total", sum(lr.get(t, 0) for t in TOOLS)))
        for t in TOOLS:
            rows.append(("sr_supported_%s" % t, sr.get(t, 0)))
        rows.append(("sr_supported_total", sum(sr.get(t, 0) for t in TOOLS)))
    rows.append(("novel_pre_sqanti", n_consensus))
    rows.append(("sqanti_pass", sqanti_pass))
    rows.append(("sqanti_removed", n_consensus - sqanti_pass))
    rows.append(("final_novel", sqanti_pass))
    split = read_exon_split(args.sqanti_classification)
    if split is not None:
        mono, multi, mono_intronic = split
        rows.append(("final_novel_multiexon", multi))
        rows.append(("final_novel_monoexon", mono))
        rows.append(("final_novel_monoexon_genic_intron", mono_intronic))
        rows.append(("pct_novel_monoexon", round(100.0 * mono / sqanti_pass, 1) if sqanti_pass else 0))
        rows.append(("pct_monoexon_genic_intron", round(100.0 * mono_intronic / mono, 1) if mono else 0))
    rows.append(("reference_transcripts", n_ref))
    rows.append(("final_total", n_ref + sqanti_pass))
    if split is not None:
        rows.append(("final_novel_EXCL_monoexon", multi))
        rows.append(("final_total_EXCL_monoexon", n_ref + multi))

    with open(args.output, "w") as out:
        out.write("metric\tvalue\n")
        for k, v in rows:
            out.write("%s\t%s\n" % (k, v))
    print("* tracking_numbers -> %s (%d metrics)" % (args.output, len(rows)))


if __name__ == "__main__":
    main()
