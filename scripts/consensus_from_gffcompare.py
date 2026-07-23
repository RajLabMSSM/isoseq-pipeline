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

import wobble


def _support(queries, tools_per_group):
    """Tool count for one query vector. tools_per_group>0 -> max present within
    any per-group chunk (per-group consensus, unioned across groups)."""
    if tools_per_group > 0:
        chunks = [queries[i:i + tools_per_group]
                  for i in range(0, len(queries), tools_per_group)]
        return max(sum(1 for q in c if q.strip() != "-") for c in chunks)
    return sum(1 for q in queries if q.strip() != "-")


def load_tracking(tracking_path):
    """TCONS id -> list of per-query (per-tool x group) columns."""
    t = {}
    with open(tracking_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            t[fields[0].split("|")[0]] = fields[4:]   # one column per input GTF
    return t


def kept_tcons(tracking, min_tools, tools_per_group=0):
    """Set of TCONS ids supported by >= min_tools tools (exact intron chain).

    tools_per_group>0 splits query columns into per-group chunks (GTFs must be
    passed to gffcompare grouped: g1 tool1,tool2,tool3, g2 tool1,...) and keeps a
    transfrag with >= min_tools in ANY single group -- per-group consensus unioned
    across groups (the old two-gffcompare chain collapsed into one pass).
    """
    return {tc: q for tc, q in tracking.items()
            if _support(q, tools_per_group) >= min_tools}


def kept_tcons_wobble(tracking, combined_path, min_tools, tools_per_group, wobble_bp):
    """Like kept_tcons but tolerant of +/- wobble_bp at splice sites: TCONS whose
    intron chains agree within wobble_bp are clustered, their per-tool support is
    POOLED across the cluster (a tool counts if it hit any member), and a cluster
    passing >= min_tools contributes one representative (most-supported, then
    longest chain). Lets near-matches reach 2/3 instead of being split & dropped,
    and collapses near-identical isoforms within the cohort.

    Returns {representative_tcons: pooled_query_vector}.
    """
    chains = wobble.intron_chains(combined_path)
    keep = {}
    for members in wobble.cluster_by_wobble(chains, wobble_bp):
        qvs = [tracking[m] for m in members if m in tracking]
        if not qvs:
            continue
        pooled = ["x" if any(qv[i].strip() != "-" for qv in qvs) else "-"
                  for i in range(len(qvs[0]))]
        if _support(pooled, tools_per_group) >= min_tools:
            rep = max(members, key=lambda m: (
                _support(tracking.get(m, []), tools_per_group),
                len(chains[m][2]), m))
            keep[rep] = pooled
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


def class_codes(combined_path):
    """TCONS id -> gffcompare class_code (from transcript lines in combined.gtf)."""
    cc = {}
    with open(combined_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            tid = transcript_id(cols[8])
            m = wobble.parse_attr(cols[8], "class_code")
            if tid:
                cc[tid] = m
    return cc


def tool_combo(vec, tool_names, tools_per_group):
    """Which tools supported a transfrag, pooled across groups (OR of each tool's
    per-group columns). Returns the sorted subset of tool_names present."""
    n_groups = max(1, len(vec) // tools_per_group) if tools_per_group else 1
    present = []
    for j, name in enumerate(tool_names):
        cols = [j + g * tools_per_group for g in range(n_groups)] if tools_per_group else [j]
        if any(c < len(vec) and vec[c].strip() != "-" for c in cols):
            present.append(name)
    return present


def write_tool_overlap(keep_vectors, classcodes, drop, tool_names, tools_per_group, out_tsv):
    """Tally NOVEL kept isoforms by which combination of tools supported them
    (the data behind a 3-set Venn). Writes a TSV of all subset counts."""
    from itertools import combinations
    counts = {}
    total = 0
    for tcons, vec in keep_vectors.items():
        cc = classcodes.get(tcons)
        if cc is None or cc in drop:          # novel only
            continue
        total += 1
        combo = tuple(tool_combo(vec, tool_names, tools_per_group))
        counts[combo] = counts.get(combo, 0) + 1
    # every non-empty subset, biggest first, for complete Venn data
    subsets = []
    for r in range(len(tool_names), 0, -1):
        subsets.extend(combinations(tool_names, r))
    with open(out_tsv, "w") as out:
        out.write("# NOVEL consensus isoforms by supporting-tool combination "
                  "(presence pooled across groups + wobble cluster)\n")
        out.write("combination\tn_isoforms\n")
        for s in subsets:
            out.write("%s\t%d\n" % ("+".join(s), counts.get(tuple(s), 0)))
        out.write("total_novel\t%d\n" % total)


def write_tool_membership(keep_vectors, classcodes, drop, tool_names, tools_per_group, out_tsv):
    """Per-isoform tool provenance: one row per kept consensus transfrag with a 0/1
    column per assembler (support pooled across groups, same rule as the Venn tally),
    plus n_tools and an is_novel flag (class_code not in `drop`). Drop-in UpSet/Venn
    membership matrix for R; novel rows reconcile with _tool_overlap.tsv."""
    with open(out_tsv, "w") as out:
        out.write("transcript_id\tclass_code\tis_novel\tn_tools\t%s\n"
                  % "\t".join(tool_names))
        for tcons, vec in sorted(keep_vectors.items()):
            cc = classcodes.get(tcons, "")
            present = tool_combo(vec, tool_names, tools_per_group)
            flags = "\t".join("1" if t in present else "0" for t in tool_names)
            is_novel = "0" if (cc and cc in drop) else "1"
            out.write("%s\t%s\t%s\t%d\t%s\n" % (tcons, cc, is_novel, len(present), flags))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracking", required=True, help="gffcompare <prefix>.tracking")
    ap.add_argument("--combined", required=True, help="gffcompare <prefix>.combined.gtf")
    ap.add_argument("--min-tools", type=int, default=2,
                    help="min number of tools a transcript must appear in (default 2)")
    ap.add_argument("--tools-per-group", type=int, default=0,
                    help="if >0, query columns split into per-group chunks of this size; "
                         "keep transfrags with >=min-tools in ANY group (per-group consensus, "
                         "unioned). GTFs must be passed to gffcompare grouped (default 0=off)")
    ap.add_argument("--junction-wobble", type=int, default=0,
                    help="if >0, cluster TCONS whose intron chains agree within this "
                         "many bp and POOL their tool support before the min-tools test "
                         "(near-matches count toward consensus instead of being split & "
                         "dropped; collapses near-identical isoforms). default 0=exact")
    ap.add_argument("--multi-exon-only", action="store_true",
                    help="drop single-exon transfrags (matched by overlap, not intron chain)")
    ap.add_argument("--overlap-out", help="TSV: NOVEL kept isoforms by tool combination (Venn data)")
    ap.add_argument("--membership-out", help="TSV: one row per kept consensus isoform with a "
                    "0/1 column per tool (per-isoform provenance; UpSet/Venn matrix for R)")
    ap.add_argument("--tool-names", default="bambu,isoquant,stringtie",
                    help="tool labels in query order, for --overlap-out (default bambu,isoquant,stringtie)")
    ap.add_argument("--novel-exclude-codes", default="=,c",
                    help="class codes counted as NON-novel for --overlap-out (default '=,c')")
    ap.add_argument("-o", "--output", required=True, help="consensus GTF out")
    args = ap.parse_args()

    tracking = load_tracking(args.tracking)
    if args.junction_wobble > 0:
        keep = kept_tcons_wobble(tracking, args.combined, args.min_tools,
                                 args.tools_per_group, args.junction_wobble)
    else:
        keep = kept_tcons(tracking, args.min_tools, args.tools_per_group)
    if args.multi_exon_only:
        multi = multi_exon_tcons(args.combined)
        keep = {k: v for k, v in keep.items() if k in multi}

    if args.overlap_out or args.membership_out:
        cc = class_codes(args.combined)
        drop = set(c for c in args.novel_exclude_codes.split(",") if c)
        names = args.tool_names.split(",")
        if args.overlap_out:
            write_tool_overlap(keep, cc, drop, names, args.tools_per_group, args.overlap_out)
        if args.membership_out:
            write_tool_membership(keep, cc, drop, names, args.tools_per_group, args.membership_out)

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

    print("* %d transfrags supported by >= %d tools%s%s -> %d transcripts written to %s"
          % (len(keep), args.min_tools,
             " (+-%dbp wobble)" % args.junction_wobble if args.junction_wobble else "",
             " (multi-exon only)" if args.multi_exon_only else "",
             n_written, args.output))


if __name__ == "__main__":
    main()
