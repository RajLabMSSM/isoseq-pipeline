#!/usr/bin/env python
"""
build_reference_with_novel.py

One quantification reference = full reference GTF (preserved verbatim) + only the
truly-novel consensus isoforms on top. The merged-consensus GTF carries a gffcompare
class_code per transcript (vs the reference); transcripts already represented in the
reference (exact '=' and contained 'c' by default) are dropped so the reference copy
is the one kept. The rest are emitted renamed <prefix>_NNNNNN with source column set
to <source>.

Gene assignment is delegated to assign_novel_genes.py: the parent gene is the
same-strand reference gene sharing the most exonic sequence, and a transcript with no
same-strand exonic overlap gets a fresh <prefix>G_NNNNNN locus rather than the nearest
gene. Inheriting the gene from gffcompare's cmp_ref, which this script used to do, put
antisense (x, s) and intron-contained (i, y) transcripts into the isoform pool of a gene
they share no sequence with -- 6.3% of the v8 novels, measured. gffcompare's choice is
still recorded in ref_gene so the two can be compared.
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assign_novel_genes import ExonIndex, assign, opener


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def base(tx):
    return tx.split(".")[0] if tx else tx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged", required=True, help="gffcompare merged consensus GTF (has class_code)")
    ap.add_argument("--reference", required=True, help="full reference GTF, passed through unchanged")
    ap.add_argument("--exclude-codes", default="=,c",
                    help="class codes already represented in the reference -> dropped (default '=,c')")
    ap.add_argument("--prefix", default="NOVEL")
    ap.add_argument("--source", default="ONT", help="GTF source column for novel records")
    ap.add_argument("--membership", help="per-TCONS tool-membership TSV from "
                    "consensus_from_gffcompare.py --membership-out (keyed on the merged "
                    "consensus transcript_id); enables --provenance-out")
    ap.add_argument("--provenance-out", help="TSV: per FINAL novel isoform (renamed id) "
                    "with its tool-support columns, joinable to the output GTF. Requires --membership")
    ap.add_argument("--keep-list", help="file of original merged transcript_ids to keep "
                    "(e.g. SQANTI-passed ids); novels not listed are dropped")
    ap.add_argument("--novel-only", action="store_true",
                    help="emit only the novel records, skip the reference passthrough")
    ap.add_argument("--no-rename", action="store_true",
                    help="keep original transcript_id/gene_id (for SQANTI input) instead "
                    "of renaming to <prefix>_*")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    drop = set(c for c in args.exclude_codes.split(",") if c)
    keep = None
    if args.keep_list:
        keep = set()
        with open(args.keep_list) as fh:
            for line in fh:
                t = line.strip()
                if t:
                    keep.add(t)

    # per-TCONS tool membership (header + 0/1 tool columns) carried onto the final ids
    member_hdr, member = None, {}
    if args.membership:
        with open(args.membership) as fh:
            member_hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                member[f[0]] = f

    # reference transcript_id (versionless) -> gene_id, for the cmp_ref_gene record, and
    # the same-strand exonic index the assignment itself is made against
    tx2gene = {}
    index = ExonIndex()
    ref_exons, ref_meta = {}, {}
    with opener(args.reference) as ref:
        for line in ref:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9:
                continue
            if c[2] == "transcript":
                tx2gene[base(attr(c[8], "transcript_id"))] = attr(c[8], "gene_id")
            elif c[2] == "exon":
                gid = attr(c[8], "gene_id")
                ref_exons.setdefault(gid, []).append((int(c[3]), int(c[4])))
                ref_meta[gid] = (c[0], c[6])
    for gid, exons in ref_exons.items():
        index.add_gene(gid, ref_meta[gid][0], ref_meta[gid][1], exons)
    index.build()
    del ref_exons, ref_meta

    # pass 1: pick the transcripts to keep, in file order so the renamed ids stay stable
    order, cc_of, cmp_of = [], {}, {}
    with opener(args.merged) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "transcript":
                continue
            cc = attr(c[8], "class_code")
            if cc is None or cc in drop:
                continue
            tid = attr(c[8], "transcript_id")
            if keep is not None and tid not in keep:
                continue
            order.append(tid)
            cc_of[tid] = cc
            cmp_of[tid] = tx2gene.get(base(attr(c[8], "cmp_ref"))) or "NA"

    wanted = set(order)
    novel_exons, novel_meta = {}, {}
    with opener(args.merged) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "exon":
                continue
            tid = attr(c[8], "transcript_id")
            if tid not in wanted:
                continue
            novel_exons.setdefault(tid, []).append((int(c[3]), int(c[4])))
            novel_meta.setdefault(tid, (c[0], c[6], attr(c[8], "gene_id")))

    placed = assign([(t, novel_meta[t][0], novel_meta[t][1], sorted(novel_exons[t]))
                     for t in order if t in novel_exons], index, prefix=args.prefix)

    tx_map = {}
    n_tx = 0
    for tid in order:
        n_tx += 1
        gid, how, bridging = placed.get(tid, (None, "unplaced", []))
        if gid is None:
            gid = novel_meta.get(tid, (None, None, None))[2]
        if args.no_rename:
            tx_map[tid] = (tid, novel_meta[tid][2], cmp_of[tid], cc_of[tid], how, bridging)
        else:
            tx_map[tid] = ("%s_%06d" % (args.prefix, n_tx), gid, cmp_of[tid],
                           cc_of[tid], how, bridging)
    n_locus = len(set(v[1] for v in tx_map.values() if v[4] == "new_locus"))

    # per-isoform tool provenance keyed on the FINAL renamed id (joins to the output GTF)
    if args.provenance_out and member_hdr:
        carry = member_hdr[3:]   # n_tools + the 0/1 tool columns
        with open(args.provenance_out, "w") as out:
            out.write("transcript_id\tgene_id\tref_gene\tclass_code\tgene_assignment"
                      "\tbridging\t%s\n" % "\t".join(carry))
            for tid, (new_tid, gid, ref_gene, cc, how, bridging) in tx_map.items():
                row = member.get(tid)
                vals = row[3:] if row else ["NA"] * len(carry)
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                          % (new_tid, gid, ref_gene, cc, how, ",".join(bridging),
                             "\t".join(vals)))

    # pass 2: reference verbatim (unless --novel-only), then the novel records
    with open(args.output, "w") as out:
        if not args.novel_only:
            with opener(args.reference) as ref:
                for line in ref:
                    out.write(line)
        with opener(args.merged) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 9 or c[2] not in ("transcript", "exon"):
                    continue
                tid = attr(c[8], "transcript_id")
                if tid not in tx_map:
                    continue
                new_tid, gid, ref_gene, cc, how, bridging = tx_map[tid]
                c[1] = args.source
                if args.no_rename:
                    c[8] = 'gene_id "%s"; transcript_id "%s";' % (gid, new_tid)
                elif c[2] == "exon":
                    c[8] = 'gene_id "%s"; transcript_id "%s";' % (gid, new_tid)
                else:
                    c[8] = ('gene_id "%s"; transcript_id "%s"; ref_gene "%s"; class_code "%s"; '
                            'gene_assignment "%s"; bridging "%s";'
                            % (gid, new_tid, ref_gene, cc, how, ",".join(bridging)))
                out.write("\t".join(c) + "\n")

    known = sum(1 for v in tx_map.values() if v[4] == "exonic_overlap")
    bridged = sum(1 for v in tx_map.values() if v[5])
    print("* reference + %d novel transcripts: %d are isoforms of known genes, %d new loci "
          "(no same-strand exonic overlap), %d bridging; dropped codes: %s"
          % (n_tx, known, n_locus, bridged, ",".join(sorted(drop))))


if __name__ == "__main__":
    main()
