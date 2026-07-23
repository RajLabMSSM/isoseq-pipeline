#!/usr/bin/env python
"""
build_reference_with_novel.py

One quantification reference = full reference GTF (preserved verbatim) + only the
truly-novel consensus isoforms on top. The merged-consensus GTF carries a gffcompare
class_code per transcript (vs the reference); transcripts already represented in the
reference (exact '=' and contained 'c' by default) are dropped so the reference copy
is the one kept. The rest are emitted renamed <prefix>_NNNNNN with source column set
to <source>.

Gene assignment: a novel transcript that overlaps a known gene (has cmp_ref) inherits
that gene's reference gene_id, so it aggregates as a NOVEL ISOFORM of the known gene.
Only transcripts with no reference overlap (intergenic, class 'u') get a fresh
<prefix>G_NNNNNN locus. (Most novel are isoforms of known genes, not new genes.)
"""
import argparse
import re


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

    # reference transcript_id (versionless) -> gene_id, to map a novel's cmp_ref to its gene
    tx2gene = {}
    with open(args.reference) as ref:
        for line in ref:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "transcript":
                continue
            tx2gene[base(attr(c[8], "transcript_id"))] = attr(c[8], "gene_id")

    # pass 1: assign ids; novel-of-known inherits parent gene, intergenic gets a new locus
    tx_map = {}
    novelg = {}
    n_tx = n_locus = 0
    with open(args.merged) as fh:
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
            n_tx += 1
            parent = tx2gene.get(base(attr(c[8], "cmp_ref")))
            if parent:
                gid, ref_gene = parent, parent
            else:
                xloc = attr(c[8], "gene_id")
                if xloc not in novelg:
                    n_locus += 1
                    novelg[xloc] = "%sG_%06d" % (args.prefix, n_locus)
                gid, ref_gene = novelg[xloc], "NA"
            if args.no_rename:
                tx_map[tid] = (tid, attr(c[8], "gene_id"), ref_gene, cc)
            else:
                tx_map[tid] = ("%s_%06d" % (args.prefix, n_tx), gid, ref_gene, cc)

    # per-isoform tool provenance keyed on the FINAL renamed id (joins to the output GTF)
    if args.provenance_out and member_hdr:
        carry = member_hdr[3:]   # n_tools + the 0/1 tool columns
        with open(args.provenance_out, "w") as out:
            out.write("transcript_id\tgene_id\tref_gene\tclass_code\t%s\n" % "\t".join(carry))
            for tid, (new_tid, gid, ref_gene, cc) in tx_map.items():
                row = member.get(tid)
                vals = row[3:] if row else ["NA"] * len(carry)
                out.write("%s\t%s\t%s\t%s\t%s\n" % (new_tid, gid, ref_gene, cc, "\t".join(vals)))

    # pass 2: reference verbatim (unless --novel-only), then the novel records
    with open(args.output, "w") as out:
        if not args.novel_only:
            with open(args.reference) as ref:
                for line in ref:
                    out.write(line)
        with open(args.merged) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 9 or c[2] not in ("transcript", "exon"):
                    continue
                tid = attr(c[8], "transcript_id")
                if tid not in tx_map:
                    continue
                new_tid, gid, ref_gene, cc = tx_map[tid]
                c[1] = args.source
                if args.no_rename:
                    c[8] = 'gene_id "%s"; transcript_id "%s";' % (gid, new_tid)
                else:
                    c[8] = ('gene_id "%s"; transcript_id "%s"; ref_gene "%s"; class_code "%s";'
                            % (gid, new_tid, ref_gene, cc))
                out.write("\t".join(c) + "\n")

    print("* reference + %d novel transcripts: %d are isoforms of known genes, %d new intergenic loci; dropped codes: %s"
          % (n_tx, n_tx - sum(1 for v in tx_map.values() if v[2] == "NA"), n_locus, ",".join(sorted(drop))))


if __name__ == "__main__":
    main()
