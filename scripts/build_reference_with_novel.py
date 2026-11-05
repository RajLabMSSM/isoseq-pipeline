#!/usr/bin/env python
"""The reference gtf verbatim plus the novel merged-consensus models (class code not in
--exclude-codes, and in --keep-list), renamed <prefix>_NNNNNN with source <source>.

The gene is SQANTI3's associated_gene; a SQANTI novelGene becomes a fresh <prefix>G_NNNNNN locus.
When several genes are named, the gene of associated_transcript leads (a fusion: the first gene),
the others go to bridging. --novel-only --no-rename writes the SQANTI input instead."""
import argparse
import gzip
import re
from collections import namedtuple

# one novel model as written out: new id, gene, reference gene of its gffcompare match, class code,
# how the gene was assigned (SQANTI structural category or "unclassified"), further genes it bridges
NovelRecord = namedtuple("NovelRecord",
                         ["new_tid", "gene_id", "ref_gene", "class_code", "gene_assignment", "bridging"])
# membership table: transcript_id, class_code, is_novel, then the per-tool columns carried into provenance
MEMBERSHIP_FIRST_CARRIED_COLUMN = 3


def opener(path):
    # open plain or gzipped text the same way
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def attr(field, key):
    # value of one gtf attribute (key "value";), or None when the key is missing
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def base(tx):
    # transcript id without its ".version" suffix
    return tx.split(".")[0] if tx else tx


def gtf_records(path, features):
    # yield the split columns of every gtf line whose feature type is in features
    with opener(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if not line.startswith("#") and len(c) >= 9 and c[2] in features:
                yield c


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--sqanti-classification")
    ap.add_argument("--exclude-codes", required=True)
    ap.add_argument("--prefix", default="NOVEL")
    ap.add_argument("--source", default="ONT")
    ap.add_argument("--membership")
    ap.add_argument("--provenance-out")
    ap.add_argument("--keep-list")
    ap.add_argument("--novel-only", action="store_true")
    ap.add_argument("--no-rename", action="store_true")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    if not args.no_rename and not args.sqanti_classification:
        raise SystemExit("--sqanti-classification is required unless --no-rename")

    # class codes to leave out, and the optional list of model ids to keep (SQANTI-passed models)
    drop = set(c for c in args.exclude_codes.split(",") if c)
    keep = None
    if args.keep_list:
        with open(args.keep_list) as fh:
            keep = {line.strip() for line in fh if line.strip()}

    # reference transcript (no version) -> gene id
    tx2gene = {base(attr(c[8], "transcript_id")): attr(c[8], "gene_id")
               for c in gtf_records(args.reference, ("transcript",))}

    # SQANTI classification: isoform -> (structural category, associated gene, associated transcript)
    sq = {}
    if args.sqanti_classification:
        with opener(args.sqanti_classification) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            isoform_col = hdr.index("isoform")
            category_col = hdr.index("structural_category")
            gene_col = hdr.index("associated_gene")
            transcript_col = hdr.index("associated_transcript")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                sq[f[isoform_col]] = (f[category_col], f[gene_col], f[transcript_col])

    # novel models to write, in file order (file order keeps the renamed ids stable), with their class
    # code and the reference gene of their gffcompare match
    order = []
    cc_of = {}
    ref_gene_of = {}
    for c in gtf_records(args.merged, ("transcript",)):
        cc = attr(c[8], "class_code")
        tid = attr(c[8], "transcript_id")
        if cc is None or cc in drop or (keep is not None and tid not in keep):
            continue
        order.append(tid)
        cc_of[tid] = cc
        ref_gene_of[tid] = tx2gene.get(base(attr(c[8], "cmp_ref"))) or "NA"
    # gene id tama gave each kept model (used when SQANTI has no call)
    merged_gene = {}
    for c in gtf_records(args.merged, ("exon",)):
        tid = attr(c[8], "transcript_id")
        if tid in cc_of:
            merged_gene.setdefault(tid, attr(c[8], "gene_id"))

    # SQANTI novelGene label -> fresh locus id, numbered in order of first use
    locus_ids = {}

    def place(tid):
        # gene, assignment label and bridged genes of one model
        category, assoc_gene, assoc_transcript = sq.get(tid, (None, None, None))
        # no SQANTI call: keep tama's gene id
        if assoc_gene is None:
            return merged_gene.get(tid), "unclassified", []
        # a SQANTI novel gene: one new locus id per novelGene label
        if assoc_gene.startswith("novelGene"):
            locus_ids.setdefault(assoc_gene, "%sG_%06d" % (args.prefix, len(locus_ids) + 1))
            return locus_ids[assoc_gene], category, []
        # one or more known genes (joined by "_"): the gene of associated_transcript leads when it is one
        # of them and the model is not a fusion; otherwise SQANTI's first gene leads. the rest are bridged
        # genes
        genes = assoc_gene.split("_")
        lead = tx2gene.get(base(assoc_transcript))
        if category != "fusion" and lead in genes:
            genes.insert(0, genes.pop(genes.index(lead)))
        return genes[0], category, genes[1:]

    # one record per kept model: renamed <prefix>_NNNNNN and placed on its gene, or as is with --no-rename
    tx_map = {}
    for n, tid in enumerate(order, start=1):
        if args.no_rename:
            tx_map[tid] = NovelRecord(tid, merged_gene[tid], ref_gene_of[tid], cc_of[tid], "unclassified", [])
        else:
            gid, how, bridging = place(tid)
            tx_map[tid] = NovelRecord("%s_%06d" % (args.prefix, n), gid, ref_gene_of[tid], cc_of[tid],
                                      how, bridging)

    # provenance table: the new record plus the per-tool membership columns of the merged model
    if args.provenance_out and args.membership:
        with open(args.membership) as fh:
            # skip transcript_id, class_code, is_novel; carry the remaining columns (n_tools and one per tool)
            carry = fh.readline().rstrip("\n").split("\t")[MEMBERSHIP_FIRST_CARRIED_COLUMN:]
            member = {}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                member[f[0]] = f
        with open(args.provenance_out, "w") as out:
            out.write("transcript_id\tgene_id\tref_gene\tclass_code\tgene_assignment\tbridging\t%s\n"
                      % "\t".join(carry))
            for tid, rec in tx_map.items():
                if tid in member:
                    vals = member[tid][MEMBERSHIP_FIRST_CARRIED_COLUMN:]
                else:
                    vals = ["NA"] * len(carry)
                out.write("\t".join([rec.new_tid, str(rec.gene_id), rec.ref_gene, rec.class_code,
                                     rec.gene_assignment, ",".join(rec.bridging)] + vals) + "\n")

    # output gtf: the reference as is (unless --novel-only), then the kept models with new ids and tags
    with open(args.output, "w") as out:
        if not args.novel_only:
            with opener(args.reference) as ref:
                for line in ref:
                    out.write(line)
        for c in gtf_records(args.merged, ("transcript", "exon")):
            tid = attr(c[8], "transcript_id")
            if tid not in tx_map:
                continue
            rec = tx_map[tid]
            c[1] = args.source
            c[8] = 'gene_id "%s"; transcript_id "%s";' % (rec.gene_id, rec.new_tid)
            if c[2] == "transcript" and not args.no_rename:
                c[8] += (' ref_gene "%s"; class_code "%s"; gene_assignment "%s"; bridging "%s";'
                         % (rec.ref_gene, rec.class_code, rec.gene_assignment, ",".join(rec.bridging)))
            out.write("\t".join(c) + "\n")

    # summary: models on known genes (not a fresh <prefix>G_ locus) and models bridging a second gene
    known = sum(1 for rec in tx_map.values() if not str(rec.gene_id).startswith(args.prefix + "G_"))
    bridged = sum(1 for rec in tx_map.values() if rec.bridging)
    print("* reference + %d novel transcripts: %d are isoforms of known genes, %d new loci "
          "(sqanti novel genes), %d spanning a second gene (bridging); dropped codes: %s"
          % (len(order), known, len(locus_ids), bridged, ",".join(sorted(drop))))


if __name__ == "__main__":
    main()
