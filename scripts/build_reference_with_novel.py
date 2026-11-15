#!/usr/bin/env python
"""Combine the reference GTF with the novel models of the merged consensus.

The reference GTF is copied unchanged. The novel models added after it are those whose class code is not
in --exclude-codes and, when --keep-list is given, whose id is on that list.
Each novel model is renamed <prefix>_NNNNNN and gets the source <source>.

Each model's gene is SQANTI3's associated_gene. A model SQANTI calls a novelGene gets a new
<prefix>G_NNNNNN locus.
When SQANTI names several genes, one gene leads and the others are recorded as bridged genes.
For a fusion, the lead is the gene that lies 5' on the model's strand. Otherwise the lead is the gene of
associated_transcript.

With --novel-only --no-rename, the script instead writes the SQANTI input: only the novel models, under
their merged ids."""
import argparse
import gzip
import re
from collections import namedtuple

# One novel model as written to the output.
# It holds the new id, the gene, the reference gene of its gffcompare match and the class code.
# It also records how the gene was assigned (the SQANTI structural category) and any further genes it bridges.
# Under --no-rename the assignment is only a placeholder, because that mode never writes it.
NovelRecord = namedtuple("NovelRecord",
                         ["new_tid", "gene_id", "ref_gene", "class_code", "gene_assignment", "bridging"])
# The membership table starts with transcript_id, class_code and is_novel.
# The columns after these are the per-tool columns that are copied into the provenance table.
MEMBERSHIP_FIRST_CARRIED_COLUMN = 3


def opener(path):
    # Open plain or gzipped text files the same way.
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def attr(field, key):
    # Return the value of one GTF attribute, written as key "value";, or None when the key is missing.
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def base(tx):
    # Return the transcript id without its ".version" suffix.
    return tx.split(".")[0] if tx else tx


def gtf_records(path, features):
    # Yield the split columns of every GTF line whose feature type is in features.
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

    # Read the class codes to leave out, and the optional list of model ids to keep.
    # The keep list holds the models that passed the SQANTI filter.
    drop = set(c for c in args.exclude_codes.split(",") if c)
    keep = None
    if args.keep_list:
        with open(args.keep_list) as fh:
            keep = {line.strip() for line in fh if line.strip()}

    # Map each reference transcript, without its version, to its gene.
    # Also record each gene's span, from its lowest transcript start to its highest transcript end.
    # The span shows which gene of a fusion lies 5'.
    tx2gene = {}
    gene_span = {}
    for c in gtf_records(args.reference, ("transcript",)):
        gene = attr(c[8], "gene_id")
        tx2gene[base(attr(c[8], "transcript_id"))] = gene
        start, end = int(c[3]), int(c[4])
        low, high = gene_span.get(gene, (start, end))
        gene_span[gene] = (min(low, start), max(high, end))

    # Read the SQANTI classification.
    # For each isoform, keep its structural category, associated gene and associated transcript.
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

    # Collect the novel models to write, in file order. File order keeps the new ids stable between runs.
    # Also record each model's class code, strand and the reference gene of its gffcompare match.
    order = []
    cc_of = {}
    ref_gene_of = {}
    strand_of = {}
    for c in gtf_records(args.merged, ("transcript",)):
        cc = attr(c[8], "class_code")
        tid = attr(c[8], "transcript_id")
        if cc is None or cc in drop or (keep is not None and tid not in keep):
            continue
        order.append(tid)
        cc_of[tid] = cc
        ref_gene_of[tid] = tx2gene.get(base(attr(c[8], "cmp_ref"))) or "NA"
        strand_of[tid] = c[6]
    # Record the gene id that TAMA gave each kept model. Under --no-rename this id is written unchanged.
    merged_gene = {}
    for c in gtf_records(args.merged, ("exon",)):
        tid = attr(c[8], "transcript_id")
        if tid in cc_of:
            merged_gene.setdefault(tid, attr(c[8], "gene_id"))

    # Give each SQANTI novelGene label a new locus id, numbered in order of first use.
    locus_ids = {}

    def five_prime_gene(genes, strand):
        # Return the gene that starts first in the model's reading direction.
        # On the plus strand this is the gene with the lowest start. On the minus strand it is the gene with
        # the highest end.
        # SQANTI sorts genes by id, so a tie goes to the lowest id.
        # Fusion partners share the model's strand. This was checked on v9, where none differ.
        if strand == "+":
            return min(genes, key=lambda gene: gene_span[gene][0])
        return max(genes, key=lambda gene: gene_span[gene][1])

    def place(tid):
        # Decide the gene, the assignment label and the bridged genes of one model.
        # Every kept model has a SQANTI call, so a missing call stops the run instead of keeping an
        # unassigned gene.
        category, assoc_gene, assoc_transcript = sq[tid]
        # A SQANTI novel gene gets one new locus id per novelGene label.
        if assoc_gene.startswith("novelGene"):
            locus_ids.setdefault(assoc_gene, "%sG_%06d" % (args.prefix, len(locus_ids) + 1))
            return locus_ids[assoc_gene], category, []
        # SQANTI can name one or more known genes, joined by "_". A fusion is led by its 5' gene.
        # Any other model is led by the gene of associated_transcript when that gene is one of them.
        # Otherwise SQANTI's first gene leads. The remaining genes are recorded as bridged.
        genes = assoc_gene.split("_")
        if category == "fusion":
            lead = five_prime_gene(genes, strand_of[tid])
        else:
            lead = tx2gene.get(base(assoc_transcript))
        if lead in genes:
            genes.insert(0, genes.pop(genes.index(lead)))
        return genes[0], category, genes[1:]

    # Build one record per kept model. Each model is renamed <prefix>_NNNNNN and placed on its gene.
    # With --no-rename the model keeps its merged id and TAMA gene.
    tx_map = {}
    for n, tid in enumerate(order, start=1):
        if args.no_rename:
            tx_map[tid] = NovelRecord(tid, merged_gene[tid], ref_gene_of[tid], cc_of[tid], "unclassified", [])
        else:
            gid, how, bridging = place(tid)
            tx_map[tid] = NovelRecord("%s_%06d" % (args.prefix, n), gid, ref_gene_of[tid], cc_of[tid],
                                      how, bridging)

    # Write the provenance table: each new record plus the per-tool membership columns of its merged model.
    if args.provenance_out and args.membership:
        with open(args.membership) as fh:
            # Skip transcript_id, class_code and is_novel.
            # Carry the remaining columns, which are n_tools and one column per tool.
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

    # Write the output GTF. The reference comes first, unchanged, unless --novel-only is set.
    # The kept models follow with their new ids and tags.
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

    # Print a summary of how many models sit on known genes rather than a new <prefix>G_ locus.
    # Also count the models that bridge a second gene.
    # The summary is skipped under --no-rename, because no genes have been assigned yet.
    if not args.no_rename:
        known = sum(1 for rec in tx_map.values() if not str(rec.gene_id).startswith(args.prefix + "G_"))
        bridged = sum(1 for rec in tx_map.values() if rec.bridging)
        print("* reference + %d novel transcripts: %d are isoforms of known genes, %d new loci "
              "(sqanti novel genes), %d spanning a second gene (bridging); dropped codes: %s"
              % (len(order), known, len(locus_ids), bridged, ",".join(sorted(drop))))


if __name__ == "__main__":
    main()
