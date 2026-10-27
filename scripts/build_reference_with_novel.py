#!/usr/bin/env python
"""The reference gtf verbatim plus the novel merged-consensus models (class code not in
--exclude-codes, and in --keep-list), renamed <prefix>_NNNNNN with source <source>.

The gene is SQANTI3's associated_gene; a SQANTI novelGene becomes a fresh <prefix>G_NNNNNN locus.
When several genes are named, the gene of associated_transcript leads (a fusion: the first gene),
the others go to bridging. --novel-only --no-rename writes the SQANTI input instead."""
import argparse
import gzip
import re


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def base(tx):
    return tx.split(".")[0] if tx else tx


def gtf_records(path, features):
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

    drop = set(c for c in args.exclude_codes.split(",") if c)
    keep = None
    if args.keep_list:
        with open(args.keep_list) as fh:
            keep = {line.strip() for line in fh if line.strip()}

    tx2gene = {base(attr(c[8], "transcript_id")): attr(c[8], "gene_id")
               for c in gtf_records(args.reference, ("transcript",))}

    sq = {}
    if args.sqanti_classification:
        with opener(args.sqanti_classification) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ii, ci, gi, ti = (hdr.index(k) for k in
                              ("isoform", "structural_category", "associated_gene", "associated_transcript"))
            for line in fh:
                f = line.rstrip("\n").split("\t")
                sq[f[ii]] = (f[ci], f[gi], f[ti])

    # file order keeps the renamed ids stable
    order, cc_of, cmp_of = [], {}, {}
    for c in gtf_records(args.merged, ("transcript",)):
        cc, tid = attr(c[8], "class_code"), attr(c[8], "transcript_id")
        if cc is None or cc in drop or (keep is not None and tid not in keep):
            continue
        order.append(tid)
        cc_of[tid] = cc
        cmp_of[tid] = tx2gene.get(base(attr(c[8], "cmp_ref"))) or "NA"
    merged_gene = {}
    for c in gtf_records(args.merged, ("exon",)):
        tid = attr(c[8], "transcript_id")
        if tid in cc_of:
            merged_gene.setdefault(tid, attr(c[8], "gene_id"))

    locus_ids = {}

    def place(tid):
        cat, ag, at = sq.get(tid, (None, None, None))
        if ag is None:
            return merged_gene.get(tid), "unclassified", []
        if ag.startswith("novelGene"):
            locus_ids.setdefault(ag, "%sG_%06d" % (args.prefix, len(locus_ids) + 1))
            return locus_ids[ag], cat, []
        genes = ag.split("_")
        lead = tx2gene.get(base(at))
        if cat != "fusion" and lead in genes:
            genes.insert(0, genes.pop(genes.index(lead)))
        return genes[0], cat, genes[1:]

    tx_map = {}
    for n, tid in enumerate(order, start=1):
        if args.no_rename:
            tx_map[tid] = (tid, merged_gene[tid], cmp_of[tid], cc_of[tid], "unclassified", [])
        else:
            gid, how, bridging = place(tid)
            tx_map[tid] = ("%s_%06d" % (args.prefix, n), gid, cmp_of[tid], cc_of[tid], how, bridging)

    if args.provenance_out and args.membership:
        with open(args.membership) as fh:
            carry = fh.readline().rstrip("\n").split("\t")[3:]
            member = {f[0]: f for f in (line.rstrip("\n").split("\t") for line in fh)}
        with open(args.provenance_out, "w") as out:
            out.write("transcript_id\tgene_id\tref_gene\tclass_code\tgene_assignment\tbridging\t%s\n"
                      % "\t".join(carry))
            for tid, (new_tid, gid, ref_gene, cc, how, bridging) in tx_map.items():
                vals = member[tid][3:] if tid in member else ["NA"] * len(carry)
                out.write("\t".join([new_tid, str(gid), ref_gene, cc, how, ",".join(bridging)] + vals) + "\n")

    with open(args.output, "w") as out:
        if not args.novel_only:
            with opener(args.reference) as ref:
                for line in ref:
                    out.write(line)
        for c in gtf_records(args.merged, ("transcript", "exon")):
            tid = attr(c[8], "transcript_id")
            if tid not in tx_map:
                continue
            new_tid, gid, ref_gene, cc, how, bridging = tx_map[tid]
            c[1] = args.source
            c[8] = 'gene_id "%s"; transcript_id "%s";' % (gid, new_tid)
            if c[2] == "transcript" and not args.no_rename:
                c[8] += (' ref_gene "%s"; class_code "%s"; gene_assignment "%s"; bridging "%s";'
                         % (ref_gene, cc, how, ",".join(bridging)))
            out.write("\t".join(c) + "\n")

    known = sum(1 for v in tx_map.values() if not str(v[1]).startswith(args.prefix + "G_"))
    bridged = sum(1 for v in tx_map.values() if v[5])
    print("* reference + %d novel transcripts: %d are isoforms of known genes, %d new loci "
          "(sqanti novel genes), %d spanning a second gene (bridging); dropped codes: %s"
          % (len(order), known, len(locus_ids), bridged, ",".join(sorted(drop))))


if __name__ == "__main__":
    main()
