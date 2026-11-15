#!/usr/bin/env python
"""Set the final gene_id of every cross-cohort novel isoform.

Each isoform takes the gene and bridging genes of its longest member, read from that member's cohort GTF.
When members are equally long, the first one in src_ids order wins.
Some longest members sit on a per-cohort novel locus (NOVELG_). These isoforms are grouped into loci
across all cohorts with gffcompare.
Each locus gets one <prefix>G_ id, numbered in gffcompare order."""
import argparse
import gzip
import os
import re
import subprocess

# Gene ids of the per-cohort novel loci carry this tag.
# build_reference_with_novel.py writes it, because the cohort runs use --prefix NOVEL.
COHORT_LOCUS_TAG = "NOVELG_"


def attr(field, key):
    # Return the value of one GTF attribute, written as key "value";, or None when the key is missing.
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def opener(path):
    # Open plain or gzipped text files the same way.
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_cohort(path, cohort, source, members):
    # Read the novel models of one cohort GTF into members, under the key <cohort>:<id>.
    # For each model, store its gene, its bridging genes and its summed exon length.
    with opener(path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 9 or c[1] != source:
                continue
            key = cohort + ":" + attr(c[8], "transcript_id")
            if c[2] == "transcript":
                members[key] = {"gene": attr(c[8], "gene_id"), "bridging": attr(c[8], "bridging") or "",
                                "length": 0}
            elif c[2] == "exon":
                members[key]["length"] += int(c[4]) - int(c[3]) + 1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--provenance", required=True)
    ap.add_argument("--cohort", action="append", required=True, metavar="NAME=PATH")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--gffcompare", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--novel-source", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--provenance-out", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    # Load every cohort's novel models, keyed the same way as the src_ids column of the provenance table.
    members = {}
    for spec in args.cohort:
        name, path = spec.split("=", 1)
        load_cohort(path, name, args.novel_source, members)

    # Read the union provenance table, which has one row per cross-cohort model.
    with open(args.provenance) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        prov = [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh]

    # Give each model the gene and bridging genes of its longest member.
    # Models whose longest member sits on a per-cohort novel locus get a cross-cohort locus id further down.
    gene = {}
    bridging = {}
    locus_tx = []
    for r in prov:
        rep = max(r["src_ids"].split("|"), key=lambda m: members[m]["length"])
        m = members[rep]
        if COHORT_LOCUS_TAG in m["gene"]:
            locus_tx.append(r["new_id"])
        else:
            gene[r["new_id"]] = m["gene"]
        bridging[r["new_id"]] = m["bridging"]

    # Group the novel-locus models into loci with gffcompare, using its XLOC ids.
    # Each locus then gets one <prefix>G_ id.
    os.makedirs(args.workdir, exist_ok=True)
    loci_gtf = os.path.join(args.workdir, "novel_loci.gtf")
    wanted = set(locus_tx)
    with opener(args.gtf) as fh, open(loci_gtf, "w") as out:
        for line in fh:
            c = line.split("\t")
            if len(c) >= 9 and c[1] == args.novel_source and attr(c[8], "transcript_id") in wanted:
                out.write(line)
    pfx = os.path.join(args.workdir, "loci")
    subprocess.run([args.gffcompare, "-r", args.reference, "-o", pfx, loci_gtf], check=True)
    xloc = {}
    with open(pfx + ".annotated.gtf") as fh:
        for line in fh:
            c = line.split("\t")
            if c[2] == "transcript":
                xloc[attr(c[8], "transcript_id")] = attr(c[8], "xloc")
    if set(xloc) != wanted:
        raise SystemExit("gffcompare placed %d of %d novel-locus isoforms" % (len(xloc), len(wanted)))
    locus_id = {x: "%sG_%06d" % (args.prefix, i)
                for i, x in enumerate(sorted(set(xloc.values())), start=1)}
    for t in locus_tx:
        gene[t] = locus_id[xloc[t]]

    # Rewrite the GTF. Reference lines are copied unchanged.
    # Novel lines get the final gene id, and novel transcript lines also get the bridging genes.
    # The old gene ids are kept so that the number of changed genes can be reported.
    old = {}
    with opener(args.gtf) as fh, open(args.output, "w") as out:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[1] != args.novel_source:
                out.write(line)
                continue
            tid = attr(c[8], "transcript_id")
            if c[2] == "transcript":
                old[tid] = attr(c[8], "gene_id")
            c[8] = re.sub(r'gene_id "[^"]*"', 'gene_id "%s"' % gene[tid], c[8], count=1)
            if c[2] == "transcript":
                c[8] = re.sub(r' bridging "[^"]*";', "", c[8]) + ' bridging "%s";' % bridging[tid]
            out.write("\t".join(c) + "\n")

    # Write the provenance table with the final gene id, and with the bridging genes as the last column.
    cols = [h for h in hdr if h != "bridging"] + ["bridging"]
    with open(args.provenance_out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in prov:
            r["gene_id"] = gene[r["new_id"]]
            r["bridging"] = bridging[r["new_id"]]
            out.write("\t".join(r[h] for h in cols) + "\n")

    changed = sum(1 for t in old if old[t] != gene[t])
    print("* %d novel isoforms, %d on %d cross-cohort novel loci, %d gene_ids changed -> %s"
          % (len(prov), len(locus_tx), len(locus_id), changed, args.output))


if __name__ == "__main__":
    main()
