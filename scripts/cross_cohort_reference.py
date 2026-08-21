#!/usr/bin/env python
"""
cross_cohort_reference.py

Build ONE cross-cohort reference = full GENCODE (verbatim) + the UNION of novel
isoforms from two or more cohorts, with only EXACT structural duplicates collapsed.
Goal: keep as many isoforms as possible when combining cohorts (e.g. knockdown +
overexpression) without exact duplicates -- so matching is gffcompare's exact
intron-chain collapse (NOT wobble; wobble is used only WITHIN a cohort to pool the
3 assemblers, where tools disagree by a few bp on identical alignments). Across
cohorts a real novel's consensus junctions tend to converge, and the goal is maximal
retention, so near-identical novels are deliberately KEPT as separate entries.

How: run gffcompare -r GENCODE over the per-cohort `_merged_consensus.gtf` files.
gffcompare collapses identical intron chains into one transfrag, assigns a class_code
vs GENCODE, and records per-query (per-cohort) presence in `.tracking` -> that gives
provenance for free. Novels = transcripts whose class_code is not in --exclude-codes
(default '=,c'). Gene assignment mirrors build_reference_with_novel.py: a novel
overlapping a known gene (has cmp_ref) inherits that gene_id (novel ISOFORM of a
known gene); intergenic novels get a fresh <prefix>G_ locus. Each novel carries a
`cohorts "tdpkd,sun"` provenance attribute.

Usage:
  cross_cohort_reference.py --reference gencode.gtf \\
      --cohort tdpkd=.../tdpkd_nanopore_merged_consensus.gtf \\
      --cohort sun=.../sun_..._merged_consensus.gtf \\
      -o cross_cohort.isoforms.gtf
"""
import argparse
import os
import re
import subprocess
import tempfile

GFFCOMPARE = "/sc/arion/projects/ad-omics/data/software/gffcompare-0.12.10.Linux_x86_64/gffcompare"


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def base(tx):
    return tx.split(".")[0] if tx else tx


def load_reference_tx2gene(ref_path):
    tx2gene = {}
    with open(ref_path) as ref:
        for line in ref:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "transcript":
                continue
            tx2gene[base(attr(c[8], "transcript_id"))] = attr(c[8], "gene_id")
    return tx2gene


def load_provenance(tracking_path, cohort_names):
    """TCONS id -> sorted list of cohorts whose query column is present."""
    prov = {}
    with open(tracking_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            tcons = f[0].split("|")[0]
            cohorts = [cohort_names[i] for i, q in enumerate(f[4:])
                       if i < len(cohort_names) and q.strip() != "-"]
            prov[tcons] = sorted(set(cohorts))
    return prov


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, help="GENCODE GTF, passed through verbatim")
    ap.add_argument("--cohort", action="append", required=True, metavar="NAME=GTF",
                    help="cohort merged_consensus GTF; repeatable. ORDER sets provenance")
    ap.add_argument("--exclude-codes", default="=,c",
                    help="class codes already in the reference -> dropped (default '=,c')")
    ap.add_argument("--prefix", default="NOVEL")
    ap.add_argument("--source", default="ONT")
    ap.add_argument("--gffcompare", default=GFFCOMPARE)
    ap.add_argument("--provenance-out",
                    help="TSV of which novels are from which cohort(s) "
                         "(default <output>.provenance.tsv)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    drop = set(c for c in args.exclude_codes.split(",") if c)

    names, gtfs = [], []
    for spec in args.cohort:
        name, path = spec.split("=", 1)
        names.append(name)
        gtfs.append(path)

    # exact intron-chain collapse across cohorts + class_code vs GENCODE + provenance
    workdir = tempfile.mkdtemp(prefix="xcohort_")
    gpref = os.path.join(workdir, "xcohort")
    subprocess.run([args.gffcompare, "-r", args.reference, "-o", gpref] + gtfs, check=True)
    combined = gpref + ".combined.gtf"
    tracking = gpref + ".tracking"

    tx2gene = load_reference_tx2gene(args.reference)
    prov = load_provenance(tracking, names)

    # pass 1: assign ids to novel transfrags
    tx_map = {}
    n_tx = n_locus = 0
    with open(combined) as fh:
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
            n_tx += 1
            parent = tx2gene.get(base(attr(c[8], "cmp_ref")))
            if parent:
                gid, ref_gene = parent, parent
            else:
                n_locus += 1
                gid, ref_gene = "%sG_%06d" % (args.prefix, n_locus), "NA"
            tx_map[tid] = ("%s_%06d" % (args.prefix, n_tx), gid, ref_gene, cc)

    # pass 2: GENCODE verbatim, then renamed novel records with provenance
    with open(args.output, "w") as out:
        with open(args.reference) as ref:
            for line in ref:
                out.write(line)
        with open(combined) as fh:
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
                c[8] = ('gene_id "%s"; transcript_id "%s"; ref_gene "%s"; class_code "%s"; cohorts "%s";'
                        % (gid, new_tid, ref_gene, cc, ",".join(prov.get(tid, []))))
                out.write("\t".join(c) + "\n")

    # provenance: which cohort(s) each novel came from (knockdown / overexpression / both)
    prov_out = args.provenance_out or (args.output + ".provenance.tsv")
    tally = {}
    with open(prov_out, "w") as out:
        out.write("transcript_id\tgene_id\tref_gene\tclass_code\tcohorts\tcategory\n")
        for tid, (new_tid, gid, ref_gene, cc) in tx_map.items():
            cohorts = prov.get(tid, [])
            category = "both" if len(cohorts) > 1 else (cohorts[0] if cohorts else "NA")
            tally[category] = tally.get(category, 0) + 1
            out.write("%s\t%s\t%s\t%s\t%s\t%s\n"
                      % (new_tid, gid, ref_gene, cc, ",".join(cohorts), category))

    multi = sum(1 for v in tx_map if len(prov.get(v, [])) > 1)
    print("* reference + %d novel transcripts (exact-dedup across cohorts): %d shared by "
          ">1 cohort, %d new intergenic loci -> %s"
          % (n_tx, multi, n_locus, args.output))
    print("* novel provenance -> %s : %s"
          % (prov_out, ", ".join("%s=%d" % (k, v) for k, v in sorted(tally.items()))))


if __name__ == "__main__":
    main()
