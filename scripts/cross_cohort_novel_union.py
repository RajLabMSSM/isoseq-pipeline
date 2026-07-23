#!/usr/bin/env python
"""
cross_cohort_novel_union.py

Combine the SQANTI-passed NOVEL isoforms from >=2 cohorts into ONE cross-cohort set
using the SAME end-aware collapse as the within-cohort consensus (consensus_endaware:
two transcripts merge only when intron chains agree within +-junction_wobble bp AND
5' ends within +-tss_tol bp AND 3' ends within +-tes_tol bp; single-exon collapses on
the two ends alone). Everything else is kept -- alternative-TSS/TES (APA) isoforms and
near-but-not-identical structures stay distinct. Then GENCODE (verbatim) is appended so
the output is a full reference+novel GTF, mirroring build_reference_with_novel.

Novel input = transcripts whose source column is the novel source tag (default ONT) in
each cohort's *_reference_plus_novel.gtf(.gz). Provenance (which cohorts contributed) is
recorded on each surviving isoform as cohorts "a,b" plus src_ids.

Usage:
  cross_cohort_novel_union.py --reference gencode.v48.gtf \\
     --cohort tdpkd=.../tdpkd_nanopore_reference_plus_novel.gtf.gz \\
     --cohort sun=.../sun_..._reference_plus_novel.gtf.gz \\
     --cohort tanaka=.../tanaka_nanopore_reference_plus_novel.gtf.gz \\
     --novel-source ONT --wobble 6 --tss-tol 50 --tes-tol 50 \\
     --prefix XNOVEL -o cross_cohort_reference_plus_novel_v48.gtf
"""
import argparse
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consensus_endaware import introns_match, cluster, choose_rep  # noqa: E402


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_novel(path, cohort, source_tag, tx, order):
    """Add this cohort's novel transcripts (source col == source_tag) to tx."""
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[1] != source_tag or c[2] not in ("transcript", "exon"):
                continue
            tid = attr(c[8], "transcript_id")
            if tid is None:
                continue
            uid = cohort + "::" + tid
            d = tx.get(uid)
            if d is None:
                d = {"chrom": c[0], "strand": c[6], "cohort": cohort, "orig_id": tid,
                     "gene_id": attr(c[8], "gene_id"), "class_code": attr(c[8], "class_code"),
                     "ref_gene": attr(c[8], "ref_gene"), "exons": [], "tline": None, "elines": []}
                tx[uid] = d
                order.append(uid)
            if c[2] == "transcript":
                d["tline"] = c
            else:
                d["exons"].append((int(c[3]), int(c[4])))
                d["elines"].append(c)


def finalize(tx):
    for d in tx.values():
        d["exons"].sort()
        ex = d["exons"]
        d["introns"] = tuple(x for i in range(len(ex) - 1) for x in (ex[i][1], ex[i + 1][0]))
        lo, hi = ex[0][0], ex[-1][1]
        d["span"] = (lo, hi)
        d["len"] = sum(e - s + 1 for s, e in ex)
        d["tss"], d["tes"] = (hi, lo) if d["strand"] == "-" else (lo, hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, help="GENCODE GTF (.gtf or .gtf.gz), appended verbatim")
    ap.add_argument("--cohort", action="append", required=True, metavar="NAME=PATH",
                    help="cohort novel GTF (reference_plus_novel); repeat per cohort")
    ap.add_argument("--novel-source", default="ONT", help="GTF source column tag of novel records")
    ap.add_argument("--wobble", type=int, default=6)
    ap.add_argument("--tss-tol", type=int, default=50)
    ap.add_argument("--tes-tol", type=int, default=50)
    ap.add_argument("--prefix", default="XNOVEL")
    ap.add_argument("--provenance-out")
    ap.add_argument("--summary-out", help="TSV funnel (tracking_numbers_final_numbers): "
                    "per-cohort input, end-aware collapse, by-combination, + GENCODE totals")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    tx, order = {}, []
    per_cohort = {}
    for spec in args.cohort:
        name, path = spec.split("=", 1)
        before = len(order)
        load_novel(path, name, args.novel_source, tx, order)
        per_cohort[name] = len(order) - before
    finalize(tx)
    print("* loaded novel isoforms: " + ", ".join("%s=%d" % (k, v) for k, v in per_cohort.items())
          + " (total %d)" % len(order))

    clusters = cluster(order, tx, args.wobble, args.tss_tol, args.tes_tol)
    print("* end-aware collapse: %d novel -> %d cross-cohort isoforms "
          "(introns +-%dbp, 5' +-%dbp, 3' +-%dbp)"
          % (len(order), len(clusters), args.wobble, args.tss_tol, args.tes_tol))

    reps = []
    for members in clusters:
        rep = choose_rep(members, tx)
        cohorts = sorted({tx[u]["cohort"] for u in members})
        src_ids = "|".join("%s:%s" % (tx[u]["cohort"], tx[u]["orig_id"])
                           for u in sorted(members, key=lambda u: tx[u]["cohort"]))
        reps.append((rep, cohorts, src_ids, len(members)))

    reps.sort(key=lambda r: (tx[r[0]]["chrom"], tx[r[0]]["span"]))

    prov = open(args.provenance_out, "w") if args.provenance_out else None
    if prov:
        prov.write("new_id\tgene_id\tclass_code\tn_cohorts\tcohorts\tn_members\tsrc_ids\n")

    ref_tx = 0
    with open(args.output, "w") as out:
        with opener(args.reference) as ref:
            for line in ref:
                if line.startswith("#"):
                    continue
                out.write(line)
                if line.split("\t", 3)[2] == "transcript":
                    ref_tx += 1
        for i, (rep, cohorts, src_ids, nmem) in enumerate(reps, start=1):
            d = tx[rep]
            nid = "%s_%08d" % (args.prefix, i)
            gid = d["gene_id"] or nid
            cc = d["class_code"] or "u"
            extra = ('gene_id "%s"; transcript_id "%s"; ref_gene "%s"; class_code "%s"; '
                     'cohorts "%s";' % (gid, nid, d["ref_gene"] or "-", cc, ",".join(cohorts)))
            tl = list(d["tline"]); tl[8] = extra
            out.write("\t".join(tl) + "\n")
            for el in d["elines"]:
                e = list(el)
                e[8] = 'gene_id "%s"; transcript_id "%s";' % (gid, nid)
                out.write("\t".join(e) + "\n")
            if prov:
                prov.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\n"
                           % (nid, gid, cc, len(cohorts), ",".join(cohorts), nmem, src_ids))
    if prov:
        prov.close()

    shared = sum(1 for r in reps if len(r[1]) > 1)
    print("* wrote %d cross-cohort novel isoforms (%d shared by >1 cohort) + GENCODE -> %s"
          % (len(reps), shared, args.output))

    if args.summary_out:
        from collections import Counter
        by_combo = Counter(",".join(r[1]) for r in reps)
        with open(args.summary_out, "w") as s:
            s.write("metric\tvalue\n")
            for name, n in per_cohort.items():
                s.write("novel_in_%s\t%d\n" % (name, n))
            s.write("novel_in_total\t%d\n" % len(order))
            s.write("cross_cohort_novel_after_endaware_collapse\t%d\n" % len(reps))
            s.write("collapsed_as_identical\t%d\n" % (len(order) - len(reps)))
            for combo, n in sorted(by_combo.items(), key=lambda kv: (-kv[1], kv[0])):
                s.write("by_combo:%s\t%d\n" % (combo, n))
            s.write("single_cohort_only\t%d\n" % (len(reps) - shared))
            s.write("shared_by_2plus_cohorts\t%d\n" % shared)
            s.write("gencode_transcripts\t%d\n" % ref_tx)
            s.write("final_total_transcripts\t%d\n" % (ref_tx + len(reps)))
        print("* summary -> %s" % args.summary_out)


if __name__ == "__main__":
    main()
