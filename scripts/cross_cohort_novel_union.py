#!/usr/bin/env python
"""Merge the cohorts' novel isoforms (source --novel-source) with tama_merge at the within-cohort
tolerances, drop single-exon novels, and append GENCODE verbatim. Gene ids are provisional;
cross_cohort_novel_genes.py sets the final ones. Cohort names carry no underscore (tama's separator)."""
import argparse
import gzip
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consensus_tama import read_bed12, run_tama, write_bed12  # noqa: E402


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_novel(path, cohort, source_tag, tx):
    n = 0
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
            uid = cohort + "_" + tid
            d = tx.get(uid)
            if d is None:
                d = {"chrom": c[0], "strand": c[6], "cohort": cohort, "orig_id": tid,
                     "gene_id": attr(c[8], "gene_id"), "class_code": attr(c[8], "class_code"),
                     "ref_gene": attr(c[8], "ref_gene"), "exons": []}
                tx[uid] = d
                n += 1
            if c[2] == "exon":
                d["exons"].append((int(c[3]), int(c[4])))
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--cohort", action="append", required=True, metavar="NAME=PATH")
    ap.add_argument("--novel-source", required=True)
    ap.add_argument("--wobble", type=int, required=True)
    ap.add_argument("--tes-tol", type=int, required=True)
    ap.add_argument("--ends", required=True, choices=["common_ends", "longest_ends"])
    ap.add_argument("--tama-merge", required=True)
    ap.add_argument("--python2", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--provenance-out", required=True)
    ap.add_argument("--summary-out", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    tx = {}
    per_cohort = {}
    for spec in args.cohort:
        name, path = spec.split("=", 1)
        if "_" in name:
            sys.exit("cohort name '%s' carries an underscore, which tama uses as its separator" % name)
        per_cohort[name] = load_novel(path, name, args.novel_source, tx)
    for d in tx.values():
        d["exons"].sort()
    n_loaded = len(tx)
    spliced = [u for u in tx if len(tx[u]["exons"]) > 1 and tx[u]["strand"] in ("+", "-")]
    n_mono = sum(1 for u in tx if len(tx[u]["exons"]) == 1)
    if len(spliced) + n_mono != n_loaded:
        sys.exit("%d spliced novel isoforms carry no strand; tama cannot place them"
                 % (n_loaded - n_mono - len(spliced)))
    print("* loaded novel isoforms: " + ", ".join("%s=%d" % (k, v) for k, v in per_cohort.items())
          + " (total %d, %d single-exon dropped)" % (n_loaded, n_mono))

    os.makedirs(args.workdir, exist_ok=True)
    filelist = os.path.join(args.workdir, "filelist.txt")
    with open(filelist, "w") as fl:
        for name in per_cohort:
            bed = os.path.join(args.workdir, name + ".bed")
            write_bed12(bed, [u for u in spliced if tx[u]["cohort"] == name], tx, lambda u: tx[u]["orig_id"])
            fl.write("%s\tno_cap\t1,1,1\t%s\n" % (bed, name))
    prefix = os.path.join(args.workdir, "tama")
    run_tama(args.python2, args.tama_merge, filelist, prefix, args.wobble, args.tes_tol, args.ends,
             len(spliced))
    models = read_bed12(prefix + ".bed")
    members = {}
    with open(prefix + "_trans_report.txt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i_all = hdr.index("all_source_trans")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            members[f[0]] = [m for m in f[i_all].split(",") if m]
    print("* tama merge: %d spliced novel -> %d cross-cohort isoforms (junctions +-%dbp, 3' +-%dbp, no_cap)"
          % (len(spliced), len(models), args.wobble, args.tes_tol))

    reps = []
    for mid, (chrom, strand, exons) in models.items():
        mem = members[mid]
        cohorts = sorted({m.split("_", 1)[0] for m in mem})
        mem = sorted(mem, key=lambda m: tx[m]["cohort"])
        longest = max(mem, key=lambda m: sum(e - s + 1 for s, e in tx[m]["exons"]))
        src_ids = "|".join("%s:%s" % (tx[m]["cohort"], tx[m]["orig_id"]) for m in mem)
        reps.append((mid, chrom, strand, exons, cohorts, src_ids, len(mem), tx[longest]))
    reps.sort(key=lambda r: (r[1], r[3][0][0], r[3][-1][1]))

    ref_tx = 0
    with open(args.output, "w") as out, open(args.provenance_out, "w") as prov:
        prov.write("new_id\tgene_id\tclass_code\tn_cohorts\tcohorts\tn_members\tsrc_ids\n")
        with opener(args.reference) as ref:
            for line in ref:
                if line.startswith("#"):
                    continue
                out.write(line)
                if line.split("\t", 3)[2] == "transcript":
                    ref_tx += 1
        for i, (_, chrom, strand, exons, cohorts, src_ids, nmem, d) in enumerate(reps, start=1):
            nid = "%s_%08d" % (args.prefix, i)
            gid = d["gene_id"] or nid
            # per-cohort locus ids restart in every cohort
            if "NOVELG_" in gid:
                gid = "%s_%s" % (d["cohort"], gid)
            cc = d["class_code"] or "u"
            out.write("\t".join([chrom, args.novel_source, "transcript", str(exons[0][0]), str(exons[-1][1]),
                                 ".", strand, ".",
                                 'gene_id "%s"; transcript_id "%s"; ref_gene "%s"; class_code "%s"; cohorts "%s";'
                                 % (gid, nid, d["ref_gene"] or "-", cc, ",".join(cohorts))]) + "\n")
            for s, e in exons:
                out.write("\t".join([chrom, args.novel_source, "exon", str(s), str(e), ".", strand, ".",
                                     'gene_id "%s"; transcript_id "%s";' % (gid, nid)]) + "\n")
            prov.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\n"
                       % (nid, gid, cc, len(cohorts), ",".join(cohorts), nmem, src_ids))

    shared = sum(1 for r in reps if len(r[4]) > 1)
    print("* wrote %d cross-cohort novel isoforms (%d shared by >1 cohort) + GENCODE -> %s"
          % (len(reps), shared, args.output))

    by_combo = Counter(",".join(r[4]) for r in reps)
    rows = [("novel_in_%s" % name, n) for name, n in per_cohort.items()]
    rows += [("novel_in_total", n_loaded), ("novel_monoexon_dropped", n_mono),
             ("cross_cohort_novel_after_endaware_collapse", len(reps)),
             ("collapsed_as_identical", len(spliced) - len(reps))]
    rows += [("by_combo:%s" % combo, n) for combo, n in sorted(by_combo.items(), key=lambda kv: (-kv[1], kv[0]))]
    rows += [("single_cohort_only", len(reps) - shared), ("shared_by_2plus_cohorts", shared),
             ("gencode_transcripts", ref_tx), ("final_total_transcripts", ref_tx + len(reps))]
    with open(args.summary_out, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.writelines("%s\t%d\n" % r for r in rows)
    print("* summary -> %s" % args.summary_out)


if __name__ == "__main__":
    main()
