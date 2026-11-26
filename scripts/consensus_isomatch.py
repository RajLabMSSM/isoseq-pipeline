#!/usr/bin/env python
"""Merge the novel models of all tools and cohorts with isomatch merge 0.8.0 (Zheng et al. 2026, bioRxiv,
doi:10.64898/2026.10.09.758114).

The spliced, stranded novel models of the union (gffcompare class not in --exclude-codes) go in as one GTF per
tool, cohort and group, and every model built by any one source is kept. isomatch runs at its defaults, except that
start sites are not used to split models. The merged models are then classified against the reference with
gffcompare (Pertea and Pertea 2020). Each merged model records which tools, cohorts and groups built it, read from
isomatch's ISOM_SRC attribute."""
import argparse
import glob
import gzip
import os
import re
import subprocess
from collections import Counter, defaultdict
from itertools import combinations

TOOLS = ["bambu", "isoquant", "stringtie"]
# class code and reference id given to a model gffcompare did not report (treated as intergenic, no reference)
UNMATCHED = ("u", "-")


def attr(field, key):
    # value of one gtf attribute (key "value";), or None when the key is missing
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def load_union(path):
    # Read the union GTF into {union id: chrom, strand, sorted exons, source tool, cohort, group and gene}.
    # A group name such as TDP occurs in more than one cohort, so the group is stored as <cohort>_<group>.
    tx = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("transcript", "exon"):
                continue
            uid = attr(c[8], "transcript_id")
            if uid is None:
                continue
            d = tx.get(uid)
            if d is None:
                cohort = attr(c[8], "src_cohort")
                d = {"chrom": c[0], "strand": c[6], "exons": [], "tool": attr(c[8], "src_tool"),
                     "cohort": cohort, "group": "%s_%s" % (cohort, attr(c[8], "src_group")),
                     "gene": attr(c[8], "gene_id") or uid}
                tx[uid] = d
            if c[2] == "exon":
                d["exons"].append((int(c[3]), int(c[4])))
    for d in tx.values():
        d["exons"].sort()
    return tx


def load_tmap(prefix):
    # gffcompare .tmap -> {query id: (class code, reference id)}; header names first, fixed columns as fallback
    cc = {}
    for f in glob.glob(prefix + ".*.tmap"):
        with open(f) as fh:
            header = fh.readline().split("\t")
            idx = {name.strip(): i for i, name in enumerate(header)}
            i_ref, i_cc, i_q = idx.get("ref_id", 1), idx.get("class_code", 2), idx.get("qry_id", 4)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) > max(i_ref, i_cc, i_q):
                    cc[p[i_q]] = (p[i_cc], p[i_ref])
    return cc


def write_gtf(path, uids, tx):
    # one isomatch input per source, so ISOM_SRC names the source of every merged member
    with open(path, "w") as out:
        for uid in sorted(uids, key=lambda u: (tx[u]["chrom"], tx[u]["exons"][0][0])):
            d = tx[uid]
            out.write(gtf_lines(uid, d["gene"], d["chrom"], d["strand"], d["exons"], source="union"))


def gtf_lines(mid, gene, chrom, strand, exons, extra="", source="isomatch"):
    ids = 'gene_id "%s"; transcript_id "%s";' % (gene, mid)
    rows = [[chrom, source, "transcript", str(exons[0][0]), str(exons[-1][1]), ".", strand, ".", ids + extra]]
    rows += [[chrom, source, "exon", str(s), str(e), ".", strand, ".", ids] for s, e in exons]
    return "".join("\t".join(r) + "\n" for r in rows)


def run_isomatch(isomatch, genome, prefix, gtfs):
    # long-read 5' ends are unreliable after rna degradation and truncated reverse transcription, so start sites do not split models
    log = prefix + ".log"
    with open(log, "w") as lf:
        result = subprocess.run([isomatch, "merge", "--ref-fa", genome, "-o", prefix,
                                 "--terminal-refine", "tes", "--terminal-refine-nc", "tes"] + gtfs,
                                stdout=lf, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        with open(log) as lf:
            raise SystemExit("isomatch merge failed, see %s:\n%s" % (log, lf.read()[-2000:]))


def read_isomatch(path):
    # the union ids of each merged model come from ISOM_SRC, whose records are S<input>:<union id>:... joined by "|"
    models, members, exons = {}, {}, defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(c) < 9:
                continue
            mid = attr(c[8], "transcript_id")
            if c[2] == "exon":
                exons[mid].append((int(c[3]), int(c[4])))
            elif c[2] == "transcript":
                if mid in models:
                    raise SystemExit("isomatch wrote transcript id %s twice" % mid)
                models[mid] = (attr(c[8], "gene_id"), c[0], c[6])
                members[mid] = [rec.split(":", 2)[1] for rec in attr(c[8], "ISOM_SRC").split("|")]
    return {mid: m + (sorted(exons[mid]),) for mid, m in models.items()}, members


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--union", required=True)
    ap.add_argument("--gffcmp-prefix", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--exclude-codes", required=True)
    ap.add_argument("--isomatch", required=True)
    ap.add_argument("--gffcompare", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--overlap-out", required=True)
    ap.add_argument("--membership-out", required=True)
    ap.add_argument("--stats-out", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    # load the union models and their gffcompare class against the reference
    drop = set(c for c in args.exclude_codes.split(",") if c)
    tx = load_union(args.union)
    tmap = load_tmap(args.gffcmp_prefix)

    # novel candidates = class not excluded; split into single-exon and spliced models
    # unstranded spliced models go into the merge too, and sqanti3 later sets them aside because it cannot classify them
    candidates = [uid for uid in tx if tmap.get(uid, UNMATCHED)[0] not in drop]
    mono = [uid for uid in candidates if len(tx[uid]["exons"]) == 1]
    novel = [uid for uid in candidates if len(tx[uid]["exons"]) > 1]

    # Count the models of each cohort and tool at each step, for the stats table. raw_all counts the union, which
    # already lacks the transcripts that carry a GENCODE id.
    cohort_names = sorted({d["cohort"] for d in tx.values()})

    def count(uids):
        return Counter((tx[u]["cohort"], tx[u]["tool"]) for u in uids)

    raw, cand, mono_n = count(tx), count(candidates), count(mono)
    with open(args.stats_out, "w") as sf:
        sf.write("cohort\ttool\traw_all\tnovel_candidates\tmonoexon_dropped\n")
        for c in cohort_names:
            for t in TOOLS:
                k = (c, t)
                sf.write("%s\t%s\t%d\t%d\t%d\n" % (c, t, raw[k], cand[k], mono_n[k]))

    # one input gtf per tool, cohort and group, written to the workdir because isomatch indexes next to its inputs
    os.makedirs(args.workdir, exist_ok=True)
    by_source = defaultdict(list)
    for uid in novel:
        by_source[(tx[uid]["tool"], tx[uid]["group"])].append(uid)
    gtfs = []
    for tool, group in sorted(by_source):
        path = os.path.join(args.workdir, "%s__%s.gtf" % (tool, group))
        write_gtf(path, by_source[(tool, group)], tx)
        gtfs.append(path)

    # merge with isomatch and read back the merged models and their members
    prefix = os.path.join(args.workdir, "isomatch")
    run_isomatch(args.isomatch, args.genome, prefix, gtfs)
    models, members = read_isomatch(prefix + ".merged.gtf.gz")

    # every spliced candidate must sit in exactly one merged model, or provenance would be lost or double counted
    consumed = Counter(u for mem in members.values() for u in mem)
    if set(consumed) != set(novel) or max(consumed.values()) > 1:
        raise SystemExit("isomatch merged %d distinct of %d input transcripts, %d listed twice"
                         % (len(set(consumed) & set(novel)), len(novel), sum(n > 1 for n in consumed.values())))

    # re-classify the merged models against the reference with gffcompare
    merged_gtf = os.path.join(args.workdir, "merged.gtf")
    with open(merged_gtf, "w") as out:
        for mid, model in models.items():
            out.write(gtf_lines(mid, *model))
    cmp_prefix = os.path.join(args.workdir, "mergedcmp")
    subprocess.run([args.gffcompare, "-r", args.reference, "-o", cmp_prefix, merged_gtf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    merged_cc = load_tmap(cmp_prefix)

    # sort the merged models by chrom, first exon start, last exon end
    def model_position(mid):
        _gene, chrom, _strand, exons = models[mid]
        return (chrom, exons[0][0], exons[-1][1])

    # Every group (condition) of every cohort in this run, for the per-group support columns.
    group_names = sorted({d["group"] for d in tx.values()})

    # Write the consensus GTF with the class code, the reference tag and the cohorts that built each model. The
    # tools, cohorts and groups of a model come from the union models that isomatch merged into it. The cohorts are
    # written in alphabetical order and joined by a comma, as in the published v9 GTF (sun, tdpkd or sun,tdpkd).
    rows = []
    order = sorted(models, key=model_position)
    with open(args.output, "w") as out:
        for mid in order:
            cc, ref_id = merged_cc.get(mid, UNMATCHED)
            mem = members[mid]
            tools = sorted({tx[u]["tool"] for u in mem})
            cohorts = sorted({tx[u]["cohort"] for u in mem})
            groups = sorted({tx[u]["group"] for u in mem})
            cohort_tools = {(tx[u]["cohort"], tx[u]["tool"]) for u in mem}
            out.write(gtf_lines(mid, *models[mid], extra=' class_code "%s"; cmp_ref "%s"; cohorts "%s";'
                                % (cc, ref_id, ",".join(cohorts))))
            rows.append({"transcript_id": mid, "class_code": cc, "tools": tools, "cohorts": cohorts,
                         "groups": groups, "cohort_tools": cohort_tools})

    # count merged models per supporting-tool combination, all three tools first
    counts = Counter(tuple(r["tools"]) for r in rows)
    with open(args.overlap_out, "w") as out:
        out.write("# NOVEL consensus isoforms by supporting-tool combination "
                  "(isomatch merge: default settings, terminal refine on 3' ends only)\n")
        out.write("combination\tn_isoforms\n")
        for k in range(len(TOOLS), 0, -1):
            for s in combinations(TOOLS, k):
                out.write("%s\t%d\n" % ("+".join(s), counts[s]))
        out.write("total_novel\t%d\n" % len(rows))
    # Write the membership table with one row per merged model. It holds a 0/1 column per tool, the number and
    # names of the cohorts, a 0/1 column per cohort and tool (<cohort>_<tool>), and the number of groups with a
    # 0/1 column per group (group_<cohort>_<group>). A model built in one cohort or one condition only can so be
    # told apart. Readers pick columns by name, and build_reference_with_novel.py carries every column after
    # is_novel into the provenance table.
    cohort_tool_cols = [(c, t) for c in cohort_names for t in TOOLS]
    with open(args.membership_out, "w") as out:
        out.write("transcript_id\tclass_code\tis_novel\tn_tools\t%s\tn_cohorts\tcohorts\t%s\tn_groups\t%s\n"
                  % ("\t".join(TOOLS), "\t".join("%s_%s" % ct for ct in cohort_tool_cols),
                     "\t".join("group_" + g for g in group_names)))
        for r in rows:
            tool_flags = "\t".join("1" if t in r["tools"] else "0" for t in TOOLS)
            cohort_tool_flags = "\t".join("1" if ct in r["cohort_tools"] else "0" for ct in cohort_tool_cols)
            group_flags = "\t".join("1" if g in r["groups"] else "0" for g in group_names)
            out.write("%s\t%s\t1\t%d\t%s\t%d\t%s\t%s\t%d\t%s\n"
                      % (r["transcript_id"], r["class_code"], len(r["tools"]), tool_flags, len(r["cohorts"]),
                         ",".join(r["cohorts"]), cohort_tool_flags, len(r["groups"]), group_flags))

    print("* isomatch merge: %d novel candidates, %d single-exon dropped, %d spliced "
          "-> %d merged models -> %s"
          % (len(candidates), len(mono), len(novel), len(rows), args.output))


if __name__ == "__main__":
    main()
