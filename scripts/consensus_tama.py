#!/usr/bin/env python
"""Cross-tool consensus with tama_merge (Kuo et al. 2020): the spliced, stranded novel models of the
union (gffcompare class not in --exclude-codes) go in as one no_cap source per tool and group; every
model from any one tool is kept. The merged models are re-classified against the reference."""
import argparse
import glob
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
    # union gtf -> {union id: chrom, strand, sorted exons, source tool and group}
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
                d = {"chrom": c[0], "strand": c[6], "exons": [],
                     "tool": attr(c[8], "src_tool"), "group": attr(c[8], "src_group")}
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


def write_bed12(path, uids, tx, name=lambda uid: uid):
    # tama_merge input: bed12, 0-based start, sorted by chrom and start; the name field is
    # "<gene>;<transcript>" (tama splits on ";", here the gene part is a dummy "g"), score 40, colour 0,0,0
    with open(path, "w") as out:
        for uid in sorted(uids, key=lambda u: (tx[u]["chrom"], tx[u]["exons"][0][0])):
            d = tx[uid]
            ex = d["exons"]
            start, end = ex[0][0] - 1, ex[-1][1]
            sizes = ",".join(str(e - s + 1) for s, e in ex)
            starts = ",".join(str(s - 1 - start) for s, e in ex)
            out.write("\t".join([d["chrom"], str(start), str(end), "g;" + name(uid), "40", d["strand"],
                                 str(start), str(end), "0,0,0", str(len(ex)), sizes, starts]) + "\n")


def read_bed12(path):
    # tama_merge output bed12 -> {merged model id: (chrom, strand, 1-based exons)}; the id is the part of the
    # name field after ";"
    models = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            start = int(c[1])
            sizes = [int(x) for x in c[10].rstrip(",").split(",")]
            starts = [int(x) for x in c[11].rstrip(",").split(",")]
            exons = [(start + s + 1, start + s + z) for s, z in zip(starts, sizes)]
            mid = c[3].split(";")[1]
            models[mid] = (c[0], c[5], exons)
    return models


def gtf_lines(mid, chrom, strand, exons, extra=""):
    # one transcript line plus its exon lines; tama names models <gene>.<n>, so the gene id is the id
    # minus its last ".<n>"
    ids = 'gene_id "%s"; transcript_id "%s";' % (mid.rsplit(".", 1)[0], mid)
    rows = [[chrom, "TAMA", "transcript", str(exons[0][0]), str(exons[-1][1]), ".", strand, ".", ids + extra]]
    rows += [[chrom, "TAMA", "exon", str(s), str(e), ".", strand, ".", ids] for s, e in exons]
    return "".join("\t".join(r) + "\n" for r in rows)


def run_tama(python2, tama_merge, filelist, prefix, wobble, tes_tol, n_input):
    # tama_merge prints errors to stdout and still exits 0
    log = prefix + ".log"
    with open(log, "w") as lf:
        subprocess.run([python2, tama_merge, "-f", filelist, "-p", prefix,
                        "-m", str(wobble), "-z", str(tes_tol),
                        "-d", "merge_dup"], check=True, stdout=lf, stderr=subprocess.STDOUT)
    with open(log) as lf:
        text = lf.read()
    if "Error" in text or "TAMA Merge has completed successfully" not in text:
        raise SystemExit("tama_merge failed, see %s:\n%s" % (log, text[-2000:]))
    # a source model can be listed under two merged models, so count distinct members
    members = set()
    with open(prefix + "_merge.txt") as mf:
        for line in mf:
            members.add(line.split("\t")[3].split(";")[1])
    if len(members) != n_input:
        raise SystemExit("tama_merge consumed %d of %d input transcripts, see %s" % (len(members), n_input, log))


def read_members(report):
    # tama trans_report -> {merged model id: union ids merged into it}; members are <source>_<union id>
    members = {}
    with open(report) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i_all = hdr.index("all_source_trans")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            members[f[0]] = [m.split("_", 1)[1] for m in f[i_all].split(",") if m]
    return members


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--union", required=True)
    ap.add_argument("--gffcmp-prefix", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--junction-wobble", type=int, required=True)
    ap.add_argument("--tes-tol", type=int, required=True)
    ap.add_argument("--exclude-codes", required=True)
    ap.add_argument("--tama-merge", required=True)
    ap.add_argument("--python2", required=True)
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

    # novel candidates = class not excluded; split into single-exon, unstranded and spliced stranded models
    candidates = [uid for uid in tx if tmap.get(uid, UNMATCHED)[0] not in drop]
    mono = [uid for uid in candidates if len(tx[uid]["exons"]) == 1]
    # tama cannot place an unstranded model
    unstranded = [uid for uid in candidates if len(tx[uid]["exons"]) > 1 and tx[uid]["strand"] not in ("+", "-")]
    novel = [uid for uid in candidates if len(tx[uid]["exons"]) > 1 and tx[uid]["strand"] in ("+", "-")]

    # per-tool counts at each step, for the stats table
    raw_by_tool = Counter(tx[u]["tool"] for u in tx)
    candidates_by_tool = Counter(tx[u]["tool"] for u in candidates)
    mono_by_tool = Counter(tx[u]["tool"] for u in mono)
    unstranded_by_tool = Counter(tx[u]["tool"] for u in unstranded)
    with open(args.stats_out, "w") as sf:
        sf.write("tool\traw_all\tnovel_candidates\tmonoexon_dropped\tunstranded_dropped\n")
        for t in TOOLS:
            sf.write("%s\t%d\t%d\t%d\t%d\n" % (t, raw_by_tool[t], candidates_by_tool[t],
                                                  mono_by_tool[t], unstranded_by_tool[t]))

    # one bed per tool + group source. tama joins source name and transcript with "_", so a source name may not
    # contain one: keep only letters and digits. two different groups could then end up with the same name
    # (cohort_1CTRL and cohort1_CTRL both become cohort1CTRL), so stop if that happens instead of merging them
    os.makedirs(args.workdir, exist_ok=True)
    source_of = {}
    for d in tx.values():
        source_of[(d["tool"], d["group"])] = re.sub(r"[^A-Za-z0-9]", "", d["tool"] + d["group"])
    if len(set(source_of.values())) != len(source_of):
        raise SystemExit("two tool + group pairs give the same tama source name: %s" % source_of)
    by_source = defaultdict(list)
    for uid in novel:
        d = tx[uid]
        by_source[source_of[(d["tool"], d["group"])]].append(uid)
    filelist = os.path.join(args.workdir, "filelist.txt")
    with open(filelist, "w") as fl:
        for src in sorted(by_source):
            bed = os.path.join(args.workdir, src + ".bed")
            write_bed12(bed, by_source[src], tx)
            # tama filelist row: bed, cap flag (no_cap = 5' ends not trusted), priority for
            # 5' end, junctions, 3' end (1,1,1 = all sources equal), source name
            fl.write("%s\tno_cap\t1,1,1\t%s\n" % (bed, src))

    # merge with tama and read back the merged models and their members
    prefix = os.path.join(args.workdir, "tama")
    run_tama(args.python2, args.tama_merge, filelist, prefix, args.junction_wobble, args.tes_tol,
             len(novel))

    models = read_bed12(prefix + ".bed")
    members = read_members(prefix + "_trans_report.txt")

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
        chrom, _strand, exons = models[mid]
        return (chrom, exons[0][0], exons[-1][1])

    # every group (condition) in this run, for the per-group support columns
    group_names = sorted({d["group"] for d in tx.values()})

    # write the consensus gtf with class code and reference tags, and note which tools and which groups built
    # each model (from the union models tama merged into it)
    rows = []
    order = sorted(models, key=model_position)
    with open(args.output, "w") as out:
        for mid in order:
            cc, ref_id = merged_cc.get(mid, UNMATCHED)
            tools = sorted({tx[u]["tool"] for u in members.get(mid, [])})
            groups = sorted({tx[u]["group"] for u in members.get(mid, [])})
            out.write(gtf_lines(mid, *models[mid], extra=' class_code "%s"; cmp_ref "%s";' % (cc, ref_id)))
            rows.append({"transcript_id": mid, "class_code": cc, "tools": tools, "groups": groups})

    # count merged models per supporting-tool combination, all three tools first
    counts = Counter(tuple(r["tools"]) for r in rows)
    with open(args.overlap_out, "w") as out:
        out.write("# NOVEL consensus isoforms by supporting-tool combination "
                  "(tama_merge: junctions +-wobble, 3' +-tes_tol, no_cap)\n")
        out.write("combination\tn_isoforms\n")
        for k in range(len(TOOLS), 0, -1):
            for s in combinations(TOOLS, k):
                out.write("%s\t%d\n" % ("+".join(s), counts[s]))
        out.write("total_novel\t%d\n" % len(rows))
    # membership table: one row per merged model, a 0/1 column per tool, then the number of groups and a 0/1
    # column per group (group_<name>), so a model built in one condition only can be told apart. the group
    # columns come last: readers pick columns by name, and build_reference_with_novel.py carries every column
    # after is_novel into the provenance table
    with open(args.membership_out, "w") as out:
        out.write("transcript_id\tclass_code\tis_novel\tn_tools\t%s\tn_groups\t%s\n"
                  % ("\t".join(TOOLS), "\t".join("group_" + g for g in group_names)))
        for r in rows:
            tool_flags = "\t".join("1" if t in r["tools"] else "0" for t in TOOLS)
            group_flags = "\t".join("1" if g in r["groups"] else "0" for g in group_names)
            out.write("%s\t%s\t1\t%d\t%s\t%d\t%s\n" % (r["transcript_id"], r["class_code"], len(r["tools"]),
                                                        tool_flags, len(r["groups"]), group_flags))

    print("* tama merge: %d novel candidates, %d single-exon and %d unstranded dropped, %d spliced "
          "-> %d merged models (junctions +-%dbp, 3' +-%dbp, no_cap) -> %s"
          % (len(candidates), len(mono), len(unstranded), len(novel), len(rows),
             args.junction_wobble, args.tes_tol, args.output))


if __name__ == "__main__":
    main()
