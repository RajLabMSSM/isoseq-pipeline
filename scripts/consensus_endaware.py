#!/usr/bin/env python
"""
consensus_endaware.py

End-aware cross-tool consensus collapse. Unlike gffcompare (which merges transcripts
by intron chain and DISCARDS 5'/3' ends), this keeps alternative-TSS/TES (APA)
isoforms distinct: two transcripts collapse into one only when their intron chains
agree within +/- junction_wobble bp AND their 5' ends agree within tss_tol bp AND
their 3' ends agree within tes_tol bp. Single-exon transcripts collapse on the two
ends alone. This is the discovery catalog for publication, so it errs toward
retention; SQANTI (downstream) is the precision filter, and LAPA (downstream) is the
internal-priming-aware APA quantifier.

Inputs:
  --union      prefixed union GTF from build_union.py (carries src_tool/src_group)
  --gffcmp-prefix   gffcompare -o prefix; we read <prefix>.*.tmap for per-transcript
                    class_code + matched reference id (cmp_ref) vs the reference
Outputs:
  -o           merged consensus GTF: one representative per collapse cluster, novel
               only (class_code not in --exclude-codes), carrying gene_id (locus),
               class_code, cmp_ref -- the exact attributes build_reference_with_novel
               consumes.
  --overlap-out / --membership-out   Venn + per-isoform tool provenance (same formats
               as consensus_from_gffcompare.py so downstream R is unchanged).

Clustering avoids transitive over-chaining on the loose end tolerances by ANCHORING:
each cluster's first member is its anchor and every other member must fall within
tolerance of THAT anchor (not merely of its neighbour), so a run of APA sites spaced
just under the tolerance cannot chain into one blob.

Brooke Friedman
"""
import argparse
import glob
import re
from collections import defaultdict


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def base(tx):
    return tx.split(".")[0] if tx else tx


def load_union(path):
    """uid -> dict(chrom, strand, exons[(s,e)...] sorted, tool, group, lines[raw])."""
    tx = {}
    order = []
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
                     "tool": attr(c[8], "src_tool"), "group": attr(c[8], "src_group"),
                     "tline": None, "elines": []}
                tx[uid] = d
                order.append(uid)
            if c[2] == "transcript":
                d["tline"] = c
            else:
                d["exons"].append((int(c[3]), int(c[4])))
                d["elines"].append(c)
    for uid, d in tx.items():
        d["exons"].sort()
        ex = d["exons"]
        d["introns"] = tuple(x for i in range(len(ex) - 1) for x in (ex[i][1], ex[i + 1][0]))
        lo, hi = ex[0][0], ex[-1][1]
        d["span"] = (lo, hi)
        d["len"] = sum(e - s + 1 for s, e in ex)
        # 5' = TSS, 3' = TES, by strand
        if d["strand"] == "-":
            d["tss"], d["tes"] = hi, lo
        else:
            d["tss"], d["tes"] = lo, hi
    return tx, order


def load_tmap(prefix):
    """qry transcript_id -> (class_code, ref_id) from gffcompare <prefix>.*.tmap."""
    files = glob.glob(prefix + ".*.tmap")
    cc = {}
    for f in files:
        with open(f) as fh:
            header = fh.readline().split("\t")
            idx = {name.strip(): i for i, name in enumerate(header)}
            i_ref = idx.get("ref_id", 1)
            i_cc = idx.get("class_code", 2)
            i_q = idx.get("qry_id", 4)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) <= max(i_ref, i_cc, i_q):
                    continue
                cc[p[i_q]] = (p[i_cc], p[i_ref])
    return cc


def introns_match(a, b, wobble):
    if len(a) != len(b):
        return False
    return all(abs(x - y) <= wobble for x, y in zip(a, b))


def cluster(uids, tx, wobble, tss_tol, tes_tol):
    """Anchor-bounded greedy clustering within (chrom, strand, n_introns) buckets."""
    buckets = defaultdict(list)
    for uid in uids:
        d = tx[uid]
        buckets[(d["chrom"], d["strand"], len(d["introns"]))].append(uid)
    clusters = []
    for key, members in buckets.items():
        members.sort(key=lambda u: (tx[u]["introns"], tx[u]["tss"], tx[u]["tes"], u))
        anchors = []   # list of (anchor_uid, [members])
        for uid in members:
            d = tx[uid]
            placed = False
            for a_uid, grp in anchors:
                a = tx[a_uid]
                if (abs(d["tss"] - a["tss"]) <= tss_tol
                        and abs(d["tes"] - a["tes"]) <= tes_tol
                        and introns_match(d["introns"], a["introns"], wobble)):
                    grp.append(uid)
                    placed = True
                    break
            if not placed:
                anchors.append((uid, [uid]))
        for _, grp in anchors:
            clusters.append(grp)
    return clusters


def support(members, tx, min_tools):
    """(passes, tools_present_pooled). Per-group tool count, max across groups."""
    per_group = defaultdict(set)
    pooled = set()
    for uid in members:
        d = tx[uid]
        per_group[d["group"]].add(d["tool"])
        pooled.add(d["tool"])
    best = max((len(s) for s in per_group.values()), default=0)
    return best >= min_tools, pooled


def choose_rep(members, tx):
    return max(members, key=lambda u: (tx[u]["len"], u))


def assign_loci(reps, tx, prefix="XLOC"):
    """Group representatives into loci by same-strand genomic overlap (single-linkage
    sweep). Returns uid -> locus_id. Only used downstream for intergenic novels."""
    loci = {}
    by_cs = defaultdict(list)
    for uid in reps:
        d = tx[uid]
        by_cs[(d["chrom"], d["strand"])].append(uid)
    n = 0
    for _, uids in by_cs.items():
        uids.sort(key=lambda u: tx[u]["span"])
        cur = []
        cur_end = None
        for uid in uids:
            s, e = tx[uid]["span"]
            if cur and s <= cur_end:
                cur.append(uid)
                cur_end = max(cur_end, e)
            else:
                if cur:
                    n += 1
                    lid = "%s_%06d" % (prefix, n)
                    for u in cur:
                        loci[u] = lid
                cur = [uid]
                cur_end = e
        if cur:
            n += 1
            lid = "%s_%06d" % (prefix, n)
            for u in cur:
                loci[u] = lid
    return loci


TOOLS = ["bambu", "isoquant", "stringtie"]


def write_overlap(rows, out_tsv):
    from itertools import combinations
    counts = defaultdict(int)
    total = 0
    for r in rows:
        total += 1
        counts[tuple(sorted(r["tools"]))] += 1
    subsets = []
    for k in range(len(TOOLS), 0, -1):
        subsets.extend(combinations(TOOLS, k))
    with open(out_tsv, "w") as out:
        out.write("# NOVEL consensus isoforms by supporting-tool combination "
                  "(end-aware collapse: introns+-wobble, 5'+-tss_tol, 3'+-tes_tol)\n")
        out.write("combination\tn_isoforms\n")
        for s in subsets:
            out.write("%s\t%d\n" % ("+".join(s), counts.get(tuple(s), 0)))
        out.write("total_novel\t%d\n" % total)


def write_membership(rows, out_tsv):
    with open(out_tsv, "w") as out:
        out.write("transcript_id\tclass_code\tis_novel\tn_tools\t%s\n" % "\t".join(TOOLS))
        for r in rows:
            flags = "\t".join("1" if t in r["tools"] else "0" for t in TOOLS)
            out.write("%s\t%s\t1\t%d\t%s\n" % (r["tid"], r["cc"], len(r["tools"]), flags))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union", required=True)
    ap.add_argument("--gffcmp-prefix", required=True)
    ap.add_argument("--min-tools", type=int, default=1)
    ap.add_argument("--junction-wobble", type=int, default=6)
    ap.add_argument("--tss-tol", type=int, default=50)
    ap.add_argument("--tes-tol", type=int, default=50)
    ap.add_argument("--exclude-codes", default="=,c",
                    help="class codes treated as already-in-reference -> not novel")
    ap.add_argument("--overlap-out")
    ap.add_argument("--membership-out")
    ap.add_argument("--stats-out", help="TSV: per-tool raw + novel-candidate counts "
                    "(pre-collapse), for tracking_numbers")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    drop = set(c for c in args.exclude_codes.split(",") if c)
    tx, _ = load_union(args.union)
    tmap = load_tmap(args.gffcmp_prefix)

    novel = [uid for uid in tx if tmap.get(uid, ("u", "-"))[0] not in drop]

    if args.stats_out:
        raw_all = defaultdict(int)
        novel_cand = defaultdict(int)
        for uid, d in tx.items():
            raw_all[d["tool"]] += 1
        for uid in novel:
            novel_cand[tx[uid]["tool"]] += 1
        with open(args.stats_out, "w") as sf:
            sf.write("tool\traw_all\tnovel_candidates\n")
            for t in TOOLS:
                sf.write("%s\t%d\t%d\n" % (t, raw_all.get(t, 0), novel_cand.get(t, 0)))

    clusters = cluster(novel, tx, args.junction_wobble, args.tss_tol, args.tes_tol)

    kept = []   # (rep_uid, members, tools)
    for members in clusters:
        ok, tools = support(members, tx, args.min_tools)
        if not ok:
            continue
        rep = choose_rep(members, tx)
        kept.append((rep, members, tools))

    reps = [r for r, _, _ in kept]
    loci = assign_loci(reps, tx)

    rows = []
    with open(args.output, "w") as out:
        for i, (rep, members, tools) in enumerate(sorted(kept, key=lambda k: (
                tx[k[0]]["chrom"], tx[k[0]]["span"])), start=1):
            tid = "TCONS_%08d" % i
            cc, ref_id = tmap.get(rep, ("u", "-"))
            gid = loci.get(rep, "XLOC_%06d" % i)
            d = tx[rep]
            tline = list(d["tline"])
            tline[8] = ('gene_id "%s"; transcript_id "%s"; class_code "%s"; cmp_ref "%s";'
                        % (gid, tid, cc, ref_id))
            out.write("\t".join(tline) + "\n")
            for el in d["elines"]:
                e = list(el)
                e[8] = 'gene_id "%s"; transcript_id "%s";' % (gid, tid)
                out.write("\t".join(e) + "\n")
            rows.append({"tid": tid, "cc": cc, "tools": tools})

    if args.overlap_out:
        write_overlap(rows, args.overlap_out)
    if args.membership_out:
        write_membership(rows, args.membership_out)

    print("* end-aware collapse: %d novel candidates -> %d clusters kept "
          "(>= %d tools; introns +-%dbp, 5' +-%dbp, 3' +-%dbp) -> %s"
          % (len(novel), len(rows), args.min_tools, args.junction_wobble,
             args.tss_tol, args.tes_tol, args.output))


if __name__ == "__main__":
    main()
