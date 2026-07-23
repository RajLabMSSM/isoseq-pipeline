#!/usr/bin/env python
"""
check_genes.py -- sanity check that expected genes are in the final reference, and
trace each diagnostic cryptic-splicing EVENT through the pipeline.

Two reports are written.

1) <out>                  per-GENE: gene, gene_id, in_reference, n_total_isoforms,
                          n_novel_isoforms  (unchanged columns, back-compatible)

2) <out>.events.tsv       per-EVENT: for each diagnostic exon (e.g. the STMN2 or
                          UNC13A cryptic exon), does the exon appear in the FINAL
                          GTF at all -- as a NOVEL isoform or as a REFERENCE one --
                          and if not, at which stage was it lost?

Why (2) exists: counting NOVEL isoforms alone silently inverts meaning whenever
GENCODE catches up on a cryptic exon. GENCODE v50 annotates the STMN2 cryptic exon
(ENST00001131000 / STMN2-207), so a PERFECTLY recovered tSTMN2 now matches the
reference (class_code "="), is discarded as a known isoform, and scores 0 novel --
which reads as failure but is success. An event is RECOVERED if its exon is present
in the final GTF, whatever the provenance.

Stages traced per event (each optional; skipped if the input is not supplied):
  assembled   -- exon present in the gffcompare-annotated union (which tools)
  consensus   -- exon survived the end-aware collapse
  sqanti      -- transcript(s) carrying the exon passed the SQANTI filter
  final       -- exon present in reference_plus_novel.gtf (novel and/or reference)

Events are given as GENE:chrom:start-end (1-based, GTF coords), matched with a
tolerance (--event-tol, default 0) on both boundaries. Two event kinds:
  * exon inclusion (default)      -- recovered when the cryptic exon is PRESENT
                                     (STMN2, UNC13A). start-end = the exon.
  * cassette-exon skip (:skip)    -- recovered when the flanking exons are joined
                                     by an intron with the cassette exon GONE
                                     (KCNQ2 exon-5 skip). start-end = that intron
                                     span; syntax GENE:chrom:start-end:skip.
"""
import argparse
import gzip
from collections import defaultdict

import wobble

attr = wobble.parse_attr


def opener(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def base(g):
    return g.split(".")[0] if g else g


def parse_events(spec):
    events = []
    for item in (spec or "").split(","):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        gene, chrom, span = parts[0], parts[1], parts[2]
        etype = parts[3] if len(parts) > 3 else "exon"
        start, end = span.split("-")
        events.append({"gene": gene, "chrom": chrom,
                       "start": int(start), "end": int(end), "type": etype})
    return events


def event_hits(path, events, tol):
    """transcript_ids (and their source/class) carrying each event.
    exon events match a single exon; skip events (start-end = the intron span that
    replaces the skipped exon) match an intron between two consecutive exons, i.e. the
    upstream exon spliced straight to the downstream exon with the cassette exon gone."""
    hits = [defaultdict(dict) for _ in events]
    if not path:
        return hits
    pad = 5000
    windows = {}
    for ev in events:
        lo, hi = ev["start"] - pad, ev["end"] + pad
        w = windows.get(ev["chrom"])
        windows[ev["chrom"]] = (min(w[0], lo), max(w[1], hi)) if w else (lo, hi)
    tx = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "exon":
                continue
            w = windows.get(c[0])
            if not w:
                continue
            s, e = int(c[3]), int(c[4])
            if e < w[0] or s > w[1]:
                continue
            tid = attr(c[8], "transcript_id") or "?"
            t = tx.get(tid)
            if t is None:
                t = tx[tid] = {"chrom": c[0], "exons": [], "meta": {
                    "source": c[1],
                    "gene_id": attr(c[8], "gene_id") or "",
                    "class_code": attr(c[8], "class_code") or "",
                    "cmp_ref": attr(c[8], "cmp_ref") or "",
                    "ref_id": attr(c[8], "reference_id") or "",
                }}
            t["exons"].append((s, e))
    for tid, t in tx.items():
        exons = sorted(t["exons"])
        introns = [(exons[k][1] + 1, exons[k + 1][0] - 1)
                   for k in range(len(exons) - 1)]
        for i, ev in enumerate(events):
            if t["chrom"] != ev["chrom"]:
                continue
            spans = introns if ev.get("type") == "skip" else exons
            for a, b in spans:
                if abs(a - ev["start"]) <= tol and abs(b - ev["end"]) <= tol:
                    hits[i][tid] = t["meta"]
                    break
    return hits


def load_sqanti(path):
    """transcript_id -> (structural_category, RTS_stage, all_canonical)"""
    info = {}
    if not path:
        return info
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(head)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if not c or not c[0]:
                continue
            info[c[0]] = (
                c[idx["structural_category"]] if "structural_category" in idx else "",
                c[idx["RTS_stage"]] if "RTS_stage" in idx else "",
                c[idx["all_canonical"]] if "all_canonical" in idx else "",
            )
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, help="GENCODE GTF (symbol -> gene_id)")
    ap.add_argument("--gtf", required=True, help="reference_plus_novel GTF to check")
    ap.add_argument("--genes", required=True, help="comma-separated gene symbols")
    ap.add_argument("--prefix", default="NOVEL", help="novel transcript_id prefix")
    ap.add_argument("--source", default="ONT", help="novel GTF source column")
    ap.add_argument("--events", default="",
                    help="GENE:chrom:start-end[,...] diagnostic exons; append :skip for a "
                         "cassette-exon skip (start-end = the skip intron)")
    ap.add_argument("--event-tol", type=int, default=0, help="bp tolerance on exon bounds")
    ap.add_argument("--union-gtf", default="", help="gffcompare-annotated union GTF")
    ap.add_argument("--consensus-gtf", default="", help="merged_consensus GTF")
    ap.add_argument("--sqanti-classification", default="", help="SQANTI classification.txt")
    ap.add_argument("--sqanti-pass-ids", default="", help="one passing transcript_id per line")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    wanted = [g.strip() for g in args.genes.split(",") if g.strip()]

    sym2gid = {}
    with opener(args.reference) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] not in ("gene", "transcript"):
                continue
            name, gid = attr(c[8], "gene_name"), attr(c[8], "gene_id")
            if name and gid:
                sym2gid.setdefault(name, set()).add(base(gid))

    total = {}
    novel = {}
    with opener(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "transcript":
                continue
            gid = base(attr(c[8], "gene_id"))
            tid = attr(c[8], "transcript_id") or ""
            total[gid] = total.get(gid, 0) + 1
            if tid.startswith(args.prefix) or c[1] == args.source:
                novel[gid] = novel.get(gid, 0) + 1

    rows, warns = [], []
    for sym in wanted:
        gids = sorted(sym2gid.get(sym, []))
        if not gids:
            rows.append((sym, "NA", "no", 0, 0))
            warns.append("%s: symbol not found in reference" % sym)
            continue
        for gid in gids:
            nt, nv = total.get(gid, 0), novel.get(gid, 0)
            rows.append((sym, gid, "yes" if nt else "no", nt, nv))

    with open(args.output, "w") as out:
        out.write("gene\tgene_id\tin_reference\tn_total_isoforms\tn_novel_isoforms\n")
        for r in rows:
            out.write("%s\t%s\t%s\t%d\t%d\n" % r)

    print("* gene sanity check -> %s" % args.output)
    for r in rows:
        print("  %-8s %-20s total=%-4d novel=%-4d" % (r[0], r[1], r[3], r[4]))

    events = parse_events(args.events)
    if not events:
        if warns:
            print("  WARN: " + "; ".join(warns))
        return

    tol = args.event_tol
    in_ref = event_hits(args.reference, events, tol)
    in_union = event_hits(args.union_gtf, events, tol)
    in_cons = event_hits(args.consensus_gtf, events, tol)
    in_final = event_hits(args.gtf, events, tol)
    sq = load_sqanti(args.sqanti_classification)

    passed = set()
    if args.sqanti_pass_ids:
        with open(args.sqanti_pass_ids) as fh:
            passed = {l.strip() for l in fh if l.strip()}

    epath = args.output + ".events.tsv"
    with open(epath, "w") as out:
        out.write("gene\tevent\tin_gencode\tassembled_by\tn_consensus\t"
                  "sqanti_kept\tsqanti_dropped\tdrop_reason\t"
                  "in_final_gtf\tfinal_via\tverdict\n")
        for i, ev in enumerate(events):
            loc = "%s:%d-%d" % (ev["chrom"], ev["start"], ev["end"])
            if ev.get("type") == "skip":
                loc += ":skip"

            tools = sorted({t.split("__")[0] for t in in_union[i]
                            if "__" in t} | {v["source"].lower()
                                             for t, v in in_union[i].items()
                                             if "__" not in t})
            kept, dropped, reasons = [], [], set()
            for tid in in_cons[i]:
                if not passed or tid in passed:
                    kept.append(tid)
                else:
                    dropped.append(tid)
                    cat, rts, canon = sq.get(tid, ("", "", ""))
                    if rts == "TRUE":
                        reasons.add("RTS_stage=TRUE")
                    if canon and canon != "canonical":
                        reasons.add("non_canonical")
                    if not rts and not canon:
                        reasons.add("filtered")

            novel_final = [t for t in in_final[i] if t.startswith(args.prefix)]
            ref_final = [t for t in in_final[i] if not t.startswith(args.prefix)]
            via = ",".join(filter(None, [
                "novel" if novel_final else "",
                "reference" if ref_final else "",
            ])) or "-"

            recovered = bool(in_final[i])
            if recovered:
                verdict = "RECOVERED"
            elif dropped:
                verdict = "LOST_AT_SQANTI"
            elif in_union[i]:
                verdict = "LOST_AT_CONSENSUS"
            else:
                verdict = "NOT_ASSEMBLED"

            out.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n" % (
                ev["gene"], loc,
                "yes" if in_ref[i] else "no",
                ",".join(tools) if tools else "-",
                len(in_cons[i]), len(kept), len(dropped),
                ";".join(sorted(reasons)) if reasons else "-",
                "yes" if recovered else "no", via, verdict,
            ))

    print("* event trace -> %s" % epath)
    with open(epath) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            print("  %-8s %-28s %s" % (c[0], c[1], c[-1]))
    if warns:
        print("  WARN: " + "; ".join(warns))


if __name__ == "__main__":
    main()
