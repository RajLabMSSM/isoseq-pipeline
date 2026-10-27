#!/usr/bin/env python
"""Sanity check of the final reference.

<out>             per gene: in_reference, total and novel isoform counts
<out>.events.tsv  per diagnostic event, traced through union -> consensus -> SQANTI -> final gtf.
                  Events are exons, GENE:chrom:start-end (1-based). RECOVERED needs an assembler
                  to have built it, since the final gtf holds the whole reference."""
import argparse
import gzip
import re
from collections import Counter


def attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def base(g):
    return g.split(".")[0] if g else g


def parse_events(spec):
    events = []
    for item in filter(None, (x.strip() for x in spec.split(","))):
        parts = item.split(":")
        start, end = parts[2].split("-")
        events.append({"gene": parts[0], "chrom": parts[1], "start": int(start), "end": int(end)})
    return events


def event_hits(path, events, tol):
    """per event: transcript_id -> source of every transcript carrying it"""
    pad = 5000
    windows = {}
    for ev in events:
        lo, hi = ev["start"] - pad, ev["end"] + pad
        w = windows.get(ev["chrom"])
        windows[ev["chrom"]] = (min(w[0], lo), max(w[1], hi)) if w else (lo, hi)
    tx = {}
    with opener(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(c) < 9 or c[2] != "exon" or c[0] not in windows:
                continue
            s, e = int(c[3]), int(c[4])
            w = windows[c[0]]
            if e < w[0] or s > w[1]:
                continue
            t = tx.setdefault(attr(c[8], "transcript_id") or "?", {"chrom": c[0], "source": c[1], "exons": []})
            t["exons"].append((s, e))
    hits = [{} for _ in events]
    for tid, t in tx.items():
        for i, ev in enumerate(events):
            if t["chrom"] == ev["chrom"] and any(
                    abs(a - ev["start"]) <= tol and abs(b - ev["end"]) <= tol for a, b in t["exons"]):
                hits[i][tid] = t["source"]
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--source", required=True)
    ap.add_argument("--events", required=True)
    ap.add_argument("--event-tol", type=int, required=True)
    ap.add_argument("--union-gtf", required=True)
    ap.add_argument("--consensus-gtf", required=True)
    ap.add_argument("--sqanti-reasons", required=True)
    ap.add_argument("--sqanti-pass-ids", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    sym2gid = {}
    with opener(args.reference) as fh:
        for line in fh:
            c = line.split("\t")
            if line.startswith("#") or len(c) < 9 or c[2] not in ("gene", "transcript"):
                continue
            name, gid = attr(c[8], "gene_name"), attr(c[8], "gene_id")
            if name and gid:
                sym2gid.setdefault(name, set()).add(base(gid))

    total, novel = Counter(), Counter()
    with opener(args.gtf) as fh:
        for line in fh:
            c = line.split("\t")
            if line.startswith("#") or len(c) < 9 or c[2] != "transcript":
                continue
            gid = base(attr(c[8], "gene_id"))
            total[gid] += 1
            if (attr(c[8], "transcript_id") or "").startswith(args.prefix) or c[1] == args.source:
                novel[gid] += 1

    rows, warns = [], []
    for sym in filter(None, (g.strip() for g in args.genes.split(","))):
        gids = sorted(sym2gid.get(sym, []))
        if not gids:
            rows.append((sym, "NA", "no", 0, 0))
            warns.append("%s: symbol not found in reference" % sym)
        rows += [(sym, gid, "yes" if total[gid] else "no", total[gid], novel[gid]) for gid in gids]
    with open(args.output, "w") as out:
        out.write("gene\tgene_id\tin_reference\tn_total_isoforms\tn_novel_isoforms\n")
        out.writelines("%s\t%s\t%s\t%d\t%d\n" % r for r in rows)
    print("* gene sanity check -> %s" % args.output)
    for r in rows:
        print("  %-8s %-20s total=%-4d novel=%-4d" % (r[0], r[1], r[3], r[4]))

    events = parse_events(args.events)
    in_ref, in_union, in_cons, in_final = (event_hits(p, events, args.event_tol) for p in
                                           (args.reference, args.union_gtf, args.consensus_gtf, args.gtf))
    with opener(args.sqanti_reasons) as fh:
        idx = fh.readline().rstrip("\n").split("\t").index("filter_reason")
        why = {c[0]: c[idx] for c in (line.rstrip("\n").split("\t") for line in fh) if len(c) > idx}
    with open(args.sqanti_pass_ids) as fh:
        passed = {line.strip() for line in fh if line.strip()}

    epath = args.output + ".events.tsv"
    with open(epath, "w") as out:
        out.write("gene\tevent\tin_gencode\tassembled_by\tn_consensus\tsqanti_kept\tsqanti_dropped\t"
                  "drop_reason\tin_final_gtf\tfinal_via\tverdict\n")
        for i, ev in enumerate(events):
            loc = "%s:%d-%d" % (ev["chrom"], ev["start"], ev["end"])
            tools = sorted({t.split("__")[0] if "__" in t else src.lower() for t, src in in_union[i].items()})
            kept = [t for t in in_cons[i] if not passed or t in passed]
            dropped = [t for t in in_cons[i] if passed and t not in passed]
            reasons = sorted({why[t] for t in dropped if t in why})
            via = ",".join(filter(None, [
                "novel" if any(t.startswith(args.prefix) for t in in_final[i]) else "",
                "reference" if any(not t.startswith(args.prefix) for t in in_final[i]) else ""])) or "-"
            if not in_union[i]:
                verdict = "NOT_ASSEMBLED"
            elif in_final[i]:
                verdict = "RECOVERED"
            elif dropped:
                verdict = "LOST_AT_SQANTI"
            else:
                verdict = "LOST_AT_CONSENSUS"
            out.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n" % (
                ev["gene"], loc, "yes" if in_ref[i] else "no", ",".join(tools) or "-",
                len(in_cons[i]), len(kept), len(dropped), ";".join(reasons) or "-",
                "yes" if in_final[i] else "no", via, verdict))

    print("* event trace -> %s" % epath)
    with open(epath) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            print("  %-8s %-28s %s" % (c[0], c[1], c[-1]))
    if warns:
        print("  WARN: " + "; ".join(warns))


if __name__ == "__main__":
    main()
