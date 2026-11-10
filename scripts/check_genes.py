#!/usr/bin/env python
"""Sanity check of the final reference.

<out>             per gene: in_reference, total and novel isoform counts
<out>.events.tsv  per diagnostic event, the isoforms of the final gtf (GENCODE plus the novels that passed
                  SQANTI) that carry it, traced through union -> consensus -> SQANTI -> final gtf.
                  Events are junctions, GENE:chrom:start-end (first and last intron base, 1-based); several
                  junctions joined by + must all be in one isoform. RECOVERED needs an assembler
                  to have built it, since the final gtf holds the whole reference."""
import argparse
import gzip
import re
from collections import Counter

# bases added on each side of the events when choosing which gtf lines to keep in memory
WINDOW_PAD = 5000


def attr(field, key):
    # value of one gtf attribute (key "value";), or None when the key is missing
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def opener(path):
    # open plain or gzipped text the same way
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def base(g):
    # id without its ".version" suffix
    return g.split(".")[0] if g else g


def parse_events(spec):
    # "GENE:chrom:start-end+start-end,..." -> list of event dicts; empty items are skipped
    events = []
    items = [x.strip() for x in spec.split(",")]
    for item in items:
        if not item:
            continue
        gene, chrom, joined = item.split(":")
        # one (start, end) per junction; a cryptic exon inside an intron has one junction on each side
        junctions = []
        for junction in joined.split("+"):
            start, end = junction.split("-")
            junctions.append((int(start), int(end)))
        events.append({"gene": gene, "chrom": chrom, "junctions": junctions, "label": chrom + ":" + joined})
    return events


def event_hits(path, events, tol):
    """per event: transcript_id -> source of every transcript carrying it"""
    # one window per chromosome spanning all of its events plus the pad
    windows = {}
    for ev in events:
        lo = min(start for start, end in ev["junctions"]) - WINDOW_PAD
        hi = max(end for start, end in ev["junctions"]) + WINDOW_PAD
        w = windows.get(ev["chrom"])
        if w:
            windows[ev["chrom"]] = (min(w[0], lo), max(w[1], hi))
        else:
            windows[ev["chrom"]] = (lo, hi)
    # exons of every transcript with an exon inside a window
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
    # a transcript carries an event when it has every junction of the event, both ends within tol
    hits = [{} for _ in events]
    for tid, t in tx.items():
        exons = sorted(t["exons"])
        # the intron between two neighboring exons runs from one exon's end + 1 to the next exon's start - 1
        introns = [(left[1] + 1, right[0] - 1) for left, right in zip(exons, exons[1:])]
        for i, ev in enumerate(events):
            if t["chrom"] == ev["chrom"] and all(
                    any(abs(a - start) <= tol and abs(b - end) <= tol for a, b in introns)
                    for start, end in ev["junctions"]):
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

    # gene symbol -> reference gene ids (no version)
    sym2gid = {}
    with opener(args.reference) as fh:
        for line in fh:
            c = line.split("\t")
            if line.startswith("#") or len(c) < 9 or c[2] not in ("gene", "transcript"):
                continue
            name, gid = attr(c[8], "gene_name"), attr(c[8], "gene_id")
            if name and gid:
                sym2gid.setdefault(name, set()).add(base(gid))

    # isoforms per gene in the final gtf, total and novel (novel = new id prefix or the novel source tag)
    total = Counter()
    novel = Counter()
    with opener(args.gtf) as fh:
        for line in fh:
            c = line.split("\t")
            if line.startswith("#") or len(c) < 9 or c[2] != "transcript":
                continue
            gid = base(attr(c[8], "gene_id"))
            total[gid] += 1
            if (attr(c[8], "transcript_id") or "").startswith(args.prefix) or c[1] == args.source:
                novel[gid] += 1

    # per-gene table: one row per gene id of each requested symbol, a placeholder row when the symbol is missing
    rows = []
    warns = []
    symbols = [g.strip() for g in args.genes.split(",")]
    for sym in symbols:
        if not sym:
            continue
        gids = sorted(sym2gid.get(sym, []))
        if not gids:
            rows.append((sym, "NA", "no", 0, 0))
            warns.append("%s: symbol not found in reference" % sym)
        for gid in gids:
            in_reference = "yes" if total[gid] else "no"
            rows.append((sym, gid, in_reference, total[gid], novel[gid]))
    with open(args.output, "w") as out:
        out.write("gene\tgene_id\tin_reference\tn_total_isoforms\tn_novel_isoforms\n")
        out.writelines("%s\t%s\t%s\t%d\t%d\n" % r for r in rows)
    print("* gene sanity check -> %s" % args.output)
    for r in rows:
        print("  %-8s %-20s total=%-4d novel=%-4d" % (r[0], r[1], r[3], r[4]))

    # find the transcripts carrying each event at every stage
    events = parse_events(args.events)
    in_ref = event_hits(args.reference, events, args.event_tol)
    in_union = event_hits(args.union_gtf, events, args.event_tol)
    in_cons = event_hits(args.consensus_gtf, events, args.event_tol)
    in_final = event_hits(args.gtf, events, args.event_tol)
    # SQANTI filter reason per dropped model, and the ids that passed
    why = {}
    with opener(args.sqanti_reasons) as fh:
        idx = fh.readline().rstrip("\n").split("\t").index("filter_reason")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) > idx:
                why[c[0]] = c[idx]
    with open(args.sqanti_pass_ids) as fh:
        passed = {line.strip() for line in fh if line.strip()}

    # event trace: one row per event, with the stage at which it was lost (or RECOVERED)
    epath = args.output + ".events.tsv"
    header = ["gene", "event", "in_gencode", "assembled_by", "n_consensus", "sqanti_kept", "sqanti_dropped",
              "drop_reason", "in_final_gtf", "final_via", "final_isoforms", "verdict"]
    # echoed to the log after writing: gene, event, verdict of each row (header included)
    log_rows = [(header[0], header[1], header[-1])]
    with open(epath, "w") as out:
        out.write("\t".join(header) + "\n")
        for i, ev in enumerate(events):
            loc = ev["label"]
            # assembler that built it: the <tool>__ prefix of union ids, else the gtf source column
            tools = sorted({t.split("__")[0] if "__" in t else src.lower() for t, src in in_union[i].items()})
            # consensus models carrying it that SQANTI kept or dropped
            kept = [t for t in in_cons[i] if t in passed]
            dropped = [t for t in in_cons[i] if t not in passed]
            reasons = sorted({why[t] for t in dropped if t in why})
            # in the final gtf: via a novel model, a reference model, or both
            via_parts = []
            if any(t.startswith(args.prefix) for t in in_final[i]):
                via_parts.append("novel")
            if any(not t.startswith(args.prefix) for t in in_final[i]):
                via_parts.append("reference")
            via = ",".join(via_parts) or "-"
            if not in_union[i]:
                verdict = "NOT_ASSEMBLED"
            elif in_final[i]:
                verdict = "RECOVERED"
            elif dropped:
                verdict = "LOST_AT_SQANTI"
            else:
                verdict = "LOST_AT_CONSENSUS"
            # the ids of the isoforms after SQANTI that carry the event
            final_ids = ",".join(sorted(in_final[i])) or "-"
            out.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n" % (
                ev["gene"], loc, "yes" if in_ref[i] else "no", ",".join(tools) or "-",
                len(in_cons[i]), len(kept), len(dropped), ";".join(reasons) or "-",
                "yes" if in_final[i] else "no", via, final_ids, verdict))
            log_rows.append((ev["gene"], loc, verdict))

    print("* event trace -> %s" % epath)
    for gene, loc, verdict in log_rows:
        print("  %-8s %-28s %s" % (gene, loc, verdict))
    if warns:
        print("  WARN: " + "; ".join(warns))


if __name__ == "__main__":
    main()
