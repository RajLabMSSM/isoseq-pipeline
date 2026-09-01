#!/usr/bin/env python
"""
Concatenate the per-tool, per-group assembler gtfs into one union gtf, giving every
transcript a globally unique id and tagging its origin.

Every transcript is kept distinct here, rather than pre-merged by gffcompare, because
gffcompare merges on intron chain alone and drops the 5' and 3' ends. Keeping them lets
consensus_endaware.py apply its own end tolerance and so retain alternative polyA
isoforms (and alternative TSS for mono-exons, the only case where the 5' end is used).

Ids are prefixed "<tool>__<group>__" since they collide across tools and groups
(STRG.1.1 and BambuTx1 both recur), and src_tool / src_group attributes are appended so
the collapser can recover tool support without gffcompare tracking columns.

With --reference, transcripts whose id is already a reference id are dropped. Bambu and
IsoQuant emit extended annotations, meaning the whole reference passed through verbatim
alongside their novel models, so the union would otherwise carry one copy of the
annotation per tool per group. Those copies are byte-identical to the reference and are
discarded downstream by class code regardless, so dropping them here only avoids the
memory cost. StringTie has no reference ids, so its known isoforms still ride through and
are still dropped by class code.
"""
import argparse
import re

SEP = "__"


def sub_attr(field, key, newval):
    pat = re.compile(key + r' "([^"]*)"')
    if pat.search(field):
        return pat.sub('%s "%s"' % (key, newval), field, count=1)
    return field


def get_attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", nargs=3, action="append", metavar=("TOOL", "GROUP", "PATH"),
                    required=True, help="repeatable: assembler GTF with its tool + group labels")
    ap.add_argument("--reference", help="reference GTF; transcripts reusing a reference "
                                        "transcript_id are dropped from the union")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    ref_ids = set()
    if args.reference:
        with open(args.reference) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.split("\t")
                if len(c) < 9 or c[2] != "transcript":
                    continue
                rid = get_attr(c[8], "transcript_id")
                if rid is not None:
                    ref_ids.add(rid)
        print("* reference: %d transcript ids -> dropped from the union" % len(ref_ids))

    n_in = 0
    n_skipped = 0
    with open(args.output, "w") as out:
        for tool, group, path in args.gtf:
            pfx = "%s%s%s%s" % (tool, SEP, group, SEP)
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 9 or c[2] not in ("transcript", "exon"):
                        continue
                    tid = get_attr(c[8], "transcript_id")
                    gid = get_attr(c[8], "gene_id")
                    if tid is None:
                        continue
                    if tid in ref_ids:
                        if c[2] == "transcript":
                            n_skipped += 1
                        continue
                    attrs = c[8]
                    attrs = sub_attr(attrs, "transcript_id", pfx + tid)
                    if gid is not None:
                        attrs = sub_attr(attrs, "gene_id", pfx + gid)
                    attrs = attrs.rstrip()
                    if not attrs.endswith(";"):
                        attrs += ";"
                    attrs += ' src_tool "%s"; src_group "%s";' % (tool, group)
                    c[8] = attrs
                    out.write("\t".join(c) + "\n")
                    if c[2] == "transcript":
                        n_in += 1
    print("* union: %d transcripts from %d assembler GTFs (%d reference passthrough dropped) -> %s"
          % (n_in, len(args.gtf), n_skipped, args.output))


if __name__ == "__main__":
    main()
