#!/usr/bin/env python
"""Concatenate the tool x group assembler gtfs into one union, ids prefixed <tool>__<group>__ and
tagged src_tool / src_group. Transcripts reusing a reference id (the passthrough of bambu and
isoquant extended annotations) are dropped."""
import argparse
import re


def get_attr(field, key):
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def sub_attr(field, key, newval):
    return re.sub(key + r' "([^"]*)"', '%s "%s"' % (key, newval), field, count=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", nargs=3, action="append", metavar=("TOOL", "GROUP", "PATH"), required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    ref_ids = set()
    with open(args.reference) as fh:
        for line in fh:
            c = line.split("\t")
            if not line.startswith("#") and len(c) >= 9 and c[2] == "transcript":
                ref_ids.add(get_attr(c[8], "transcript_id"))
    ref_ids.discard(None)
    print("* reference: %d transcript ids -> dropped from the union" % len(ref_ids))

    n_in = n_skipped = 0
    with open(args.output, "w") as out:
        for tool, group, path in args.gtf:
            pfx = "%s__%s__" % (tool, group)
            with open(path) as fh:
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if line.startswith("#") or len(c) < 9 or c[2] not in ("transcript", "exon"):
                        continue
                    tid = get_attr(c[8], "transcript_id")
                    gid = get_attr(c[8], "gene_id")
                    if tid is None:
                        continue
                    if tid in ref_ids:
                        n_skipped += c[2] == "transcript"
                        continue
                    attrs = sub_attr(c[8], "transcript_id", pfx + tid)
                    if gid is not None:
                        attrs = sub_attr(attrs, "gene_id", pfx + gid)
                    attrs = attrs.rstrip()
                    if not attrs.endswith(";"):
                        attrs += ";"
                    c[8] = attrs + ' src_tool "%s"; src_group "%s";' % (tool, group)
                    out.write("\t".join(c) + "\n")
                    n_in += c[2] == "transcript"
    print("* union: %d transcripts from %d assembler GTFs (%d reference passthrough dropped) -> %s"
          % (n_in, len(args.gtf), n_skipped, args.output))


if __name__ == "__main__":
    main()
