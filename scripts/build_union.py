#!/usr/bin/env python
"""Concatenate the tool x group assembler gtfs into one union, ids prefixed <tool>__<group>__ and
tagged src_tool / src_group. Transcripts reusing a reference id (the passthrough of bambu and
isoquant extended annotations) are dropped."""
import argparse
import re

# GOAL = take the 12 input GTF files, remove transcripts whose IDs already exist in GENCODE, give every remaining transcript a unique identifier, 
# record where it came from, and concatenate everything into one GTF.

# gtf has 9 tab-separated columns; column 3 (index 2) is the feature type, column 9 (index 8) the attributes
GTF_N_COLUMNS = 9
GTF_FEATURE = 2
GTF_ATTRIBUTES = 8


def get_attr(field, key):
    # value of one gtf attribute (key "value";), or None when the key is missing
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def sub_attr(field, key, newval):
    # replace the first value of one gtf attribute
    return re.sub(key + r' "([^"]*)"', '%s "%s"' % (key, newval), field, count=1)


def main():
    # read the options: one --gtf TOOL GROUP PATH triple per assembler output
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", nargs=3, action="append", metavar=("TOOL", "GROUP", "PATH"), required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    # collect every reference transcript id, so the passthrough copies can be dropped below
    ref_ids = set()
    with open(args.reference) as fh:
        for line in fh:
            cols = line.split("\t")
            if not line.startswith("#") and len(cols) >= GTF_N_COLUMNS and cols[GTF_FEATURE] == "transcript":
                ref_ids.add(get_attr(cols[GTF_ATTRIBUTES], "transcript_id"))
    ref_ids.discard(None)
    print("* reference: %d transcript ids -> dropped from the union" % len(ref_ids))

    # write the transcript and exon lines of every assembler gtf, ids made unique by the tool + group prefix
    n_in = 0
    n_skipped = 0
    with open(args.output, "w") as out:
        for tool, group, path in args.gtf:
            prefix = "%s__%s__" % (tool, group)
            with open(path) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if line.startswith("#") or len(cols) < GTF_N_COLUMNS:
                        continue
                    if cols[GTF_FEATURE] not in ("transcript", "exon"):
                        continue
                    tid = get_attr(cols[GTF_ATTRIBUTES], "transcript_id")
                    gid = get_attr(cols[GTF_ATTRIBUTES], "gene_id")
                    if tid is None:
                        continue
                    # a reference id = an annotated transcript passed through by the tool: drop it, count it once
                    if tid in ref_ids:
                        if cols[GTF_FEATURE] == "transcript":
                            n_skipped += 1
                        continue
                    # prefix the ids and add the source tags
                    attrs = sub_attr(cols[GTF_ATTRIBUTES], "transcript_id", prefix + tid)
                    if gid is not None:
                        attrs = sub_attr(attrs, "gene_id", prefix + gid)
                    attrs = attrs.rstrip()
                    if not attrs.endswith(";"):
                        attrs += ";"
                    cols[GTF_ATTRIBUTES] = attrs + ' src_tool "%s"; src_group "%s";' % (tool, group)
                    out.write("\t".join(cols) + "\n")
                    if cols[GTF_FEATURE] == "transcript":
                        n_in += 1
    print("* union: %d transcripts from %d assembler GTFs (%d reference passthrough dropped) -> %s"
          % (n_in, len(args.gtf), n_skipped, args.output))


if __name__ == "__main__":
    main()
