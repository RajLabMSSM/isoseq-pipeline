#!/usr/bin/env python
"""Combine the tool x cohort x group assembler GTFs into one union GTF, without the transcripts that carry a GENCODE id.
Each kept transcript gets the id <tool>__<cohort>__<group>__<own id> and the tags src_tool, src_cohort and src_group."""
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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", nargs=4, action="append", metavar=("TOOL", "COHORT", "GROUP", "PATH"), required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    # every reference transcript id
    ref_ids = set()
    with open(args.reference) as fh:
        for line in fh:
            cols = line.split("\t")
            if not line.startswith("#") and len(cols) >= GTF_N_COLUMNS and cols[GTF_FEATURE] == "transcript":
                ref_ids.add(get_attr(cols[GTF_ATTRIBUTES], "transcript_id"))
    ref_ids.discard(None)
    print("* reference: %d transcript ids" % len(ref_ids))

    # removed transcripts are counted once each, on their transcript line
    n_in = 0
    n_ref_id = 0
    n_reference_id = 0
    skipped = set()
    with open(args.output, "w") as out:
        for tool, cohort, group, path in args.gtf:
            prefix = "%s__%s__%s__" % (tool, cohort, group)
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
                    # bambu and isoquant pass the reference through under its own GENCODE transcript_id
                    if tid in ref_ids:
                        if cols[GTF_FEATURE] == "transcript":
                            n_ref_id += 1
                        continue
                    # stringtie -G keeps its STRG id and names the matched GENCODE transcript in reference_id, on the
                    # transcript line only, so its exon lines are dropped through the remembered (file, id) pair
                    if cols[GTF_FEATURE] == "transcript" and get_attr(cols[GTF_ATTRIBUTES], "reference_id") in ref_ids:
                        n_reference_id += 1
                        skipped.add((path, tid))
                        continue
                    if (path, tid) in skipped:
                        continue
                    attrs = sub_attr(cols[GTF_ATTRIBUTES], "transcript_id", prefix + tid)
                    if gid is not None:
                        attrs = sub_attr(attrs, "gene_id", prefix + gid)
                    attrs = attrs.rstrip()
                    if not attrs.endswith(";"):
                        attrs += ";"
                    cols[GTF_ATTRIBUTES] = attrs + ' src_tool "%s"; src_cohort "%s"; src_group "%s";' % (tool, cohort, group)
                    out.write("\t".join(cols) + "\n")
                    if cols[GTF_FEATURE] == "transcript":
                        n_in += 1
    print("* union: %d transcripts from %d assembler GTFs; removed %d with a GENCODE transcript_id and %d with a "
          "GENCODE reference_id -> %s" % (n_in, len(args.gtf), n_ref_id, n_reference_id, args.output))


if __name__ == "__main__":
    main()
