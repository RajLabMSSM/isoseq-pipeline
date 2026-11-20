#!/usr/bin/env python
"""Combine the assembler GTFs of every cohort into one union GTF and remove the transcripts that carry a GENCODE id.

Each input is one tool x cohort x group GTF. Every kept transcript gets the id <tool>__<cohort>__<group>__<own id>
and the tags src_tool, src_cohort and src_group, so each model can be traced back to the tool, cohort and
condition that built it.

A transcript carries a GENCODE id in one of two ways. Bambu and IsoQuant write the whole reference into their
extended annotations under the reference transcript_id. StringTie -G keeps its own STRG id and names the matching
reference transcript in a reference_id attribute. Both kinds are removed here."""
import argparse
import re

# A GTF line has 9 tab-separated columns. Column 3 (index 2) is the feature type and column 9 (index 8) holds the
# attributes.
GTF_N_COLUMNS = 9
GTF_FEATURE = 2
GTF_ATTRIBUTES = 8


def get_attr(field, key):
    # Return the value of one GTF attribute, written as key "value";, or None when the key is missing.
    m = re.search(key + r' "([^"]*)"', field)
    return m.group(1) if m else None


def sub_attr(field, key, newval):
    # Replace the first value of one GTF attribute.
    return re.sub(key + r' "([^"]*)"', '%s "%s"' % (key, newval), field, count=1)


def main():
    # Read the options. Each assembler output is given as --gtf TOOL COHORT GROUP PATH.
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", nargs=4, action="append", metavar=("TOOL", "COHORT", "GROUP", "PATH"), required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    # Collect every reference transcript id, so that transcripts carrying one can be removed below.
    ref_ids = set()
    with open(args.reference) as fh:
        for line in fh:
            cols = line.split("\t")
            if not line.startswith("#") and len(cols) >= GTF_N_COLUMNS and cols[GTF_FEATURE] == "transcript":
                ref_ids.add(get_attr(cols[GTF_ATTRIBUTES], "transcript_id"))
    ref_ids.discard(None)
    print("* reference: %d transcript ids" % len(ref_ids))

    # Write the transcript and exon lines of every assembler GTF. The tool, cohort and group prefix keeps ids unique.
    # The removed transcripts are counted once each, on their transcript line, separately for the two id kinds.
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
                    # A transcript whose own id is a GENCODE id is the reference copy that Bambu and IsoQuant pass
                    # through, so it is removed.
                    if tid in ref_ids:
                        if cols[GTF_FEATURE] == "transcript":
                            n_ref_id += 1
                        continue
                    # StringTie names the GENCODE transcript it matched in reference_id, and only on the transcript
                    # line. On the v9 tdpkd CTRL output all 224,734 such models were gffcompare class =, so removing
                    # them here changes no result. Their ids are remembered per file so that their exon lines,
                    # which carry no reference_id, are removed too.
                    if cols[GTF_FEATURE] == "transcript" and get_attr(cols[GTF_ATTRIBUTES], "reference_id") in ref_ids:
                        n_reference_id += 1
                        skipped.add((path, tid))
                        continue
                    if (path, tid) in skipped:
                        continue
                    # Prefix the ids and add the source tags.
                    attrs = sub_attr(cols[GTF_ATTRIBUTES], "transcript_id", prefix + tid)
                    if gid is not None:
                        attrs = sub_attr(attrs, "gene_id", prefix + gid)
                    attrs = attrs.rstrip()
                    if not attrs.endswith(";"):
                        attrs += ";"
                    cols[GTF_ATTRIBUTES] = attrs + ' src_tool "%s"; src_cohort "%s"; src_group "%s";' % (tool, cohort,
                                                                                                         group)
                    out.write("\t".join(cols) + "\n")
                    if cols[GTF_FEATURE] == "transcript":
                        n_in += 1
    print("* union: %d transcripts from %d assembler GTFs; removed %d with a GENCODE transcript_id and %d with a "
          "GENCODE reference_id -> %s" % (n_in, len(args.gtf), n_ref_id, n_reference_id, args.output))


if __name__ == "__main__":
    main()
