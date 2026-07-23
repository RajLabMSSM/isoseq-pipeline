#!/usr/bin/env python
"""
nanostat_mqc.py -- render the collated NanoStat table into MultiQC custom content.

MultiQC has no NanoStat module, so long-read stats are invisible in multiqc_report.html
unless injected as "custom content" (*_mqc.tsv is picked up generically, no module needed).

DESIGN -- readable view, complete record:
  * TWO tables, one per stage (raw / aligned), rather than one wide table with every metric
    prefixed raw_/aligned_. Splitting by stage halves the width and keeps NanoStat's OWN
    metric names verbatim, so the report reads like NanoStat's output rather than a
    reshaped version of it.
  * Only the SUMMARY metrics are shown. NanoStat's per-read top-5 extremes (longest_read,
    highest_Q_read) and the count/megabase variants of the Q bins are NOT displayed -- they
    are single-read trivia / redundant with the percentage, and 95 columns is unreadable.
    They are NOT lost: the collated TSV holds every field, and the per-sample
    <sample>.nanostat.{raw,aligned}.tsv files are untouched on disk.
  * The aligned table gains pct_reads_aligned (aligned/raw) -- the stage-to-stage drop is
    the metric that exposed the tdpkd yield deficit, and it exists in neither stage alone.

Columns are looked up by name, so a missing metric degrades to NA rather than shifting the
table.
"""
import argparse


SHOW = {
    "raw": [
        "number_of_reads",
        "number_of_bases",
        "mean_read_length",
        "median_read_length",
        "read_length_stdev",
        "n50",
        "mean_qual",
        "median_qual",
        "pct_reads_above_q10",
        "pct_reads_above_q15",
    ],
    "aligned": [
        "number_of_reads",
        "number_of_bases",
        "number_of_bases_aligned",
        "fraction_bases_aligned",
        "average_identity",
        "median_identity",
        "mean_read_length",
        "median_read_length",
        "read_length_stdev",
        "n50",
        "mean_qual",
        "median_qual",
        "pct_reads_above_q10",
        "pct_reads_above_q15",
    ],
}

DESC = {
    "raw": "NanoStat on the raw basecalled FASTQ, before pychopper. Note mean vs median read "
           "length differ substantially (the distribution is right-skewed) -- state which one "
           "you are quoting.",
    "aligned": "NanoStat on the aligned BAM. pct_reads_aligned = aligned / raw reads, so it "
               "captures losses at BOTH pychopper and minimap2. average_identity is the "
               "cleanest read-accuracy / chemistry metric.",
}


def header(stage):
    return (
        "# id: 'nanostat_%s'\n"
        "# section_name: 'NanoStat -- %s reads'\n"
        "# description: '%s'\n"
        "# plot_type: 'table'\n"
        "# pconfig:\n"
        "#     id: 'nanostat_%s_table'\n"
        "#     namespace: 'NanoStat'\n"
        "#     scale: false\n"
    ) % (stage, stage, DESC[stage], stage)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="collated nanostat TSV")
    ap.add_argument("-o", "--output", required=True,
                    help="output prefix; writes <prefix>_raw_mqc.tsv and <prefix>_aligned_mqc.tsv")
    args = ap.parse_args()

    with open(args.input) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(head)}
        rows = {}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < len(head):
                continue
            rows[(c[idx["sample"]], c[idx["stage"]])] = c

    samples = sorted({s for s, _ in rows})

    def get(sample, stage, field):
        r = rows.get((sample, stage))
        if not r or field not in idx:
            return "NA"
        return r[idx[field]]

    for stage in ("raw", "aligned"):
        cols = SHOW[stage]
        path = "%s_%s_mqc.tsv" % (args.output, stage)
        with open(path, "w") as out:
            out.write(header(stage))
            extra = "\tpct_reads_aligned" if stage == "aligned" else ""
            out.write("Sample\t" + "\t".join(cols) + extra + "\n")
            for s in samples:
                vals = [get(s, stage, c) for c in cols]
                if stage == "aligned":
                    try:
                        raw = float(get(s, "raw", "number_of_reads"))
                        aln = float(get(s, "aligned", "number_of_reads"))
                        vals.append("%.1f" % (100.0 * aln / raw) if raw else "NA")
                    except ValueError:
                        vals.append("NA")
                out.write(s + "\t" + "\t".join(vals) + "\n")
        print("* %s -> %s (%d samples, %d cols)"
              % (stage, path, len(samples), len(cols) + (1 if stage == "aligned" else 0)))


if __name__ == "__main__":
    main()
