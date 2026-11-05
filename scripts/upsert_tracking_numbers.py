#!/usr/bin/env python3
"""Upsert one isoseq build into the project tracking table: per-cohort funnel rows plus the
cross-cohort rows under every cohort, replacing that cohort's older isoseq rows.

usage: upsert_tracking_numbers.py <tracking_numbers.tsv> <run_code> <union_tracking.tsv> <cohort>=<tracking.tsv> ..."""
import csv
import fcntl
import os
import sys

HEADER = ["cohort", "metric", "value"]
# metric name prefixes this pipeline writes: rows of an updated cohort starting with one of these are
# replaced; rows written by other pipelines are left as they are
OWN_PREFIXES = ("raw_", "novel_", "consensus_", "collapsed_", "overlap_", "sqanti_", "lr_specific",
                "sr_supported", "final_", "reference_transcripts", "pct_", "crosscohort_",
                "cross_cohort_", "sr_reference_", "gencode_", "by_combo", "single_cohort",
                "shared_by", "isoseq_run_code")


def read_metrics(path):
    # (metric, value) pairs of a metric<TAB>value table; skip the header, short rows and the cohort row
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return [(r[0], r[1]) for r in rows[1:] if len(r) >= 2 and r[0] != "cohort"]


def main(tracking, run_code, union, cohort_args):
    # new rows: each cohort's funnel, then the cross-cohort rows, then the run code
    cohorts = dict(a.split("=", 1) for a in cohort_args)
    union_rows = read_metrics(union)
    new = []
    for cohort, path in cohorts.items():
        new += [(cohort, m, v) for m, v in read_metrics(path)]
        new += [(cohort, m, v) for m, v in union_rows]
        new.append((cohort, "isoseq_run_code", run_code))

    # lock the shared table so two cohort runs cannot write it at the same time
    with open(tracking + ".lock", "w") as lock:
        fcntl.lockf(lock, fcntl.LOCK_EX)
        # read the existing table
        rows = []
        if os.path.exists(tracking):
            with open(tracking) as fh:
                reader = csv.reader(fh, delimiter="\t")
                if next(reader) != HEADER:
                    raise SystemExit("unexpected header in " + tracking)
                rows = [r for r in reader if len(r) == 3]
        # drop this pipeline's old rows for the updated cohorts
        kept = [r for r in rows if not (r[0] in cohorts and r[1].startswith(OWN_PREFIXES))]
        dropped = len(rows) - len(kept)
        # updated cohorts first (their other rows, then the new rows), then every other cohort untouched
        out = []
        for cohort in cohorts:
            out += [r for r in kept if r[0] == cohort]
            out += [list(r) for r in new if r[0] == cohort]
        out += [r for r in kept if r[0] not in cohorts]
        # write to a temporary file and swap it in, so a crash never leaves a half-written table
        tmp = tracking + ".tmp"
        with open(tmp, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(HEADER)
            w.writerows(out)
        os.replace(tmp, tracking)
        fcntl.lockf(lock, fcntl.LOCK_UN)
    sys.stderr.write(f"{tracking}: {dropped} isoseq rows replaced by {len(new)}, {len(out)} total\n")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
