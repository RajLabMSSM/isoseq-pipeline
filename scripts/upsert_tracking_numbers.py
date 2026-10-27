#!/usr/bin/env python3
"""Upsert one isoseq build into the project tracking table: per-cohort funnel rows plus the
cross-cohort rows under every cohort, replacing that cohort's older isoseq rows.

usage: upsert_tracking_numbers.py <tracking_numbers.tsv> <run_code> <union_tracking.tsv> <cohort>=<tracking.tsv> ..."""
import csv
import fcntl
import os
import sys

HEADER = ["cohort", "metric", "value"]
OWN_PREFIXES = ("raw_", "novel_", "consensus_", "collapsed_", "overlap_", "sqanti_", "lr_specific",
                "sr_supported", "final_", "reference_transcripts", "pct_", "crosscohort_",
                "cross_cohort_", "sr_reference_", "gencode_", "by_combo", "single_cohort",
                "shared_by", "isoseq_run_code")


def read_metrics(path):
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return [(r[0], r[1]) for r in rows[1:] if len(r) >= 2 and r[0] != "cohort"]


def main(tracking, run_code, union, cohort_args):
    cohorts = dict(a.split("=", 1) for a in cohort_args)
    union_rows = read_metrics(union)
    new = []
    for cohort, path in cohorts.items():
        new += [(cohort, m, v) for m, v in read_metrics(path)]
        new += [(cohort, m, v) for m, v in union_rows]
        new.append((cohort, "isoseq_run_code", run_code))

    lock = open(tracking + ".lock", "w")
    fcntl.lockf(lock, fcntl.LOCK_EX)
    rows = []
    if os.path.exists(tracking):
        with open(tracking) as fh:
            reader = csv.reader(fh, delimiter="\t")
            if next(reader) != HEADER:
                raise SystemExit("unexpected header in " + tracking)
            rows = [r for r in reader if len(r) == 3]
    kept = [r for r in rows if not (r[0] in cohorts and r[1].startswith(OWN_PREFIXES))]
    dropped = len(rows) - len(kept)
    out = []
    for cohort in cohorts:
        out += [r for r in kept if r[0] == cohort]
        out += [list(r) for r in new if r[0] == cohort]
    out += [r for r in kept if r[0] not in cohorts]
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
