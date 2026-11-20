#!/usr/bin/env bash
# Launch every cohort's isoseq run (discovery), then start the pooled stage once all of them finish.
# An LSF job dependency makes the pooled stage wait for the cohort runs.
#   ./build_reference.sh       full build
#   ./build_reference.sh -r    resume (--rerun-triggers mtime)
#   ./build_reference.sh -x    pooled stage only
set -euo pipefail

REPO=/sc/arion/projects/als-omics/brooke-phd-thesis
XC=$REPO/pipelines/isoseq-pipeline/cross_cohort
COHORTS=(
  "$REPO/cohorts/tdpkd/isoseq-pipeline:tdp_nanopore_config.yaml"
  "$REPO/cohorts/sun_tdp_overexpression/isoseq-pipeline:sun_tdp_overexpression_nanopore_config.yaml"
)

RESUME=""
CROSS_ONLY=""
while getopts 'rxh' flag; do
  # -r resumes the runs and -x runs only the pooled stage.
  # Any other flag prints the usage lines at the top of this file.
  case "$flag" in
    r) RESUME="-r" ;;
    x) CROSS_ONLY="1" ;;
    *) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
  esac
done
shift $((OPTIND - 1))
[ $# -eq 0 ] || { echo "usage: build_reference.sh [-r] [-x]" >&2; exit 1; }

launch_cross() {
  # Start the pooled stage straight away when no dependency is given.
  # Otherwise submit a small LSF job that waits (-w) for the cohort runs to finish and then starts it.
  cd "$XC"
  mkdir -p logs
  if [ -z "$1" ]; then
    snakejob_HPC -s Snakefile -c cross_cohort_config.yaml $RESUME
  else
    bsub -P acc_als-omics -q premium -W 1:00 -n 1 -R "rusage[mem=4000]" \
      -w "$1" -J cross_cohort_trigger \
      -o logs/cross_cohort_trigger.stdout -e logs/cross_cohort_trigger.stderr \
      "cd $XC && snakejob_HPC -s Snakefile -c cross_cohort_config.yaml $RESUME"
    echo "cross_cohort will launch once these finish: $1"
  fi
}

if [ -n "$CROSS_ONLY" ]; then
  launch_cross ""
  exit 0
fi

dep=""
for entry in "${COHORTS[@]}"; do
  # Submit each cohort run and add its LSF job id to one dependency expression.
  # The expression reads done(<id1>) && done(<id2>), so the pooled-stage trigger waits for all of them.
  # Each entry is written as <cohort folder>:<config file>, so it is split at the colon.
  dir="${entry%%:*}"
  cfg="${entry##*:}"
  cd "$dir"
  # Use the real path to the Snakefile, not the cohort's symlink.
  # workflow.basedir must point to this repo, so that the workflows and scripts are found.
  out=$(snakejob_HPC -s "$REPO/pipelines/isoseq-pipeline/Snakefile" -c "$cfg" $RESUME)
  echo "$out"
  # Take the job id from the last "Job <id> is submitted" line that LSF printed.
  id=$(echo "$out" | grep -oP 'Job <\K[0-9]+' | tail -1)
  [ -n "$id" ] || { echo "ERROR: no job id captured for $dir" >&2; exit 1; }
  if [ -z "$dep" ]; then
    dep="done($id)"
  else
    dep="$dep && done($id)"
  fi
done

launch_cross "$dep"
