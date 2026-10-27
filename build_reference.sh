#!/usr/bin/env bash
# Launch every cohort's isoseq run, then the cross-cohort union once all of them finish (LSF dependency).
#   ./build_reference.sh       full build
#   ./build_reference.sh -r    resume (--rerun-triggers mtime)
#   ./build_reference.sh -x    cross-cohort union only
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
  case "$flag" in
    r) RESUME="-r" ;;
    x) CROSS_ONLY="1" ;;
    *) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
  esac
done
shift $((OPTIND - 1))
[ $# -eq 0 ] || { echo "usage: build_reference.sh [-r] [-x]" >&2; exit 1; }

launch_cross() {
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
  dir="${entry%%:*}"; cfg="${entry##*:}"
  cd "$dir"
  # the real path, not the cohort's symlink: workflow.basedir must be this repo (misc/sqanti_rules.json)
  out=$(snakejob_HPC -s "$REPO/pipelines/isoseq-pipeline/Snakefile" -c "$cfg" $RESUME)
  echo "$out"
  id=$(echo "$out" | grep -oP 'Job <\K[0-9]+' | tail -1)
  [ -n "$id" ] || { echo "ERROR: no job id captured for $dir" >&2; exit 1; }
  dep="${dep:+$dep && }done($id)"
done

launch_cross "$dep"
