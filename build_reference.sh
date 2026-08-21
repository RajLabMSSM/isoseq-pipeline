#!/usr/bin/env bash
# build_reference.sh — end-to-end cross-cohort reference build in one command.
#
#   ./build_reference.sh          full build: per-cohort assembly (all cohorts) then cross_cohort
#   ./build_reference.sh -r       same, but resume (--rerun-triggers mtime; finished stages not redone)
#   ./build_reference.sh -x       cross_cohort ONLY (per-cohort already complete)
#
# Launches each cohort's per-cohort isoseq assembly, captures its head-job id, then submits
# cross_cohort with an LSF dependency `done(id1) && done(id2) && ...` so the merge fires
# automatically once every cohort finishes successfully. cross_cohort emits the GTF, the
# transcriptome FASTA (rule cross_cohort_fasta), and provenance in that single run.
# Edit COHORTS to match the cohorts listed in cross_cohort/cross_cohort_config.yaml.
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
    h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
  esac
done

launch_cross() {
  local dep="$1"
  cd "$XC"
  mkdir -p logs
  if [ -z "$dep" ]; then
    snakejob_HPC -s Snakefile -c cross_cohort_config.yaml $RESUME
  else
    bsub -P acc_als-omics -q premium -W 1:00 -n 1 -R "rusage[mem=4000]" \
      -w "$dep" -J cross_cohort_trigger \
      -o logs/cross_cohort_trigger.stdout -e logs/cross_cohort_trigger.stderr \
      "cd $XC && snakejob_HPC -s Snakefile -c cross_cohort_config.yaml $RESUME"
    echo "cross_cohort will auto-launch once dependencies clear: $dep"
  fi
}

if [ -n "$CROSS_ONLY" ]; then
  echo "== cross_cohort only =="
  launch_cross ""
  exit 0
fi

dep=""
for entry in "${COHORTS[@]}"; do
  dir="${entry%%:*}"; cfg="${entry##*:}"
  cd "$dir"
  out=$(snakejob_HPC -s Snakefile -c "$cfg" $RESUME)
  echo "$out"
  id=$(echo "$out" | grep -oP 'Job <\K[0-9]+' | tail -1)
  [ -z "$id" ] && { echo "ERROR: no job id captured for $dir" >&2; exit 1; }
  echo "-> per-cohort $(basename "$(dirname "$dir")") head job $id"
  dep="${dep:+$dep && }done($id)"
done

launch_cross "$dep"
