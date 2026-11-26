# isoseq-pipeline

## Overview

Builds a long-read isoform reference, `isoforms_<run_code>`, from Oxford Nanopore cDNA of two TDP-43 cohorts: knockdown (`tdpkd`) and overexpression (`sun`). Three assemblers discover isoforms per cohort; a pooled stage merges the novel models, filters them with SQANTI3 and adds the survivors to GENCODE v50.

**Status:** `v10` alignment and QC ran on 2026-11-10; assembly and the pooled stage have not run (on hold). `v9` is the published reference. The CHM13 (T2T) configs are paused.

## Quick start

From the thesis repository:

```
# one cohort (sun: cohorts/sun_tdp_overexpression/isoseq-pipeline, sun_tdp_overexpression_nanopore_config.yaml)
cd cohorts/tdpkd/isoseq-pipeline
sh ~/bin/snakejob -n -s ../../../pipelines/isoseq-pipeline/Snakefile -c tdp_nanopore_config.yaml   # dry run
snakejob_HPC -s ../../../pipelines/isoseq-pipeline/Snakefile -c tdp_nanopore_config.yaml           # submit

# pooled stage, after both cohorts
cd pipelines/isoseq-pipeline/cross_cohort
snakejob_HPC -s Snakefile -c cross_cohort_config.yaml

# or all of it (-x: pooled only, -r: resume)
pipelines/isoseq-pipeline/build_reference.sh
```

## Workflow

Per cohort (`Snakefile`):

1. **Read preparation:** pychopper orients reads and trims primers.
2. **Alignment:** minimap2 to GRCh38, guided by GENCODE v50 and STAR short-read introns.
3. **QC:** samtools, Picard, pychopper and NanoStat, collected in MultiQC.
4. **Discovery per condition:** StringTie (`--mix` with short reads, `-G` GENCODE-guided), IsoQuant and Bambu.

Pooled (`cross_cohort/Snakefile`):

5. **Combine** the 12 assembler GTFs (3 tools x 2 cohorts x 2 conditions) and remove known models.
6. **Merge** with isomatch.
7. **Filter** with SQANTI3.
8. **Build the reference:** GENCODE v50 plus the passing novels, renamed `XNOVEL_*`, with their SQANTI3 gene and a `cohorts` attribute.
9. **Check** six TDP-43 positive-control junctions at every step and count models per step.

## Filtering

**Removed before the merge:** transcripts with a GENCODE id; gffcompare `=` or `c` models (same or contained intron chain, removed again after the merge); novel single-exon models (no junction to validate); unstranded models (SQANTI3 discards them).

**Merge:** isomatch 0.8.0 (Zheng et al. 2026, bioRxiv, doi:10.64898/2026.10.09.758114) at defaults, ignoring 5' ends (`--terminal-refine tes --terminal-refine-nc tes`). Models merge when they share the exon count and canonical junctions exactly, with non-canonical junctions within 3 bp and 3' ends within 50 bp. A model from any one tool is kept.

**SQANTI3 6.0.1:** default rules, with STAR junctions from all 18 short-read libraries of both cohorts (`-c`). A novel isoform passes if it shows no intra-priming (< 60% A after its 3' end), no RT-switching, and meets the junction rule:

| Junctions | Short-read support needed |
|---|---|
| All canonical | None |
| Any non-canonical | >= 3 unique reads on its least-covered junction, so on every junction |

## Outputs

| File | Contents |
|---|---|
| `references/novel_isoform_references/<run_code>/isoforms_<run_code>.gtf.gz`, `.fa.gz` | Final reference |
| `.../novel_isoform_discovery_by_tool_cohort_and_condition.tsv` | Tools, cohorts and conditions behind each novel isoform |
| `.../isoforms_<run_code>_sanity_genes.tsv` | Positive controls per step (UNC13A, STMN2, KCNQ2, HYOU1, NUP93, XPNPEP1) |
| `.../isoform_filtering_summary*.tsv` | Model counts per step |
| `.../sqanti/`, `.../consensus/` | SQANTI3 classification; merged models and tool overlap |
| `cohorts/<cohort>/isoseq-pipeline/results/<run_code>/` | CRAM alignments, MultiQC, assembler GTFs |

## Technical details

| Cohort | Long reads | Short reads |
|---|---|---|
| tdpkd | 11 | 10 (TDP-5, CTRL-4 excluded: `exclude_samples`) |
| sun | 8 aligned, 7 for discovery (GFP_1: `exclude_from_discovery`) | 8 |

- **Settings:** pychopper kit `PCS109` (tdpkd) or `PCB114` (sun); minimap2 `splice`, k = 15; short-read introns need >= 2 unique reads; IsoQuant `sensitive_ont`; Bambu NDR 0.5. Conditions are sample-name prefixes (`TDP_1` belongs to `TDP`).
- **Software:** pychopper 2.7.10, minimap2 2.31, samtools 1.10, NanoStat 1.6.0, Picard, MultiQC 1.25.1, StringTie 3.0.3, IsoQuant 3.13.0, Bambu 3.12.1, gffcompare 0.12.10, isomatch 0.8.0, SQANTI3 6.0.1, gffread.
- **Environments:** conda `isoseq-pipeline`, `pychopper_env`, `multiqc`, `isoquant313`, `bambu_env`, `sqanti3_v6`; minimap2, StringTie and gffcompare paths are set in the workflow files; isomatch is the release binary at `/hpc/users/friedb09/tools/isomatch-0.8.0/isomatch`.
