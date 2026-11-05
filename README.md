# isoseq-pipeline

Builds the long-read isoform atlas (`isoforms_v9`): per-cohort isoform discovery from ONT cDNA reads,
then one cross-cohort union with GENCODE.

## Run

Each cohort folder (`cohorts/<cohort>/isoseq-pipeline/`) holds its config, metadata and logs; the code
lives here. `run_code` in the config names the output folder.

```
cd cohorts/tdpkd/isoseq-pipeline
sh ~/bin/snakejob -n -s ../../../pipelines/isoseq-pipeline/Snakefile -c tdp_nanopore_config.yaml   # dry run
snakejob_HPC -s ../../../pipelines/isoseq-pipeline/Snakefile -c tdp_nanopore_config.yaml           # submit
```

`snakejob_HPC -n` submits; dry runs go through `snakejob -n`. Then the cross-cohort union:

```
cd pipelines/isoseq-pipeline/cross_cohort
snakejob_HPC -s Snakefile -c cross_cohort_config.yaml
```

`build_reference.sh` launches both cohorts and chains the union behind them (`-x`: union only).

Metadata columns: `sample`, `condition`, `long_read_fastq`, `short_read_bam_path`. Groups are
sample-name prefixes (`TDP_1` -> `TDP`).

## Steps (v9)

1. pychopper 2.7.10 (`-m edlib`; kit PCS109 tdpkd, PCS111 sun), polyA kept
2. minimap2 2.31 `-ax splice` (all cohorts, default k15 from the index), `--secondary=no -G 500k`,
   `--junc-bed` = GENCODE v50 junctions + the cohort's STAR junctions (>= 2 unique reads, canonical).
   tdpkd drops two short-read libraries via `exclude_samples`, here and in the SQANTI junction coverage:
   TDP-5 (aberrant differentiation) and CTRL-4 (low-depth outlier). The published v9 build (2026-09-21)
   dropped them from the SQANTI coverage only; its tdpkd junc-bed used all 12 STAR libraries (508,354
   short-read junctions, vs 500,324 with the two dropped)
3. Per group, on the pooled long-read BAM:
   - StringTie 3.0.3 `--mix -G` (short-read BAM first; no `-t`)
   - IsoQuant 3.13.0 `default_ont`, `--illumina_bam`, `--clean_start`, `--read_group file_name`
   - Bambu 3.12.1, NDR 0.7
4. Union of tool x group models; models in GENCODE (`=`, `c`), single-exon and unstranded models
   dropped; `tama_merge -m 6 -z 100 -e common_ends -d merge_dup`, one `no_cap` source per tool and
   group (no tool vote); gffcompare 0.12.10 for class codes
5. SQANTI3 6.0.1 QC with the cohort's `SJ.out.tab`, CAGE peaks and polyA motifs; `sqanti3_filter.py
   rules -e` with `misc/sqanti_rules.json`
6. GENCODE v50 plus the passed novels (`NOVEL_*`, gene from SQANTI's `associated_gene`), GTF + FASTA
7. Cross-cohort: the cohorts' novels merged by `tama_merge` at the same tolerances, single-exon
   dropped, gene ids set by `cross_cohort_novel_genes.py` (`XNOVEL_*`, `XNOVELG_*` loci)

QC: NanoStat, samtools, RNA-SeQC, Picard into MultiQC 1.25.1, plus `<data_code>_longread_qc_summary.tsv`
and `<data_code>_read_length_hist.tsv.gz`.

Conda environments: `isoseq-pipeline` (`isoseq-pipeline.yml`), `pychopper_env`, `bambu_env`,
`isoquant313`, `sqanti3_v6`, `tama_env` (Python 2), `multiqc`.

## Outputs

`results/<run_code>/consensus/`: `<data_code>.isoforms.{gtf,fa}.gz`, `_novel_provenance.tsv`,
`_tool_overlap.tsv`, `_tracking_numbers.tsv` (the funnel), `_sanity_genes.tsv{,.events.tsv}` (the STMN2
and UNC13A cryptic exons traced stage by stage). Cross-cohort: `references/novel_isoform_references/v9/`.
