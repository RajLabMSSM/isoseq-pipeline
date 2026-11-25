# isoseq-pipeline

Builds the long-read isoform atlas (`isoforms_<run_code>`): per-cohort isoform discovery from ONT cDNA
reads, then one pooled stage over both cohorts that removes GENCODE transcripts, merges, classifies against
GENCODE v50, filters with SQANTI3 and appends GENCODE. The code and configs are set for `v10`, whose alignment and QC
ran on 2026-11-10 and whose discovery and pooled stage have not run yet; the last completed build is `v9`.

## Run

Each cohort folder (`cohorts/<cohort>/isoseq-pipeline/`) holds its config, metadata and logs; the code
lives here. `run_code` in the config names the output folder.

```
cd cohorts/tdpkd/isoseq-pipeline
sh ~/bin/snakejob -n -s ../../../pipelines/isoseq-pipeline/Snakefile -c tdp_nanopore_config.yaml   # dry run
snakejob_HPC -s ../../../pipelines/isoseq-pipeline/Snakefile -c tdp_nanopore_config.yaml           # submit
```

`snakejob_HPC -n` submits; dry runs go through `snakejob -n`. Then the pooled stage, once both cohort runs
have finished:

```
cd pipelines/isoseq-pipeline/cross_cohort
snakejob_HPC -s Snakefile -c cross_cohort_config.yaml
```

`cross_cohort_config.yaml` names the two cohort configs; the pooled stage reads the genome, annotation,
SQANTI evidence, assembler output paths and junction files from them, and stops when a cohort config has
another `run_code`. `cross_cohort_config_chm13.yaml` is the T2T copy (no tracking-table rows).
`build_reference.sh` launches both cohorts and chains the pooled stage behind them (`-x`: pooled stage only).

Metadata columns: `sample`, `condition`, `long_read_fastq`, `short_read_bam_path`. Groups are
sample-name prefixes (`TDP_1` -> `TDP`).

## Steps (v10)

1. pychopper 2.7.10 (`-m edlib`; kit PCS109 tdpkd, PCB114 sun), polyA kept
2. minimap2 2.31 `-ax splice` (all cohorts, default k15 from the index), `--secondary=no -G 500k`,
   `--junc-bed` = GENCODE v50 junctions + the cohort's STAR junctions (>= 2 unique reads, canonical).
   tdpkd drops two short-read libraries via `exclude_samples`, here, in the StringTie short-read pool
   and in the SQANTI junction coverage: TDP-5 (aberrant differentiation) and CTRL-4 (low-depth outlier).
   CTRL4 long reads stay in, so tdpkd has 11 long-read and 10 short-read libraries. TDP2 long reads
   combine both of its sequencing runs. sun leaves GFP_1 long reads out of the group pools via
   `exclude_from_discovery` (aligned and QC'd, short reads kept), so sun pools 7 long-read libraries
3. Long-read alignments are merged per group with `samtools merge`, the short-read STAR alignments
   likewise. Per group, on the pooled long-read BAM:
   - StringTie 3.0.3 `--mix -G` (pooled short-read BAM first; no `-t`)
   - IsoQuant 3.13.0 `sensitive_ont`, `--complete_genedb --clean_start --large_output none`, long reads only
   - Bambu 3.12.1, NDR 0.5, default discovery gates
Steps 4 onwards run once, in `cross_cohort/Snakefile`, on both cohorts together.

4. The 12 assembler GTFs (3 tools x 2 groups x 2 cohorts) combined; transcripts carrying a GENCODE id
   removed (Bambu and IsoQuant reference passthrough by `transcript_id`, StringTie `-G` matches by
   `reference_id`); gffcompare 0.12.10 against GENCODE v50, models of class `=` or `c`, single-exon and
   unstranded models dropped; isomatch 0.8.0 `merge --terminal-refine tes --terminal-refine-nc tes`, otherwise
   at its defaults (Zheng et al. 2026, bioRxiv, doi:10.64898/2026.10.09.758114), one input GTF per tool,
   cohort and group (no tool vote); the merged models classified against GENCODE v50 with gffcompare again,
   `=` and `c` dropped again
5. SQANTI3 6.0.1 QC with the `SJ.out.tab` files of both cohorts (10 tdpkd + 8 sun libraries), CAGE peaks
   and polyA motifs; `sqanti3_filter.py rules -e` with SQANTI3's default rules (intra-priming, RT switching, and canonical
   junctions or `min_cov >= 3` short reads summed over the libraries of both cohorts)
6. SQANTI-browser 1.1.2 (Paniagua et al. 2026, bioRxiv, doi:10.64898/2026.05.25.727625) turns
   `_corrected.gtf` and `_RulesFilter_classification.txt` into a UCSC track hub (hg38, chromosome sizes
   from the genome `.fai`, the SQANTI CAGE peaks as a track, `--hub-name isoforms_<run_code>`); `hubCheck`
   validates it. Each model's details page carries its SQANTI columns, `filter_result` included
   (`Isoform` kept, `Artifact` removed)
7. GENCODE v50 plus the passed novels (`XNOVEL_*`, gene from SQANTI's `associated_gene`, SQANTI novel
   genes as `XNOVELG_*` loci; a fusion is led by its 5' gene; `cohorts "tdpkd"`, `"sun"` or `"sun,tdpkd"`
   per novel transcript), GTF + FASTA

QC: samtools, Picard and pychopper into MultiQC 1.25.1; NanoStat in `<data_code>_longread_qc_summary.tsv`; plus `<data_code>_read_length_hist.tsv.gz`.

Conda environments: `isoseq-pipeline` (`isoseq-pipeline.yml`), `pychopper_env`, `bambu_env`,
`isoquant313`, `sqanti3_v6`, `multiqc`, `sqanti-browser` (the tool's
`environment.yml`, then `pip install -e` of a clone of `github.com/ConesaLab/SQANTI-browser` at
commit `8257746`). isomatch is not on bioconda; the pipeline runs the v0.8.0 Linux release binary
(`isomatch-v0.8.0-linux-x86_64.tar.gz` from `github.com/zhengxinchang/isomatch/releases`).

## Outputs

Per cohort, `results/<run_code>/`: alignments, long-read QC (`<data_code>_longread_qc_summary.tsv`,
`_read_length_hist.tsv.gz`, `multiqc/`) and the assembler GTFs (`bambu/`, `isoquant/`, `stringtie3/`).

Pooled, `references/novel_isoform_references/<run_code>/`: `isoforms_<run_code>.{gtf,fa}.gz`,
`cross_cohort_novel_provenance.tsv` (one row per novel model: gene, class code, tools, cohorts, `<cohort>_<tool>`
and group flags), `tracking_numbers_final_numbers.tsv` (pooled funnel, `pooled_` metrics) and
`tracking_numbers_<cohort>.tsv` (each cohort's share), `isoforms_<run_code>_sanity_genes.tsv{,.events.tsv}`
(the six TDP-43 positive controls, each found by its junction and traced stage by stage, with the tools and
cohorts that built it: knockdown UNC13A and STMN2 cryptic exons and KCNQ2 exon 5 skip; overexpression HYOU1,
NUP93 and XPNPEP1 exon skips), `consensus/` (merge tables) and `sqanti/`.

UCSC track hub: `sqanti/isoforms_<run_code>_ucsc_hub/` holds only the novel models that
reached SQANTI, kept and removed. UCSC reads hubs only over public HTTP(S), so copy the whole folder
to a public web server (GitHub raw URL, lab server, S3) and load the URL of its `hub.txt` under
My Data > Track Hubs > Connected Hubs. `hg38/trix.ix` (the attribute search index) was 124 MB for tdpkd v9,
above GitHub's 100 MB file limit. Models on unplaced scaffolds keep GENCODE names (`GL000195.1`),
which UCSC hg38 does not know (`chrUn_GL000195v1`), so UCSC does not show them.
