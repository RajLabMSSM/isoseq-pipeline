# isoseq-pipeline

Builds the long-read isoform atlas (`isoforms_<run_code>`): per-cohort isoform discovery from ONT cDNA
reads, then one cross-cohort union with GENCODE. The code and configs are set for `v10`, which has not
been run yet; the last completed build is `v9`.

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

`run_code` in `cross_cohort_config.yaml` sets the cohort inputs and the output folder and names, so a
build stops on missing inputs instead of reusing an older one. `build_reference.sh` launches both
cohorts and chains the union behind them (`-x`: union only).

Metadata columns: `sample`, `condition`, `long_read_fastq`, `short_read_bam_path`. Groups are
sample-name prefixes (`TDP_1` -> `TDP`).

## Steps (v10)

1. pychopper 2.7.10 (`-m edlib`; kit PCS109 tdpkd, PCB114 sun), polyA kept
2. minimap2 2.31 `-ax splice` (all cohorts, default k15 from the index), `--secondary=no -G 500k`,
   `--junc-bed` = GENCODE v50 junctions + the cohort's STAR junctions (>= 2 unique reads, canonical).
   tdpkd drops two short-read libraries via `exclude_samples`, here, in the StringTie short-read pool
   and in the SQANTI junction coverage: TDP-5 (aberrant differentiation) and CTRL-4 (low-depth outlier).
   CTRL4 long reads stay in, so tdpkd has 11 long-read and 10 short-read libraries. TDP2 long reads
   combine both of its sequencing runs
3. Long-read alignments are merged per group with `samtools merge`, the short-read STAR alignments
   likewise. Per group, on the pooled long-read BAM:
   - StringTie 3.0.3 `--mix -G` (pooled short-read BAM first; no `-t`)
   - IsoQuant 3.13.0 `sensitive_ont`, `--complete_genedb --clean_start --large_output none`, long reads only
   - Bambu 3.12.1, NDR 0.5, default discovery gates
4. Union of tool x group models; models in GENCODE (`=`, `c`), single-exon and unstranded models
   dropped; `tama_merge -m 0 -z 100 -d merge_dup`, one `no_cap` source per tool and
   group (no tool vote); gffcompare 0.12.10 for class codes
5. SQANTI3 6.0.1 QC with the cohort's `SJ.out.tab`, CAGE peaks and polyA motifs; `sqanti3_filter.py
   rules -e` with `misc/sqanti_rules.json`
6. SQANTI-browser 1.1.2 (Paniagua et al. 2026, bioRxiv, doi:10.64898/2026.05.25.727625) turns
   `_corrected.gtf` and `_RulesFilter_classification.txt` into a UCSC track hub (hg38, chromosome sizes
   from the genome `.fai`, the SQANTI CAGE peaks as a track, `--hub-name <data_code>`); `hubCheck`
   validates it. Each model's details page carries its SQANTI columns, `filter_result` included
   (`Isoform` kept, `Artifact` removed)
7. GENCODE v50 plus the passed novels (`NOVEL_*`, gene from SQANTI's `associated_gene`; a fusion is
   led by its 5' gene), GTF + FASTA
8. Cross-cohort: the cohorts' novels merged by `tama_merge` at the same tolerances, single-exon
   dropped, gene ids set by `cross_cohort_novel_genes.py` (`XNOVEL_*`, `XNOVELG_*` loci)

QC: samtools, Picard and pychopper into MultiQC 1.25.1; NanoStat in `<data_code>_longread_qc_summary.tsv`; plus `<data_code>_read_length_hist.tsv.gz`.

Conda environments: `isoseq-pipeline` (`isoseq-pipeline.yml`), `pychopper_env`, `bambu_env`,
`isoquant313`, `sqanti3_v6`, `tama_env` (Python 2), `multiqc`, `sqanti-browser` (the tool's
`environment.yml`, then `pip install -e` of a clone of `github.com/ConesaLab/SQANTI-browser` at
commit `8257746`).

## Outputs

`results/<run_code>/consensus/`: `<data_code>.isoforms.{gtf,fa}.gz`, `_novel_provenance.tsv`,
`_tool_overlap.tsv`, `_tracking_numbers.tsv` (the funnel), `_sanity_genes.tsv{,.events.tsv}` (the six TDP-43
positive controls, each found by its junction and traced stage by stage: knockdown UNC13A and STMN2 cryptic exons
and KCNQ2 exon 5 skip; overexpression HYOU1, NUP93 and XPNPEP1 exon skips). Cross-cohort:
`references/novel_isoform_references/<run_code>/isoforms_<run_code>.{gtf,fa}.gz`.

UCSC track hub: `results/<run_code>/sqanti/<data_code>_ucsc_hub/` holds only the novel models that
reached SQANTI, kept and removed. UCSC reads hubs only over public HTTP(S), so copy the whole folder
to a public web server (GitHub raw URL, lab server, S3) and load the URL of its `hub.txt` under
My Data > Track Hubs > Connected Hubs. `hg38/trix.ix` (the attribute search index) was 124 MB for tdpkd v9,
above GitHub's 100 MB file limit. Models on unplaced scaffolds keep GENCODE names (`GL000195.1`),
which UCSC hg38 does not know (`chrUn_GL000195v1`), so UCSC does not show them.
