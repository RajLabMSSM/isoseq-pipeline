Long-read isoform discovery pipeline. Assembles a transcriptome from ONT or PacBio reads with three
assemblers, takes the end-aware union of what they find, filters it with SQANTI3, and emits a
reference plus novel annotation for downstream quantification.

See links to instruction sections below:

- [How to make a metadata sheet](#how-to-make-a-metadata-sheet)
- [How to run the pipeline](#how-to-run-the-pipeline)
- [What the pipeline does](#what-the-pipeline-does)
- [Outputs](#outputs)
- [Combining cohorts](#combining-cohorts)

The only things you should need to edit are the config yaml and the metadata sheet. The pipeline code
lives here once; a cohort directory holds only its config, metadata and logs.

### How to make a metadata sheet

A tab separated `.tsv` (or `.xlsx`). The necessary columns are

`sample` `long_read_fastq`

and, if you are running StringTie in hybrid mode or IsoQuant with short-read correction,

`short_read_bam_path`

#### What to put in each column

`sample` - the name you want the sample to have. No spaces, don't start with a number. The group a
sample belongs to is taken from the **prefix** of this name, so `TDP_1` and `TDP_2` both land in group
`TDP`. The groups themselves are listed in the config.

`long_read_fastq` - full path to the sample's long reads. Gzipped is fine.

`short_read_bam_path` - full path to the matching short-read BAM. Leave as `NA` if the sample has none;
samples without one are skipped when the short-read BAMs are pooled.

For PacBio, supply `flnc_bam_path` instead of `long_read_fastq`, or `subreads_xml` if you are starting
from subreads and want the pipeline to run ccs itself.

sample | long_read_fastq | short_read_bam_path
-- | -- | --
GFP_1 | /path/to/GFP_1.fastq.gz | /path/to/GFP_1.bam
GFP_2 | /path/to/GFP_2.fastq.gz | /path/to/GFP_2.bam
TDP_1 | /path/to/TDP_1.fastq.gz | /path/to/TDP_1.bam
TDP_2 | /path/to/TDP_2.fastq.gz | /path/to/TDP_2.bam

### How to run the pipeline

Run from the cohort's own directory, pointing at the Snakefile here. `run_code` in the config names the
output folder, so bumping it starts a clean run without touching the previous one.

**Dry run first.**

```
sh ~/bin/snakejob -n -s Snakefile -c <cohort>_config.yaml -a acc_als-omics
```

**Submit to the cluster.**

```
snakejob_HPC -s Snakefile -c <cohort>_config.yaml
```

Add `-r` to resume. This passes `--rerun-triggers mtime` to snakemake, so completed upstream stages are
compared on timestamps only and are not redone because a rule's text changed.

**NB: `snakejob_HPC -n` does not dry run, it submits.** The `-n` is swallowed by the wrapper. Use the
inner `snakejob` as above, which runs snakemake locally with no `--cluster` and so cannot submit.

#### Minimal dependencies

Conda environments, all expected to exist already: `isoseq-pipeline` (the default, set by
`shell.prefix`), `bambu_env`, `isoquant313`, `sqanti3_v6`, `pychopper_env`. StringTie, gffcompare and
minimap2 are called by absolute path.

### What the pipeline does

1. Orient and trim full-length cDNA reads with pychopper (ONT cDNA only)
2. Build a combined `--junc-bed` from the annotation and the short-read STAR junctions
3. Align with minimap2, or pbmm2 for PacBio, and store alignments as CRAM
4. QC with samtools, RNA-SeQC, Picard and NanoStat, collected into one MultiQC report
5. Pool each group's alignments into a single BAM, the shared input to all three assemblers
6. Discover isoforms per group with Bambu, IsoQuant and StringTie3
7. Take the end-aware union of the three, keeping alternative polyA isoforms distinct
8. Filter the novel isoforms with SQANTI3 QC and the rules filter
9. Emit the reference plus the surviving novel isoforms, as GTF and transcript FASTA

Assembly is per group rather than per sample, so a lowly expressed isoform seen once in several samples
is still assembled. Steps 7 and 8 split the work deliberately: the union errs towards keeping things,
and SQANTI is the only precision filter.

### Outputs

Everything lands in `<out_folder>/<run_code>/`. The files you are likely to want:

file | what
-- | --
`consensus/<data_code>.isoforms.gtf.gz` | the reference plus novel annotation
`consensus/<data_code>.isoforms.fa.gz` | matching transcript fasta
`consensus/<data_code>_tracking_numbers.tsv` | the funnel, raw per tool through to final novel
`consensus/<data_code>_novel_provenance.tsv` | which tools support each novel isoform
`consensus/<data_code>_tool_overlap.tsv` | overlap between the three assemblers
`consensus/<data_code>_sanity_genes.tsv` | are the expected genes present, and how many novel isoforms landed on them
`consensus/<data_code>_sanity_genes.tsv.events.tsv` | per known cryptic event, is it in the final gtf and if not which stage dropped it
`consensus/<data_code>_longread_specific.tsv` | novel isoforms carrying an intron absent from the short-read junctions
`multiqc/multiqc_report.html` | QC across all samples

Counting novel isoforms alone is not a valid recovery test, since an event already annotated in the
reference is scored as zero novel. Use the events table for that.

### Combining cohorts

Once every cohort has finished, `cross_cohort/` merges them into one annotation with the same collapse
tolerances used within a cohort.

```
cd cross_cohort
sh ~/bin/snakejob -n -s Snakefile -c cross_cohort_config.yaml -a acc_als-omics
snakejob_HPC -s Snakefile -c cross_cohort_config.yaml
```

The collapse settings in `cross_cohort_config.yaml` must match the per-cohort config. If they disagree,
the cross-cohort key is stricter than the within-cohort one and a single isoform is split across
several entries.

Three of its outputs are not named after the run and are overwritten in place:
`cross_cohort_novel_provenance.tsv`, `cross_cohort_global_provenance.tsv` and
`tracking_numbers_final_numbers.tsv`. Back them up before re-running if you want to keep the previous
versions.
