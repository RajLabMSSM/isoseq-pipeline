# run before every shell command. snakemake runs bash in strict mode (set -u), and conda activate
# reads PS1, which a batch job does not set, so PS1 is set empty first. then the conda shell
# functions are loaded, every cluster module is unloaded (module purge) so none shadows a tool of the
# environment, and the pipeline's conda environment is activated
shell.prefix(
    'export PS1=""; '
    'ml anaconda3; '
    'CONDA_BASE=$(conda info --base); '
    'source $CONDA_BASE/etc/profile.d/conda.sh; '
    'module purge; '
    'conda activate isoseq-pipeline; '
)
import os
import pandas as pd

data_code       = config['data_code']
out_folder      = config['out_folder'] + config['run_code'] + "/"
genome          = config['ref_genome'] + ".fa"
ref_gtf         = out_folder + data_code + ".reference_annotation.gtf"
junction_folder = config['junction_folder']
excluded_samples = set(x for x in config.get('exclude_samples', '').split(',') if x)

# sample sheet: one row per sample with its long-read fastq and short-read bam
metaDF        = pd.read_csv(config['metadata'], sep='\t')
samples       = metaDF['sample'].astype(str).tolist()

# lookup sample -> {column: value}, used by the input functions of the workflows
metadata_dict = metaDF.set_index("sample").to_dict(orient="index")

# assemblers run per group; membership is by sample-name prefix (GFP_1 -> GFP)
groups        = config['groups']
sample_groups = {g: [s for s in samples if s.startswith(g)] for g in groups}

# a sample matching no group would be left out of every assembly, and one matching two prefixes
# (e.g. TDP and TDPKD) would be pooled twice, so stop unless each sample has exactly one group
for s in samples:
    matched = [g for g in groups if s.startswith(g)]
    if len(matched) != 1:
        raise ValueError("sample %s matches %d of the groups %s, expected exactly one" % (s, len(matched), groups))


# pooled long-read bam of one group (input of bambu, isoquant and stringtie)
def pooled_bam(group):
    return out_folder + "pooled/" + data_code + "_%s.pooled.bam" % group


# pooled short-read bam of one group (only used for stringtie --mix)
def short_bam(group):
    return out_folder + "stringtie3/%s/%s_%s.pooled.short.bam" % (group, data_code, group)


# The cohort run covers read preparation, alignment, long-read QC and discovery with the three assemblers.
# Everything after discovery (removal of GENCODE transcripts, merging, gffcompare classification, SQANTI and the
# final reference) runs once on both cohorts together in cross_cohort/Snakefile.
include: "workflows/qc_pipeline.smk"
include: "workflows/bambu_pipeline.smk"
include: "workflows/isoquant_pipeline.smk"
include: "workflows/stringtie_pipeline.smk"


# The final outputs of a cohort run: the alignments, the long-read QC, and the GTF of each assembler for each group,
# which are the inputs of the pooled stage.
rule all:
    input:
        expand(out_folder + "{sample}/alignment/{sample}.aligned.cram",      sample=samples),
        expand(out_folder + "{sample}/alignment/{sample}.aligned.cram.crai", sample=samples),
        out_folder + data_code + "_read_length_hist.tsv.gz",
        out_folder + data_code + "_longread_qc_summary.tsv",
        out_folder + "multiqc/multiqc_report.html",
        expand(out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf", group=groups),
        expand(out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf",
               group=groups),
        expand(out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf", group=groups),


rule decompress_ref_gtf:
    input:
        config['ref_gtf']
    output:
        temp(ref_gtf)
    shell:
        "gunzip -c {input} > {output}"


rule pool_group_bams:
    input:
        bams = lambda wc: [out_folder + "%s/alignment/%s.aligned.cram" % (s, s) for s in sample_groups[wc.group]]
    output:
        bam = temp(out_folder + "pooled/" + data_code + "_{group}.pooled.bam"),
        bai = temp(out_folder + "pooled/" + data_code + "_{group}.pooled.bam.bai")
    threads: 4
    shell:
        """
        samtools merge -f -@ {threads} --reference {genome} {output.bam} {input.bams}
        samtools index {output.bam}
        """
