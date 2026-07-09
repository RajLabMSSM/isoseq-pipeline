# Snakefile
# Brooke Friedman
shell.prefix(
    'export PS1=""; '
    'ml anaconda3; '
    'CONDA_BASE=$(conda info --base); '
    'source $CONDA_BASE/etc/profile.d/conda.sh; '
    'module purge; '
    'conda activate isoseq-pipeline; '
)
import pandas as pd
import glob
import os

# config
metadata        = config['metadata']
prep            = config['prep']           
pacbio_input    = config.get('pacbio_input', 'flnc')
run_code        = config['run_code']
out_folder      = config['out_folder'] + run_code + "/"
data_code       = config['data_code']
ref_genome      = config['ref_genome']            # path WITHOUT .fa
ref_gtf         = config['ref_gtf']
ref_fasta       = config['ref_fasta']             
junctions_bed   = config['junctions_bed']
junction_folder = config['junction_folder']       # short-read SJ.out.tab files

# short-read junctions -> minimap2 --junc-bed (nanopore only; pbmm2 ignores it)
use_short_read_junctions   = config.get('use_short_read_junctions', True)
short_read_junc_min_uniq   = config.get('short_read_junc_min_uniq', 2)
short_read_junc_canonical_only = config.get('short_read_junc_canonical_only', True)

# pychopper: orient/trim full-length ONT cDNA reads before alignment (cdna only)
use_pychopper  = config.get('use_pychopper', True)
pychopper_kit  = config.get('pychopper_kit', "PCS111")
genes_gtf       = config['genes_gtf']             
reflat_file     = config['reflat_file']           

referenceFa  = ref_genome + ".fa"
referenceGTF = ref_gtf

# metadata
if metadata.endswith(".tsv"):
    metaDF = pd.read_csv(metadata, sep='\t')
elif metadata.endswith(".xlsx"):
    metaDF = pd.read_excel(metadata)
else:
    raise ValueError("metadata must be .tsv or .xlsx: " + metadata)

samples       = metaDF['sample'].astype(str).tolist()
metadata_dict = metaDF.set_index("sample").T.to_dict()

# groups = Join & Call, membership by sample-name prefix (GFP_1 -> GFP)
groups = config.get('groups', [])
sample_groups = {g: [s for s in samples if s.startswith(g)] for g in groups}
for g, ss in sample_groups.items():
    if not ss:
        print("WARNING: no samples matched group '%s'" % g)

print("Prep:", prep)
print("Samples:", samples)
print("Groups:", sample_groups)

# a group's aligned BAMs (pool inputs)
def group_bams(group):
    return [out_folder + "%s/alignment/%s.aligned.bam" % (s, s) for s in sample_groups[group]]

# per-group pooled BAM: the single Join & Call input shared by all three assemblers
def pooled_bam(group):
    return out_folder + "pooled/" + data_code + "_%s.pooled.bam" % group

def pooled_bai(group):
    return out_folder + "pooled/" + data_code + "_%s.pooled.bam.bai" % group

# per-sample FLNC for the pacbio aligner: CCS output (subreads) or metadata (flnc)
def pacbio_flnc(sample):
    if pacbio_input == "subreads":
        return out_folder + "%s/flnc/%s.flnc.bam" % (sample, sample)
    return metadata_dict[sample]["flnc_bam_path"]

# pacbio cDNA primers for lima/refine (CCS path)
pacbio_primers  = config.get('pacbio_primers', "reference/NEB_primers_01_2019.fa")

# minimap2 index: prebuilt from config, else built by create_index
if config.get('ref_genome_index'):
    mmi = config['ref_genome_index']
else:
    mmi = ref_genome + ".mmi"

# sub-workflows (no shell.prefix / config / rule all of their own)
include: "workflows/ccs_pipeline.smk"          # only needed for pacbio
include: "workflows/qc_pipeline.smk"
include: "workflows/bambu_pipeline.smk"
include: "workflows/isoquant_pipeline.smk"
include: "workflows/stringtie_pipeline.smk"
include: "workflows/consensus_pipeline.smk"
include: "workflows/sqanti_pipeline.smk"


rule all:
    input:
        expand(out_folder + "{sample}/alignment/{sample}.aligned.bam",     sample=samples),
        expand(out_folder + "{sample}/alignment/{sample}.aligned.bam.bai", sample=samples),
        out_folder + data_code + "_read_lengths_collated.tsv.gz",
        out_folder + data_code + "_nanostat_collated.tsv",
        out_folder + "multiqc/multiqc_report.html",

        expand(out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf",                group=groups),
        expand(out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf", group=groups),
        expand(out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf",                     group=groups),
        out_folder + "consensus/" + data_code + "_reference_plus_novel.gtf.gz",
        out_folder + "consensus/" + data_code + "_reference_plus_novel.fa.gz",
        out_folder + "consensus/" + data_code + "_tool_overlap.tsv",
        out_folder + "consensus/" + data_code + "_novel_provenance.tsv",
        out_folder + "consensus/" + data_code + "_sanity_genes.tsv",
        out_folder + "sqanti/" + data_code + "_filter_sqanti_classification.tsv",


# collapse reference GTF to one transcript per gene (RNA-SeQC input)
rule collapseAnnotationGlobal:
    input:  referenceGTF
    output: genes_gtf
    params:
        script = "scripts/collapse_annotation.py",
        python = "/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/python"
    shell:
        "{params.python} {params.script} {input} {output}"


# pooled BAM (Join & Call): merge a group's aligned BAMs -> one sorted/indexed BAM.
rule pool_group_bams:
    input:
        bams = lambda wc: group_bams(wc.group)
    output:
        bam = temp(out_folder + "pooled/" + data_code + "_{group}.pooled.bam"),
        bai = temp(out_folder + "pooled/" + data_code + "_{group}.pooled.bam.bai")
    threads: 4
    shell:
        "conda activate isoseq-pipeline; "
        "samtools merge -f -@ {threads} {output.bam} {input.bams}; "
        "samtools index {output.bam}"
