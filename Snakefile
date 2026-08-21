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

#####
## Config
#####

metadata     = config['metadata']
prep         = config['prep']
pacbio_input = config.get('pacbio_input', 'flnc')
run_code     = config['run_code']
out_folder   = config['out_folder'] + run_code + "/"
data_code    = config['data_code']
ref_genome   = config['ref_genome']   # path without the .fa suffix
ref_fasta    = config['ref_fasta']
reflat_file  = config['reflat_file']
junctions_bed   = config['junctions_bed']
junction_folder = config['junction_folder']   # short-read SJ.out.tab files

# reference gtf is gzipped in config; rule decompress_ref_gtf writes the plain copy
# that gffcompare, stringtie and collapse_annotation need
ref_gtf_source = config['ref_gtf']
assert ref_gtf_source.endswith(".gz"), "config 'ref_gtf' must point to a gzipped .gtf.gz"
ref_gtf = out_folder + data_code + ".reference_annotation.gtf"

referenceFa  = ref_genome + ".fa"
referenceGTF = ref_gtf

# one transcript per gene, for RNA-SeQC. built from the reference if the config path is missing
genes_gtf_source = config['genes_gtf']
genes_gtf        = out_folder + data_code + ".genes.gtf"

# short-read junctions passed to minimap2 --junc-bed. nanopore only, pbmm2 ignores it
use_short_read_junctions       = config.get('use_short_read_junctions', True)
short_read_junc_min_uniq       = config.get('short_read_junc_min_uniq', 2)
short_read_junc_canonical_only = config.get('short_read_junc_canonical_only', True)

# pychopper orients and trims full-length ONT cDNA reads before alignment
use_pychopper = config.get('use_pychopper', True)
pychopper_kit = config.get('pychopper_kit', "PCS111")
pychopper_q   = config.get('pychopper_q', None)   # unset means autotune the cutoff

pacbio_primers = config.get('pacbio_primers', "reference/NEB_primers_01_2019.fa")

#####
## Samples and groups
#####

if metadata.endswith(".tsv"):
    metaDF = pd.read_csv(metadata, sep='\t')
elif metadata.endswith(".xlsx"):
    metaDF = pd.read_excel(metadata)
else:
    raise ValueError("metadata must be .tsv or .xlsx: " + metadata)

samples       = metaDF['sample'].astype(str).tolist()
metadata_dict = metaDF.set_index("sample").T.to_dict()

# assemblers are run per group, not per sample. membership is by sample-name prefix (GFP_1 -> GFP)
groups = config.get('groups', [])
sample_groups = {g: [s for s in samples if s.startswith(g)] for g in groups}
for g, ss in sample_groups.items():
    if not ss:
        print("WARNING: no samples matched group '%s'" % g)

print("Prep:", prep)
print("Samples:", samples)
print("Groups:", sample_groups)


def group_bams(group):
    return [out_folder + "%s/alignment/%s.aligned.cram" % (s, s) for s in sample_groups[group]]


# the single pooled input shared by all three assemblers
def pooled_bam(group):
    return out_folder + "pooled/" + data_code + "_%s.pooled.bam" % group


def pooled_bai(group):
    return out_folder + "pooled/" + data_code + "_%s.pooled.bam.bai" % group


# per-sample FLNC for the pacbio aligner, either from ccs or straight from metadata
def pacbio_flnc(sample):
    if pacbio_input == "subreads":
        return out_folder + "%s/flnc/%s.flnc.bam" % (sample, sample)
    return metadata_dict[sample]["flnc_bam_path"]


# minimap2 index, prebuilt from config or built by create_index
if config.get('ref_genome_index'):
    mmi = config['ref_genome_index']
else:
    mmi = ref_genome + ".mmi"

#####
## Sub-workflows. these carry no shell.prefix, config or rule all of their own
#####

include: "workflows/ccs_pipeline.smk"   # pacbio only
include: "workflows/qc_pipeline.smk"
include: "workflows/bambu_pipeline.smk"
include: "workflows/isoquant_pipeline.smk"
include: "workflows/stringtie_pipeline.smk"
include: "workflows/consensus_pipeline.smk"
include: "workflows/sqanti_pipeline.smk"
# include: "workflows/lapa_pipeline.smk"   # apa branch, currently disabled

######
## Main pipeline definition
######

rule all:
    input:
        expand(out_folder + "{sample}/alignment/{sample}.aligned.cram",      sample=samples),
        expand(out_folder + "{sample}/alignment/{sample}.aligned.cram.crai", sample=samples),
        out_folder + data_code + "_read_lengths_collated.tsv.gz",
        out_folder + data_code + "_nanostat_collated.tsv",
        out_folder + "multiqc/multiqc_report.html",

        expand(out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf", group=groups),
        expand(out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf", group=groups),
        expand(out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf", group=groups),

        out_folder + "consensus/" + data_code + ".isoforms.gtf.gz",
        out_folder + "consensus/" + data_code + ".isoforms.fa.gz",
        out_folder + "consensus/" + data_code + "_tool_overlap.tsv",
        out_folder + "consensus/" + data_code + "_novel_provenance.tsv",
        out_folder + "consensus/" + data_code + "_sanity_genes.tsv",
        out_folder + "consensus/" + data_code + "_longread_specific.tsv",
        out_folder + "consensus/" + data_code + "_tracking_numbers.tsv",
        out_folder + "sqanti/" + data_code + "_filter_sqanti_classification.tsv",


rule decompress_ref_gtf:
    input:
        ref_gtf_source
    output:
        temp(ref_gtf)
    shell:
        "gunzip -c {input} > {output}"


# the run-local copy keeps this out of the shared reference folder
_genes_prepare_cmd = ("gunzip -c {input} > {output}"
                      if genes_gtf_source.endswith(".gz") else "cp {input} {output}")

if os.path.exists(genes_gtf_source):

    rule prepare_genes_gtf:
        input:
            genes_gtf_source
        output:
            temp(genes_gtf)
        shell:
            _genes_prepare_cmd

else:

    rule collapseAnnotationGlobal:
        input:
            referenceGTF
        output:
            temp(genes_gtf)
        params:
            script = "scripts/collapse_annotation.py",
            python = "/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/python"
        shell:
            "{params.python} {params.script} {input} {output}"


rule pool_group_bams:
    input:
        bams = lambda wc: group_bams(wc.group)
    output:
        bam = temp(out_folder + "pooled/" + data_code + "_{group}.pooled.bam"),
        bai = temp(out_folder + "pooled/" + data_code + "_{group}.pooled.bam.bai")
    params:
        ref = referenceFa
    threads: 4
    shell:
        """
        conda activate isoseq-pipeline
        samtools merge -f -@ {threads} --reference {params.ref} {output.bam} {input.bams}
        samtools index {output.bam}
        """
