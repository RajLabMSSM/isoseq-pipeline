
# Snakefile
# isoseq-pipeline — ONT / PacBio long-read RNA-seq isoform discovery
# (quite literally not isoseq3; used for ONT nanopore long-read RNA-seq. alas!)
#
# Strategy: Join & Call, pooled PER GROUP (Jetzinger et al. 2025) so that
# condition-specific rare novel isoforms (e.g. TDP) are not diluted by the other
# group. QC -> {bambu, isoquant, stringtie3} all run from this single Snakefile in
# one `snakemake` submission; Snakemake sequences them by file dependency.
# Brooke Friedman

# ── single global shell prefix (conda setup only) ─────────────────────────────
# Mirrors old/Snakefile: ONE prefix. It sets conda up and activates the default
# QC env (isoseq-pipeline). Downstream rules switch env inline with
# `conda activate ...`. shell.prefix is GLOBAL in Snakemake (only one can be
# active), so the included workflow .smk files do NOT define their own — they
# inherit this one and override per-rule where needed.
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

# ── config (parsed ONCE here; the included .smk inherit every name below) ──────
metadata        = config['metadata']
prep            = config['prep']                  # nanopore_cdna | nanopore_direct | pacbio
run_code        = config['run_code']
out_folder      = config['out_folder'] + run_code + "/"
data_code       = config['data_code']
ref_genome      = config['ref_genome']            # path WITHOUT .fa
ref_gtf         = config['ref_gtf']
ref_fasta       = config['ref_fasta']             # transcript fasta (downstream combine steps)
junctions_bed   = config['junctions_bed']
junction_folder = config['junction_folder']       # short-read SJ.out.tab files for SQANTI
genes_gtf       = config['genes_gtf']             # collapsed single-tx-per-gene GTF (RNA-SeQC)
reflat_file     = config['reflat_file']           # Picard refFlat

referenceFa  = ref_genome + ".fa"
referenceGTF = ref_gtf

# ── metadata ──────────────────────────────────────────────────────────────────
if metadata.endswith(".tsv"):
    metaDF = pd.read_csv(metadata, sep='\t')
elif metadata.endswith(".xlsx"):
    metaDF = pd.read_excel(metadata)
else:
    raise ValueError("metadata must be .tsv or .xlsx: " + metadata)

samples       = metaDF['sample'].astype(str).tolist()
metadata_dict = metaDF.set_index("sample").T.to_dict()

# ── groups = the Join & Call unit ─────────────────────────────────────────────
# Membership: sample name starts with a group string in config['groups']
# (e.g. GFP_1 -> GFP, TDP1 -> TDP, CTRL3 -> CTRL).
groups = config.get('groups', [])
sample_groups = {g: [s for s in samples if s.startswith(g)] for g in groups}
for g, ss in sample_groups.items():
    if not ss:
        print("WARNING: no samples matched group '%s'" % g)

print("Prep:", prep)
print("Samples:", samples)
print("Groups:", sample_groups)

# helper: aligned BAMs belonging to a group (used by all three assemblers)
def group_bams(group):
    return [out_folder + "%s/alignment/%s.aligned.bam" % (s, s) for s in sample_groups[group]]

# ── minimap2 index ────────────────────────────────────────────────────────────
# Use a prebuilt index if given in config, else build it (create_index in qc).
if config.get('ref_genome_index'):
    mmi = config['ref_genome_index']
else:
    mmi = ref_genome + ".mmi"

# ── sub-workflows ─────────────────────────────────────────────────────────────
# Each .smk has NO shell.prefix / config parsing / rule all of its own; it inherits
# everything above and contributes only its tool-specific constants + rules.
include: "workflows/qc_pipeline.smk"
include: "workflows/bambu_pipeline.smk"
include: "workflows/isoquant_pipeline.smk"
include: "workflows/stringtie_pipeline.smk"


# ══════════════════════════════════════════════════════════════════════════════
# RULE ALL  (single target rule for the whole pipeline)
# Assembler-only scope for now: QC + bambu/isoquant/stringtie3 per-group outputs.
# Downstream SQANTI/filter/combine are out of scope here (see git history / old/).
# ══════════════════════════════════════════════════════════════════════════════
rule all:
    input:
        # ── QC ───────────────────────────────────────────────────────────────
        expand(out_folder + "{sample}/alignment/{sample}.aligned.bam",     sample=samples),
        expand(out_folder + "{sample}/alignment/{sample}.aligned.bam.bai", sample=samples),
        out_folder + data_code + "_read_lengths_collated.tsv.gz",
        out_folder + "multiqc/multiqc_report.html",

        # ── Bambu  (per-group Join & Call) ───────────────────────────────────
        expand(out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf", group=groups),
        expand(out_folder + "bambu/{group}/" + data_code + "_{group}_counts_transcript.txt",    group=groups),

        # ── IsoQuant  (per-group Join & Call) ────────────────────────────────
        expand(out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf",                group=groups),
        expand(out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.transcript_model_grouped_counts.tsv",    group=groups),

        # ── StringTie3  (per-group Join & Call: pooled-BAM assembly) ──────────
        expand(out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf", group=groups),


# ══════════════════════════════════════════════════════════════════════════════
# REFERENCE PREP  (run once per reference genome)
# ══════════════════════════════════════════════════════════════════════════════
rule collapseAnnotationGlobal:
    """
    Collapse reference GTF to one transcript per gene (required by RNA-SeQC).
    Output path mirrors what qc_pipeline.smk expects as genes_gtf.
    """
    input:  referenceGTF
    output: genes_gtf
    params:
        script = "scripts/collapse_annotation.py",
        python = "/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/python"
    shell:
        "{params.python} {params.script} {input} {output}"
