# isoquant_pipeline.smk
# IsoQuant 3.4.1 (Prjibelski et al. 2023) transcript discovery, run PER GROUP =
# Join & Call (Jetzinger et al. 2025). All of a group's aligned BAMs go into one
# IsoQuant run: model construction pools every read (it does not partition discovery
# by read group), while `--read_group file_name` still gives per-sample counts.
#
# Platform is config-driven via `prep`: ONT -> --data_type nanopore (k14/splice);
# PacBio -> --data_type pacbio_ccs + --fl_data (Jetzinger used --fl_data for PacBio
# only). Included by Snakefile (no shell.prefix / config / rule all here).

isoquant         = "/sc/arion/projects/ad-omics/data/software/IsoQuant/isoquant.py"
isoquant_threads = 16

if prep.startswith("pacbio"):
    isoquant_data_type = "pacbio_ccs"
    isoquant_fl_flag   = "--fl_data"     # full-length, PacBio HiFi
else:
    isoquant_data_type = "nanopore"
    isoquant_fl_flag   = ""              # ONT: no --fl_data


rule isoquant_run:
    """
    Join & Call within a group. IsoQuant writes to
    {out}/isoquant/{group}/{data_code}_{group}/{data_code}_{group}.* (the inner
    folder is IsoQuant's --prefix/experiment name).
    """
    input:
        bams = lambda wc: group_bams(wc.group),
        ref  = referenceFa,
        gtf  = ref_gtf
    output:
        counts = out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.transcript_model_grouped_counts.tsv",
        tpm    = out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.transcript_model_grouped_tpm.tsv",
        gtf    = out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf"
    params:
        outdir = out_folder + "isoquant/{group}/",
        prefix = data_code + "_{group}"
    threads: isoquant_threads
    shell:
        "conda activate isoquant; ml samtools; ml bedtools; "
        "{isoquant} --data_type {isoquant_data_type} {isoquant_fl_flag} "
        "--bam {input.bams} "
        "--reference {input.ref} "
        "--complete_genedb --genedb {input.gtf} "
        "--output {params.outdir} --prefix {params.prefix} "
        "--read_group file_name "
        "--threads {threads}"
