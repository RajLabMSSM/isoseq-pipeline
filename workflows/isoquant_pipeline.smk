# isoquant_pipeline.smk
# IsoQuant 3.13.0 (Prjibelski et al. 2023) discovery, per group (Join & Call).
# Platform via prep: ONT --data_type nanopore; PacBio pacbio_ccs + --fl_data.

isoquant         = "isoquant"        # on PATH from the isoquant313 env (bioconda)
isoquant_threads = 16

if prep.startswith("pacbio"):
    isoquant_data_type = "pacbio_ccs"
    isoquant_fl_flag   = "--fl_data"     # full-length, PacBio HiFi
else:
    isoquant_data_type = "nanopore"
    isoquant_fl_flag   = ""

# --read_group file_name gives per-sample counts; discovery still pools all reads
rule run_isoquant:
    input:
        bams = lambda wc: pooled_bam(wc.group),
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
        "conda activate isoquant313; ml samtools; ml bedtools; "
        "{isoquant} --data_type {isoquant_data_type} {isoquant_fl_flag} "
        "--bam {input.bams} "
        "--reference {input.ref} "
        "--complete_genedb --genedb {input.gtf} "
        "--output {params.outdir} --prefix {params.prefix} "
        "--read_group file_name "
        "--threads {threads}"
