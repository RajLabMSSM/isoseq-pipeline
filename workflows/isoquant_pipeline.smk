# IsoQuant (Prjibelski et al. 2023) discovery, run once per group
# pacbio uses pacbio_ccs with --fl_data, everything else runs as nanopore

isoquant         = "isoquant"   # on PATH from the isoquant313 conda env
isoquant_threads = 16

if prep.startswith("pacbio"):
    isoquant_data_type = "pacbio_ccs"
    isoquant_fl_flag   = "--fl_data"
    _default_strategy  = "default_pacbio"
else:
    isoquant_data_type = "nanopore"
    isoquant_fl_flag   = ""
    _default_strategy  = "default_ont"

isoquant_model_strategy = config.get("isoquant_model_strategy", _default_strategy)

# report novel single-exon transcripts. off by default for nanopore, so only passed when set
isoquant_report_novel_unspliced = config.get("isoquant_report_novel_unspliced", None)

# polish long-read splice sites against the short-read bam already pooled for stringtie --mix
isoquant_illumina = config.get("isoquant_illumina_correction", False)


def isoquant_run_inputs(wildcards):
    d = {
        "bams": pooled_bam(wildcards.group),
        "bai":  pooled_bai(wildcards.group),
        "ref":  referenceFa,
        "gtf":  ref_gtf,
    }
    if isoquant_illumina:
        short = out_folder + "stringtie3/%s/%s_%s.pooled.short.bam" % (wildcards.group, data_code, wildcards.group)
        d["short_bam"] = short
        d["short_bai"] = short + ".bai"
    return d


# --clean_start stops isoquant reusing another group's converted gene database
rule run_isoquant:
    input:
        unpack(isoquant_run_inputs)
    output:
        counts = out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.discovered_transcript_grouped_file_name_counts.tsv",
        tpm    = out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.discovered_transcript_grouped_file_name_tpm.tsv",
        gtf    = temp(out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf")
    params:
        outdir    = out_folder + "isoquant/{group}/",
        prefix    = data_code + "_{group}",
        illumina  = lambda wc: ("--illumina_bam " + out_folder + "stringtie3/%s/%s_%s.pooled.short.bam" % (wc.group, data_code, wc.group)) if isoquant_illumina else "",
        unspliced = ("--report_novel_unspliced " + str(isoquant_report_novel_unspliced).lower()) if isoquant_report_novel_unspliced is not None else ""
    threads: isoquant_threads
    shell:
        """
        conda activate isoquant313
        ml samtools
        ml bedtools
        {isoquant} --data_type {isoquant_data_type} {isoquant_fl_flag} \
        --bam {input.bams} \
        {params.illumina} \
        --reference {input.ref} \
        --complete_genedb --genedb {input.gtf} --clean_start \
        --model_construction_strategy {isoquant_model_strategy} \
        {params.unspliced} \
        --output {params.outdir} --prefix {params.prefix} \
        --read_group file_name \
        --threads {threads}
        rm -f {params.outdir}*.db {params.outdir}{params.prefix}/{params.prefix}.read_assignments.tsv.gz
        """
