# isoquant discovery on the pooled long-read bam of one group (Prjibelski et al. 2023). --illumina_bam
# gives isoquant the group's short reads to confirm and correct splice junctions. --clean_start stops
# isoquant reusing another group's converted gene database. after the run the converted gene
# database (*.db) and the large per-read assignment table are deleted to save space
rule run_isoquant:
    input:
        bams      = lambda wc: pooled_bam(wc.group),
        bai       = lambda wc: pooled_bam(wc.group) + ".bai",
        ref       = genome,
        gtf       = ref_gtf,
        short_bam = lambda wc: short_bam(wc.group),
        short_bai = lambda wc: short_bam(wc.group) + ".bai"
    output:
        counts = out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.discovered_transcript_grouped_file_name_counts.tsv",
        tpm    = out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.discovered_transcript_grouped_file_name_tpm.tsv",
        gtf    = temp(out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf")
    params:
        outdir = out_folder + "isoquant/{group}/",
        prefix = data_code + "_{group}"
    threads: 16
    shell:
        """
        conda activate isoquant313
        ml samtools
        ml bedtools
        isoquant --data_type nanopore \
        --bam {input.bams} \
        --illumina_bam {input.short_bam} \
        --reference {input.ref} \
        --complete_genedb --genedb {input.gtf} --clean_start \
        --model_construction_strategy {config[isoquant_model_strategy]} \
        --output {params.outdir} --prefix {params.prefix} \
        --read_group file_name \
        --threads {threads}
        rm -f {params.outdir}*.db {params.outdir}{params.prefix}/{params.prefix}.read_assignments.tsv.gz
        """
