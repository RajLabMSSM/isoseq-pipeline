# Bambu (Chen et al. 2023) novel isoform discovery, run once per group

bambu_cores = 8                              # bambu 3.12.1 and its own R live in the bambu_env conda env
bambu_ndr   = config.get("bambu_ndr", 0.5)   # novel discovery rate, higher gives more novel transcripts


# the annotation object is built once and reused by every group
rule bambu_create_annotation:
    input:
        gtf = ref_gtf
    output:
        rdata = out_folder + "bambu/bambu_annotation.RData"
    params:
        script = "scripts/bambu_annotation.R"
    shell:
        """
        conda activate bambu_env
        Rscript {params.script} -i {input.gtf} -o {output.rdata}
        """


# discovery is deliberately generous here, the consensus union and SQANTI do the filtering
rule run_bambu:
    input:
        bams = lambda wc: pooled_bam(wc.group),
        bai  = lambda wc: pooled_bai(wc.group),
        anno = out_folder + "bambu/bambu_annotation.RData"
    output:
        gtf    = out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf",
        counts = out_folder + "bambu/{group}/" + data_code + "_{group}_counts_transcript.txt",
        rdata  = out_folder + "bambu/{group}/" + data_code + "_{group}_bambu.RData"
    params:
        prefix = out_folder + "bambu/{group}/" + data_code + "_{group}",
        script = "scripts/bambu_run.R"
    threads: bambu_cores
    shell:
        """
        conda activate bambu_env
        Rscript {params.script} \
        --cores {threads} \
        --fasta {referenceFa} \
        --anno {input.anno} \
        --prefix {params.prefix} \
        --NDR {bambu_ndr} \
        {input.bams}
        """
