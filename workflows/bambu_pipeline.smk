rule bambu_create_annotation:
    input:
        gtf = ref_gtf
    output:
        rdata = out_folder + "bambu/bambu_annotation.RData"
    shell:
        """
        conda activate bambu_env
        Rscript scripts/bambu_annotation.R -i {input.gtf} -o {output.rdata}
        """


rule run_bambu:
    input:
        bams = lambda wc: pooled_bam(wc.group),
        bai  = lambda wc: pooled_bam(wc.group) + ".bai",
        anno = out_folder + "bambu/bambu_annotation.RData"
    output:
        gtf    = out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf",
        counts = out_folder + "bambu/{group}/" + data_code + "_{group}_counts_transcript.txt",
        rdata  = out_folder + "bambu/{group}/" + data_code + "_{group}_bambu.RData"
    params:
        prefix = out_folder + "bambu/{group}/" + data_code + "_{group}"
    threads: 8
    shell:
        """
        conda activate bambu_env
        Rscript scripts/bambu_run.R \
        --cores {threads} \
        --fasta {genome} \
        --anno {input.anno} \
        --prefix {params.prefix} \
        --NDR {config[bambu_ndr]} \
        {input.bams}
        """
