# convert the reference gtf once into bambu's annotation object (prepareAnnotations)
rule bambu_create_annotation:
    input:
        gtf = ref_gtf,
        # the script is an input, so editing it re-runs this step
        script = "scripts/bambu_annotation.R"
    output:
        rdata = out_folder + "bambu/bambu_annotation.RData"
    shell:
        """
        conda activate bambu_env
        Rscript {input.script} -i {input.gtf} -o {output.rdata}
        """


# discovery and quantification on the pooled bam of one group (Chen et al. 2023). NDR = novel discovery rate: the share
# of unannotated models among all candidates ranked at or above a model's score. one bam runs on one core
rule run_bambu:
    input:
        bams = lambda wc: pooled_bam(wc.group),
        bai  = lambda wc: pooled_bam(wc.group) + ".bai",
        anno = out_folder + "bambu/bambu_annotation.RData",
        # the script is an input, so editing it re-runs this step
        script = "scripts/bambu_run.R"
    output:
        gtf    = out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf",
        counts = out_folder + "bambu/{group}/" + data_code + "_{group}_counts_transcript.txt",
        rdata  = out_folder + "bambu/{group}/" + data_code + "_{group}_bambu.RData"
    params:
        prefix = out_folder + "bambu/{group}/" + data_code + "_{group}"
    threads: 1
    shell:
        """
        conda activate bambu_env
        Rscript {input.script} \
        --cores {threads} \
        --fasta {genome} \
        --anno {input.anno} \
        --prefix {params.prefix} \
        --NDR {config[bambu_ndr]} \
        {input.bams}
        """
