# bambu_pipeline.smk
# Bambu (Chen et al. 2023) discovery, per group (Join & Call). NDR knob: higher
# -> more novel; see scripts/bambu_run.R.

bambu_cores = 8          # bambu 3.12.1 + its own R, from the bambu_env conda env

# reference annotation, shared across groups
rule bambu_create_annotation:
    input:
        gtf = ref_gtf
    output:
        rdata = out_folder + "bambu/bambu_annotation.RData"
    params:
        script = "scripts/bambu_annotation.R"
    shell:
        "conda activate bambu_env; "
        "Rscript {params.script} -i {input.gtf} -o {output.rdata}"


# NDR pinned to 0.5 (permissive): discover generously, then the >=2/3 consensus +
# SQANTI restore precision. Lower toward 0.1 for a high-confidence standalone catalog.
rule run_bambu:
    input:
        bams = lambda wc: pooled_bam(wc.group),
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
        "conda activate bambu_env; "
        "Rscript {params.script} --cores {threads} --fasta {referenceFa} "
        "--anno {input.anno} --prefix {params.prefix} --NDR 0.5 {input.bams}"
