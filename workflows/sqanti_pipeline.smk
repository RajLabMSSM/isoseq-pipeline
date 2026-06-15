# sqanti_pipeline.smk
# SQANTI3 6.0.1 QC + lab rules-filter on the per-group consensus GTF (2nd
# precision filter; junctions vs cohort STAR SJ.out.tab). Env sqanti3_v6 (bioconda).

sqanti_threads = 12


# -c globs ALL cohort SJ.out.tab (pattern passed literally for SQANTI to expand)
rule sqanti_qc:
    input:
        gtf = out_folder + "consensus/{group}/" + data_code + "_{group}_consensus.gtf",
        ref = ref_gtf,
        fa  = referenceFa
    output:
        classification = out_folder + "sqanti/{group}/" + data_code + "_{group}_classification.txt",
        corrected_gtf  = out_folder + "sqanti/{group}/" + data_code + "_{group}_corrected.gtf",
        corrected_fa   = out_folder + "sqanti/{group}/" + data_code + "_{group}_corrected.fasta"
    params:
        outdir = out_folder + "sqanti/{group}/",
        prefix = data_code + "_{group}",
        junc   = junction_folder + "/*SJ.out.tab"
    threads: sqanti_threads
    shell:
        "conda activate sqanti3_v6; "
        "sqanti3_qc.py --isoforms {input.gtf} --refGTF {input.ref} --refFasta {input.fa} "
        "-c '{params.junc}' "
        "-d {params.outdir} -o {params.prefix} -t {threads}"


# counts-free (quantification is downstream)
rule sqanti_filter:
    input:
        classification = out_folder + "sqanti/{group}/" + data_code + "_{group}_classification.txt",
        fasta          = out_folder + "sqanti/{group}/" + data_code + "_{group}_corrected.fasta",
        gff            = out_folder + "sqanti/{group}/" + data_code + "_{group}_corrected.gtf"
    output:
        gtf            = out_folder + "sqanti/{group}/" + data_code + "_{group}_filter_sqanti.cds.gtf",
        classification = out_folder + "sqanti/{group}/" + data_code + "_{group}_filter_sqanti_classification.tsv",
        fasta          = out_folder + "sqanti/{group}/" + data_code + "_{group}_filter_sqanti.fasta"
    params:
        script     = "scripts/filter_sqanti.R",
        out_prefix = out_folder + "sqanti/{group}/" + data_code + "_{group}"
    shell:
        "conda activate bambu_env; "
        "Rscript {params.script} "
        "--sqanti {input.classification} --fasta {input.fasta} --gff {input.gff} "
        "--output {params.out_prefix}"
