# sqanti_pipeline.smk
# SQANTI3 6.0.1 QC + lab rules-filter on the merged per-cohort consensus GTF (2nd
# precision filter; junctions vs cohort STAR SJ.out.tab). Env sqanti3_v6 (bioconda).
# One run per cohort -> one filtered novel-isoform GTF + FASTA per cohort.

sqanti_threads = 12


# -c globs ALL cohort SJ.out.tab (pattern passed literally for SQANTI to expand)
rule sqanti_qc:
    input:
        gtf = out_folder + "consensus/" + data_code + "_merged_consensus.gtf",
        ref = ref_gtf,
        fa  = referenceFa
    output:
        classification = out_folder + "sqanti/" + data_code + "_classification.txt",
        corrected_gtf  = out_folder + "sqanti/" + data_code + "_corrected.gtf",
        corrected_fa   = out_folder + "sqanti/" + data_code + "_corrected.fasta"
    params:
        outdir = out_folder + "sqanti/",
        prefix = data_code,
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
        classification = out_folder + "sqanti/" + data_code + "_classification.txt",
        fasta          = out_folder + "sqanti/" + data_code + "_corrected.fasta",
        gff            = out_folder + "sqanti/" + data_code + "_corrected.gtf"
    output:
        gtf            = out_folder + "sqanti/" + data_code + "_filter_sqanti.cds.gtf",
        classification = out_folder + "sqanti/" + data_code + "_filter_sqanti_classification.tsv",
        fasta          = out_folder + "sqanti/" + data_code + "_filter_sqanti.fasta"
    params:
        script     = "scripts/filter_sqanti.R",
        out_prefix = out_folder + "sqanti/" + data_code
    shell:
        "conda activate bambu_env; "
        "Rscript {params.script} "
        "--sqanti {input.classification} --fasta {input.fasta} --gff {input.gff} "
        "--output {params.out_prefix}"
