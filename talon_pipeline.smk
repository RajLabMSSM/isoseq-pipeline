import os

out_folder = "TALON/"
data_code = "test"
talon_config = "test_talon_config.csv"
GTF = "/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.annotation.gtf"
GTF_code = "GENCODE_v30"
talon_threads = "1"
genome_string = "hg38"

prefix = os.path.join(out_folder, data_code )

rule all:
    input:
        prefix + ".db",
        prefix + "_talon_abundance.tsv",
        prefix + "_talon_observedOnly.gtf"



rule create_db:
    output:
        prefix + ".db"
    shell:
        "talon_initialize_database --f {GTF} --g {genome_string} --a {GTF_code} --o {prefix} "


# input Isoseq metadata

# PBMM2 alignment

# add MD flag to BAM files

# create TALON config CSV

rule run_talon:
    input:
        meta = talon_config,
        talon_db = prefix + ".db"
    params:
        threads = talon_threads
    output:
       prefix + "_read_annot.tsv" 
    shell:
        "ml bedtools;"
        "talon --f {input.meta} --db {input.talon_db} --build {genome_string} --threads {params.threads} "
        " --o {prefix} "

rule talon_abundance:
    input:
        annot = prefix + "_read_annot.tsv",
        talon_db = prefix + ".db"
    output:
        prefix + "_talon_abundance.tsv"
    shell:
        "talon_abundance --db {input.talon_db} -b {genome_string} -a {GTF_code} --o={prefix}"
        


rule create_GTF:
    input:
        annot = prefix + "_read_annot.tsv",
        talon_db = prefix + ".db"
    output:
        prefix + "_talon_observedOnly.gtf"
    shell:
        "talon_create_GTF --db={input.talon_db} -b {genome_string} -a {GTF_code} --o={prefix} --observed"
    
        
