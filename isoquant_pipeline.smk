import pandas as pd
import os

shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoquant; ml samtools; ml bedtools;")

reference_fasta = config["ref_genome"] + ".fa"
ref_gtf = config["ref_gtf"]
metadata = config["metadata"]
data_code = config["data_code"]
out_folder = config["out_folder"] 

# read in metadata
meta_df = pd.read_excel(metadata)
samples = meta_df['sample']

metadata_dict = meta_df.set_index("sample").T.to_dict()
#isoquant = "/sc/arion/projects/ad-omics/data/software/IsoQuant/isoquant.py"

# isoquant specific params
matching_strategy = "precise"
model_strategy = "default_ccs"
isoquant_threads = 1
rule all:
    input: 
        out_folder + "isoquant/combined_gene_counts.tsv"

# full pipeline would have to align BAM files first
rule align_flnc_bam:
    output:
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    run:
        input_bam = metadata_dict[wildcards.sample]["flnc_bam_path"]

        shell("pbmm2 align --sort -j {pbmm2_threads} --sort-threads 4 -m 3G --preset=ISOSEQ \
                 --log-level INFO --unmapped {mmi} {input_bam} | \
              samtools calmd -b - {genome} > {output.bam}")

# index bams
rule index_flnc_bam:
    input:
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output:
        bai = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam.bai"
    shell:
        "samtools index {input.bam}"

# get bam statistics
rule samtools:
    input:
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output:
        flagstat =  out_folder + "{sample}/samtools/{sample}.flagstat.txt",
        idxstat =  out_folder + "{sample}/samtools/{sample}.idxstat.txt"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat};"
        "samtools idxstats {input.bam} > {output.idxstat} "

# write list of BAM files for isoquant
# currently insists that each line is a sample but each sample should be separated by a blank line:
rule write_config:
    input:
        bams = expand( out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam", sample = samples)
    output:
        config_csv = out_folder + "isoquant/" + data_code + "_isoquant_config.csv"
    run:
        import csv
        # pandas magic -  name, sample description, platform, sam file        
        bam_paths = [os.path.abspath(i) for i in input.bams]
        with open(output.config_csv, 'w') as filehandle:
            for listitem in bam_paths:
                filehandle.write('%s\n\n' % listitem)

# create gene db from GTF
rule create_db:
    input:
        gtf = ref_gtf
    output:
        db = os.path.splitext(ref_gtf)[0] + ".db"
    shell:
        "{isoquant} --genedb {input.gtf} --complete_genedb "

# run isoquant
rule run_isoquant:
    input:
        db = "/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.annotation.db",
        bams = expand( out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam", sample = samples),
        bais = expand( out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam.bai", sample = samples),
        bam_list = out_folder + "isoquant/" + data_code + "_isoquant_config.csv"
    params:
        labels = expand( "{sample}", sample = samples )
    output:
        out_folder + "isoquant/combined_gene_counts.tsv"
    shell:
        "isoquant.py -d pacbio_ccs --polya_trimmed --fl_data --bam_list {input.bam_list} "
        "--reference {reference_fasta} --genedb {input.db} --output {out_folder}isoquant/ "
        "--labels {params.labels} " 
        "--sqanti_output --check_canonical --threads {isoquant_threads} " 
        "--matching_strategy {matching_strategy} "
        "--model_construction_strategy {model_strategy} "
