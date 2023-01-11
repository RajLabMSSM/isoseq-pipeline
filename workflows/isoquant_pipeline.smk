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
isoquant_threads = 16

rule all:
    input: 
        out_folder + "isoquant/" + data_code + "/combined_gene_counts.tsv"

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
                filehandle.write('%s\n' % listitem)

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
        out_folder + "isoquant/" + data_code + "/combined_gene_counts.tsv"
    shell:
        "isoquant.py --data_type pacbio_ccs --bam_list {input.bam_list} "
        "--read_group file_name "
        "--reference {reference_fasta} "
        "--genedb {input.db} --output {out_folder}isoquant/{data_code}/ "
        #"--labels {params.labels} " 
        #"--sqanti_output  "
        "--threads {isoquant_threads} " 
        #"--matching_strategy {matching_strategy} "
        #"--model_construction_strategy {model_strategy} "
