# ISOQUANT PIPELINE
# Jack Humphrey 2021-2024
import pandas as pd
import os
isoquant="/sc/arion/projects/ad-omics/data/software/IsoQuant/isoquant.py"
shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoquant; ml samtools; ml bedtools;")

ref_fasta = config["ref_genome"] + ".fa"
ref_gtf = config["ref_gtf"]
metadata = config["metadata"]
data_code = config["data_code"]
out_folder = config["out_folder"] 
# path to YAML sample metadata
sample_config = config["sample_config"]

# read in metadata
meta_df = pd.read_excel(metadata)
samples = meta_df['sample']

metadata_dict = meta_df.set_index("sample").T.to_dict()

# isoquant specific params
matching_strategy = "precise"
model_strategy = "default_ccs"
isoquant_threads = 16
sqanti_threads = "8"
R_VERSION = "R/4.2.2"

isoquant_prefix = out_folder + "isoquant/" + data_code + "/" + data_code + "/" + data_code
miss_prefix = out_folder + "isoquant/filter1_missingness/" + data_code
sqanti_prefix = out_folder + "isoquant/SQANTI/" + data_code
filter_prefix = out_folder + "isoquant/filter2_sqanti/" + data_code

junctionFolder = "/sc/arion/projects/als-omics/microglia_isoseq/short_read/junctions"

rule all:
    input: 
        filter_prefix + "_filter_sqanti.cds.sorted.gtf.gz",
        miss_prefix + "_miss_counts.csv"
        #out_folder + "isoquant/" + data_code + "/combined_gene_counts.tsv"

rule run_isoquant:
    input:
        bams = expand( out_folder + "{sample}/alignment/{sample}.aligned.bam", sample = samples),
    output:
        counts = isoquant_prefix + ".transcript_model_grouped_counts.tsv",
        tpm = isoquant_prefix + ".transcript_model_grouped_tpm.tsv",
        gtf = isoquant_prefix + ".extended_annotation.gtf" 
        #results/isoquant/test_isoquant4/test_isoquant4/test_isoquant4.transcript_model_grouped_counts.tsv
    shell:
        "{isoquant} --data_type pacbio_ccs --yaml {sample_config} "
        "--reference {ref_fasta} "
        "--complete_genedb "
        "--genedb {ref_gtf} --output {out_folder}isoquant/{data_code}/ "
        "--read_group file_name "

        #"--sqanti_output  "
        "--threads {isoquant_threads} " 
        #"--matching_strategy {matching_strategy} "
        #"--model_construction_strategy {model_strategy} "

## FILTER MISSINGNESS 
# remove all transcripts with greater than X% missingness 
# currently - present in at least 2 samples 
# remove monoexonic transcripts to reduce overhead for SQANTI 
# input must be simple GFF - long GTF lines break rtracklayer 
rule filter_isoquant: 
    input: 
        counts = isoquant_prefix + ".transcript_model_grouped_counts.tsv", 
        tpm = isoquant_prefix + ".transcript_model_grouped_tpm.tsv", 
        gtf = isoquant_prefix + ".extended_annotation.gtf" 
    output: 
        counts = miss_prefix + "_miss_counts.csv", 
        tpm = miss_prefix + "_miss_tpm.csv", 
        gtf = miss_prefix + "_miss.gtf" 
    params: 
        script = "scripts/bambu_filter.R", 
        prefix = miss_prefix 
    shell: 
        "ml {R_VERSION};" 
        "Rscript {params.script} --matrix {input.counts} --gff {input.gtf} --prefix {params.prefix} --remove_monoexons --tpm" 


# run SQANTI using filtered GTF
rule SQANTI:
    input:
        gtf = miss_prefix + "_miss.gtf",
        abundance = miss_prefix + "_miss_counts.csv"
    output:
        out = sqanti_prefix + "_classification.txt",
        gtf = sqanti_prefix  + "_corrected.gtf",
        fasta = sqanti_prefix + "_corrected.fasta",
        gff = sqanti_prefix + "_corrected.gtf.cds.gff"
    params:
        sample = data_code + "_isoquant",
        outDir = out_folder + "isoquant/SQANTI/",
        nCores = sqanti_threads,
        nChunks = 8,
        software= "/sc/arion/projects/ad-omics/data/software",
        #junctions = junctionArgs,
        junctions = "\'" + junctionFolder + "/*SJ.out.tab\'",
        #genome = referenceFa + ".fa",
        intropolis = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified",
        cage = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/hg38.cage_peak_phase1and2combined_coord.bed",
        polya = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/human.polyA.list.txt",
        isoAnnotGFF = "/sc/arion/projects/ad-omics/data/references/hg38_reference/RefSeq/Homo_sapiens_GRCh38_RefSeq_78.gff3"
    shell:
        "conda activate SQANTI3.env; module purge;"
        "export PYTHONPATH=$PYTHONPATH:{params.software}/cDNA_Cupcake/sequence;"
        "export PYTHONPATH=$PYTHONPATH:{params.software}/cDNA_Cupcake/;"
        "python {params.software}/SQANTI3/sqanti3_qc.py -t {params.nCores} "
        " --dir {params.outDir} "
        " --out {params.sample} "
        " -c {params.junctions} " #optional at this point
        " --cage_peak {params.cage} --polyA_motif_list {params.polya} "
        " --gtf {input.gtf} "
        #" --isoAnnotLite --gff3 {params.isoAnnotGFF}"
        " {ref_gtf} {ref_fasta} "

## filter SQANTI
rule filter_sqanti:
    input:
        counts = miss_prefix + "_miss_counts.csv",
        tpm = miss_prefix + "_miss_tpm.csv",
        gff = sqanti_prefix + "_corrected.gtf.cds.gff",
        sqanti = sqanti_prefix + "_classification.txt",
        fasta = sqanti_prefix + "_corrected.fasta"
    output:
        counts = filter_prefix + "_filter_sqanti_counts.csv",
        tpm = filter_prefix + "_filter_sqanti_tpm.csv",
        gff = filter_prefix + "_filter_sqanti.cds.gtf",
        sqanti = filter_prefix + "_filter_sqanti_classification.tsv",
        fasta =  filter_prefix + "_filter_sqanti.fasta"
    params:
        script = "scripts/filter_sqanti.R"
    shell:
        "ml {R_VERSION}; Rscript {params.script} --input {miss_prefix} --counts {input.counts} --output {filter_prefix} --sqanti {input.sqanti} --fasta {input.fasta} --gff {input.gff}"



# sort and tabix index final GFF
rule indexGFF:
    input:
        gtf = filter_prefix + "_filter_sqanti.cds.gtf",
    output:
        gtf = filter_prefix + "_filter_sqanti.cds.sorted.gtf.gz",
        index = filter_prefix + "_filter_sqanti.cds.sorted.gtf.gz.tbi"
    params:
        gff3sort = "/sc/arion/projects/ad-omics/data/software/gff3sort/gff3sort.pl"
    shell:
        "ml tabix;"
        "{params.gff3sort} {input.gtf} | bgzip > {output.gtf};"
        "tabix {output.gtf} "


