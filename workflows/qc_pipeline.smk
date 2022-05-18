
# merge FLNC BAMs and reports
# collapse with Cupcake
# demultiplex to get transcript counts

#cupcake_path = "/sc/arion/projects/ad-omics/data/software/cDNA_Cupcake"
rrna = "/sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/gencode.v30.rRNA.interval.list"
shell.prefix('export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoseq-pipeline; ')
#shell.prefix('export PS1="";source activate isoseq-pipeline;ml R/3.6.0;')

import pandas as pd
import xlrd
import glob

metadata = config['metadata']

reflat_file =  "/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.primary_assembly.annotation.reflat"
genes_file = "/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.annotation.gtf.genes"

# file either TSV or XLSX
if ".tsv" in metadata:
    meta_df = pd.read_csv(metadata, sep = '\t')
if ".xlsx" in metadata:
    meta_df = pd.read_excel(metadata)

prep = config["prep"]
print(prep)
samples = meta_df['sample']
#sampleFA = meta_df['fasta_path']
#sampleREP = meta_df['cluster_report_path']

metadata_dict = meta_df.set_index("sample").T.to_dict()

#print(samples)
#rawFolder = config['rawFolder']
out_folder = config['out_folder']
data_code = config['data_code']
ref_genome = config['ref_genome']

out_folder = os.path.join(out_folder, data_code) + "/"

mmi = ref_genome + ".mmi"
genome = ref_genome + ".fa"

ref_gtf = config['ref_gtf']

pbmm2_threads = "8"
#chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "M"]

if prep == "pacbio":
    output_bam = out_folder + "merged_bam/" + data_code + ".flnc.aligned.md.bam"

if prep == "nanopore_direct":
    output_bam = out_folder + "{sample}/alignment" + "/{sample}.aligned.bam"
    minimap_string = "-ax splice -uf -k14" 

rule all:
    input:
        out_folder + "multiqc/multiqc_report.html"
        #expand(output_bam
        #out_folder + "merged_bam/" + data_code + ".flnc.aligned.md.bam"
        #out_folder + "cupcake/" + data_code + "/" + data_code + ".cupcake.collapsed.gff"
      #out_folder + "cupcake/" + data_code + ".cupcake.collapsed.gff"
        #out_folder + "flnc_bam/all_samples.flnc.aligned.bam",
      #out_folder + "cupcake/" + data_code + ".demux_fl_count.csv"
      #out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt",
      #out_folder + "SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz", 
      #out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",

rule create_index:
    input:
        genome
    output:
        mmi
    shell:
        "pbmm2 index {input} {output}"

# PACBIO STEPS 

# PBMM2 alignment
# # add MD flag to reads for TALON
rule align_flnc_bam:
    input:
        mmi
    output:
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    run:
        input_bam = metadata_dict[wildcards.sample]["flnc_bam_path"]
        shell("pbmm2 align --sort -j {pbmm2_threads} --sort-threads 4 -m 3G --preset=ISOSEQ \
        --log-level INFO --unmapped {mmi} {input_bam} | \
        samtools calmd -b - {genome} > {output.bam}")


# merge aligned flnc.bam files
rule merge_flnc_bams:
    input:
        expand(out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam", sample = samples)
    params:
        tmp = out_folder + "merged_bam/" + data_code + ".flnc.aligned.tmp.bam"
    output:
        out_folder + "merged_bam/" + data_code + ".flnc.aligned.md.bam"
    shell:
        "samtools merge {params.tmp} {input};"
        "samtools view -bh -F 4 {params.tmp} > {output};"
        "rm {params.tmp};"
        "samtools index {output}"

## NANOPORE SPECIFIC STEPS
rule nanopore_alignment:
     output:
         out_folder + "{sample}/alignment" + "/{sample}.aligned.sam"
     run:
        fastq_file = metadata_dict[wildcards.sample]["fastq"]

        shell( "minimap2 -ax splice -uf -k14 {genome} {fastq_file} > {output} ")

## convert to bam and coordinate sort
rule sam_to_bam:
        input:
            out_folder + "{sample}/alignment" + "/{sample}.aligned.sam"
        params:
            tmp = out_folder + "{sample}/alignment" + "/{sample}.tmp.bam"
        output:
            output_bam
            #out_folder + "{sample}/alignment" + "/{sample}.aligned.bam"
        shell:
            'ml samtools;'
            'samtools view -bh {input} > {params.tmp};'
            'samtools sort -o {output} {params.tmp};'


## ALL QC STEPS


# QC
rule rnaseqc:
    input:
        geneGTF = genes_file, #ref_gtf + ".genes",
        bam = output_bam
        #bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    params:
        out =  out_folder + "{sample}/qc/"
    output:
         out_folder + "{sample}/qc/{sample}.metrics.tsv"
    shell:
        "ml rnaseqc;"
        "rnaseqc {input.geneGTF} {input.bam} {params.out} "
        " --sample={wildcards.sample} "
        " --unpaired --coverage --verbose --mapping-quality 0 --base-mismatch=1000 --detection-threshold=1"

# FASTQC
rule fastqc:
    input: 
        bam = output_bam
        #bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output: 
        out_folder + "{sample}/qc/{sample}.aligned_fastqc.html"
    shell:
        "ml fastqc;"
        "fastqc --threads 8 --outdir={out_folder}/{wildcards.sample}/qc/ --format bam {input.bam}"

# Picard
rule picard:
    input:
        bam = output_bam,
        #bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam",
        reflat = reflat_file
    output:
        out_folder + "{sample}/qc/{sample}.RNASeqMetrics"
    shell:
        "ml picard;"
        "java -jar $PICARD CollectRnaSeqMetrics "
        "I={input.bam} O={output} "
        "REF_FLAT={reflat_file} "
        "STRAND=FIRST_READ_TRANSCRIPTION_STRAND "
        "RIBOSOMAL_INTERVALS={rrna} "
        "VALIDATION_STRINGENCY=LENIENT " 

# Samtools
rule samtools:
    input:
        output_bam
        #out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output:
        idx = out_folder + "{sample}/qc/{sample}.idxstat.txt",
        flag = out_folder +"{sample}/qc/{sample}.flagstat.txt"
        #bai = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam.bai"
    shell:
        "ml samtools;"
        "samtools index {input};"
        "samtools flagstat {input} > {output.flag};"
        "samtools idxstat {input} > {output.idx}"


# multiqc, version 1.8.dev0 works with rnaseqc outputs
rule multiQC:
    input:
        expand(out_folder + "{sample}/qc/{sample}.RNASeqMetrics", sample = samples),
        expand(out_folder + "{sample}/qc/{sample}.metrics.tsv", sample = samples),
        expand(out_folder + "{sample}/qc/{sample}.flagstat.txt", sample = samples),
        expand(out_folder + "{sample}/qc/{sample}.idxstat.txt", sample = samples),
        expand(out_folder + "{sample}/qc/{sample}.aligned_fastqc.html", sample = samples)
    output:
         out_folder + "multiqc/multiqc_report.html"
    shell:
        "export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8;"
        "conda activate snakemake;"
        "multiqc -f --outdir {out_folder}/multiqc/ {out_folder}" 

