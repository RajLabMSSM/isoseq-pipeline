
# merge FLNC BAMs and reports
# collapse with Cupcake
# demultiplex to get transcript counts

#cupcake_path = "/sc/arion/projects/ad-omics/data/software/cDNA_Cupcake"

shell.prefix('export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoseq-pipeline; ')
#shell.prefix('export PS1="";source activate isoseq-pipeline;ml R/3.6.0;')

import pandas as pd
import xlrd
import glob

metadata = config['metadata']

reflat_file =  "/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.primary_assembly.annotation.reflat"

# file either TSV or XLSX
if ".tsv" in metadata:
    metaDF = pd.read_csv(metadata, sep = '\t')
if ".xlsx" in metadata:
    metaDF = pd.read_excel(metadata)

samples = metaDF['sample']
sampleFA = metaDF['fasta_path']
sampleREP = metaDF['cluster_report_path']

metaDF = metaDF.set_index("sample")

#print(samples)
#rawFolder = config['rawFolder']
out_folder = config['out_folder']
data_code = config['data_code']
ref_genome = config['ref_genome']

mmi = ref_genome + ".mmi"
genome = ref_genome + ".fa"

ref_gtf = config['ref_gtf']

#chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "M"]

rule all:
    input:
        out_folder + "multiqc/multiqc_report.html",
        out_folder + "merged_bam/" + data_code + ".flnc.aligned.md.bam"
        #out_folder + "cupcake/" + data_code + "/" + data_code + ".cupcake.collapsed.gff"
      #out_folder + "cupcake/" + data_code + ".cupcake.collapsed.gff"
        #out_folder + "flnc_bam/all_samples.flnc.aligned.bam",
      #out_folder + "cupcake/" + data_code + ".demux_fl_count.csv"
      #out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt",
      #out_folder + "SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz", 
      #out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",

# PBMM2 alignment
# # add MD flag to reads for TALON
rule align_flnc_bam:
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

# QC
rule rnaseqc:
    input:
        geneGTF = ref_gtf + ".genes",
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
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
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output: 
        out_folder + "{sample}/qc/{sample}.aligned.md_fastqc.html"
    shell:
        "ml fastqc;"
        "fastqc --outdir={out_folder}/{wildcards.sample}/qc/ --format bam {input.bam}"

# Picard
rule picard:
    input:
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam",
        reflat = reflat_file
    output:
        out_folder + "{sample}/qc/{sample}.RNASeqMetrics"
    shell:
        "ml picard;"
        "java -jar $PICARD CollectRnaSeqMetrics "
        "I={input.bam} O={output} "
        "REF_FLAT={reflat_file} "
        "STRAND=FIRST_READ_TRANSCRIPTION_STRAND "
        "VALIDATION_STRINGENCY=LENIENT " 

# Samtools
rule samtools:
    input:
        out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output:
        idx = out_folder + "{sample}/qc/{sample}.idxstat.txt",
        flag = out_folder +"{sample}/qc/{sample}.flagstat.txt",
        bai = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam.bai"
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
        expand(out_folder + "{sample}/qc/{sample}.aligned.md_fastqc.html", sample = samples)
    output:
         out_folder + "multiqc/multiqc_report.html"
    shell:
        "export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8;"
        "multiqc -f --outdir {out_folder}/multiqc/ {out_folder}" 

