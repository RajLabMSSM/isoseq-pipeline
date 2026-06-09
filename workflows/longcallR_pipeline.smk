# LongCallR - phase long-read RNA-seq for allele-specific tests

# Jack Humphrey 2026
import pandas as pd
import os

lcr_threads = 16

shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoseq-pipeline;")

longcallr = "/sc/arion/projects/ad-omics/data/software/longcallR/target/release/longcallR"
asj_to_bed = "/sc/arion/projects/ad-omics/data/software/longcallR/allele_specific/asj_to_bed.py"

rediportal = "/sc/arion/projects/als-omics/microglia_isoseq/isoseq-pipeline/lcr_test/TABLE1_hg38_v3.txt.gz"
region_string = ""

ref_genome = config["ref_genome"] + ".fa"
ref_gtf = config["ref_gtf"]
metadata = config["metadata"]
data_code = config["data_code"]
outFolder = config["out_folder"]
prep = config["prep"]

gwas = "/sc/arion/projects/ad-omics/data/references/GWAS/Bellenguez_AD/Bellenguez_2021.processed.tsv.gz"

# read in metadata
meta_df = pd.read_excel(metadata)
samples = meta_df['sample']

# for testing
#samples = "16-078_2_MFG"
#region_string = "-r chr16:74821372-89210978"
#outFolder = "test_lcr_editing"

metadata_dict = meta_df.set_index("sample").T.to_dict()

if prep == "pacbio":
    input_bam = "results/{sample}/alignment" + "/{sample}.aligned.bam"
    longcallr_string = "hifi-masseq"

if prep == "nanopore_direct":
    input_bam = "results/{sample}/alignment" + "/{sample}.aligned.bam"
    minimap_string = "-ax splice -uf -k14"
    longcallr_string = "ont-drna"

if prep == "nanopore_cdna":
    input_bam = "results/{sample}/alignment" + "/{sample}.aligned.bam"
    minimap_string = "-ax splice"
    longcallr_string = "ont-cdna"



rule all:
    input:
        expand(outFolder + "/{sample}/phased/{sample}.{test}.tsv", sample = samples, test = ["ase","asj","asediting"] ),
        outFolder + "/GWAS/Bellenguez.ase_gwas_joint.tsv",
        outFolder + "/GWAS/Bellenguez.asj_gwas_joint.tsv",
        outFolder + "/GWAS/Bellenguez.asediting_gwas_joint.tsv",
        outFolder + "longcallR/" + data_code + ".phasing_rates.tsv"

rule longcallR:
    input:
        bam = input_bam
    output: "{outFolder}/{sample}/phased/{sample}.phased.bam"
    params:
        prefix = "{outFolder}/{sample}/phased/{sample}"    
    run:
        shell("{longcallr} -b {input.bam} -f {ref_genome} -p {longcallr_string} -t {lcr_threads} -o {params.prefix} {region_string}")
        shell("ml samtools; samtools index {output}")

# for each phased bam output table on how many reads were able to be phased
rule get_phasing_rate:
    input:
        expand(outFolder + "/{sample}/phased/{sample}.phased.bam", sample = samples)
    output:
        outFolder + "longcallR/" + data_code + ".phasing_rates.tsv"
    params:
        script = "scripts/get_phasing_rate.py"
    run:
        shell("ml samtools; python {params.script} --results-dir {outFolder} --out {output} -t {lcr_threads}")
rule split_bam:
    input:
        bam =  "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        h1 = "{outFolder}/{sample}/phased/{sample}.phased.hap1.bam",
        h2 = "{outFolder}/{sample}/phased/{sample}.phased.hap2.bam"
    run:
        shell("ml samtools;\
            samtools view -h -b -d HP:1 {input.bam} > {output.h1}\
            samtools index {output.h1}\
            samtools view -h -b -d HP:2 {input.bam} > {output.h2}\
            samtools index {output.h2}")
"{sample}/phased/{sample}.ase.tsv"
rule allele_specific_splicing:
    input:
        bam = "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        bed = "{outFolder}/{sample}/phased/{sample}.asj.0.05.bed",
        tsv = "{outFolder}/{sample}/phased/{sample}.asj.tsv"
    params:
        prefix = "{outFolder}/{sample}/phased/{sample}"
    run:
        shell("{longcallr} asj -a {ref_gtf} -b {input.bam} -f {ref_genome} -o {params.prefix} -t {lcr_threads}")
        shell("{asj_to_bed} {output.tsv} 0.05 > {output.bed}")

rule allele_specific_expression:
    input:
        bam = "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        ase = "{outFolder}/{sample}/phased/{sample}.ase.tsv"
    params:
        prefix = "{outFolder}/{sample}/phased/{sample}"
    run:
        shell("{longcallr} ase -a {ref_gtf} -b {input.bam} -o {params.prefix} -t {lcr_threads}")

rule allele_specific_editing:
    input:
        bam = "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        ased = "{outFolder}/{sample}/phased/{sample}.asediting.tsv"
    params:
        script = "scripts/longcallR-asediting_v13.py",
        prefix = "{outFolder}/{sample}/phased/{sample}"
    run:
        shell("python {params.script} \
        -b {input.bam} \
        -r {rediportal} \
        -o {params.prefix} \
        -a {ref_gtf} \
        -t {lcr_threads} \
        --min_coverage 10 \
        --min_editing_rate 0.01")

#rule assemble_{outFolder}:
# in R
#x <- list.files(pattern = "TSPAN14.*asj.tsv", recursive = TRUE); names(x) <- dirname(x)
#d <- map_df(x, read_tsv, .id = "sample", col_types = "cccnnnnnnnllc") %>% arrange(P_value)

#x <- list.files(pattern = "TSPAN14.*asediting.tsv", recursive = TRUE); names(x) <- dirname(x)
#d2 <- map_df(x, read_tsv, .id = "sample", col_types = "cnccnnnnnnnnnnnl") %>% arrange(P_value)

# GWAS testing
# for a GWAS, extract genome-wide significant SNPs
# overlap with SNPs used to phase RNA
# for each ASE, ASJ, ASED, find features (genes, junctions, edSites) within those haplotypes
# get all AS features the right way around, with respect to the GWAS SNP
# compare allelic FCs and do binomial test to get P-value of consistent direct of allelic bias
# do inverse-variance weighting meta-analysis when present in 2 or more samples
# scripts]$ python longcallR-gwas-joint.py -h
#
rule integrate_gwas:
    input:
        ase = expand(outFolder + "/{sample}/phased/{sample}.ase.tsv", sample = samples),
        asj = expand(outFolder + "/{sample}/phased/{sample}.asj.tsv", sample = samples),
        ased = expand(outFolder + "/{sample}/phased/{sample}.asediting.tsv", sample = samples)
    output:
        ase = outFolder + "/GWAS/Bellenguez.ase_gwas_joint.tsv",
        asj = outFolder + "/GWAS/Bellenguez.asj_gwas_joint.tsv",
        ased = outFolder + "/GWAS/Bellenguez.asediting_gwas_joint.tsv"
    params:
        script = "scripts/longcallR-gwas-joint_v3.py"
    run:
        shell("python {params.script} \
            --results-dir {outFolder} \
            --mode ase \
            --gwas {gwas} \
            --min-samples 1 \
            --out {outFolder}/GWAS/Bellenguez")
        shell("python {params.script} \
            --results-dir {outFolder} \
            --mode asj \
            --min-samples 1 \
            --gwas {gwas} \
            --out {outFolder}/GWAS/Bellenguez")
        shell("python {params.script} \
            --results-dir {outFolder} \
            --mode asediting \
            --min-samples 1 \
            --gwas {gwas} \
            --out {outFolder}/GWAS/Bellenguez")


