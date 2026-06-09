## PHASING LONG READ RNA-SEQ
## Jack Humphrey

## Use longcallR to phase BAMs 

## output list of phased SNPs and indels (VCF)
## output phased BAM file
## split BAM by haplotype
import pandas as pd


longcallD = "/sc/arion/projects/ad-omics/data/software/longcallD-v0.0.4_x64-linux/longcallD"
stringtie = "/sc/arion/projects/ad-omics/data/software/stringtie-3.0.0/stringtie"
#stringtie = "/sc/arion/projects/ad-omics/data/software/stringtie-2.2.1.compiled/stringtie"
prep = config["prep"]

if prep == "nanopore_cdna":
    phase_mode = "--ont"
else:
    phase_mode = ""

out_folder = config["out_folder"]
genome = config["ref_genome"] + ".fa"
metadata = config["metadata"]
meta_df = pd.read_excel(metadata)
samples = meta_df['sample']
print(samples)
metadata_dict = meta_df.set_index("sample").T.to_dict()

conditions = meta_df['condition'].unique()
print(meta_df)
isoform_gtf = config["isoform_gtf"]

haplotypes = ["h1", "h2", "h3"]

rule all:
    input:
        out_folder + "phased/stringtie_filter_fpkm.csv"
        #expand(out_folder + "phased/{condition}.phased.{haplotype}.bam", condition = conditions, haplotype = haplotypes),
        #expand(out_folder + "phased/stringtie/{condition}.phased.{haplotype}/t_data.ctab", condition = conditions, haplotype = haplotypes)
        #expand(out_folder + "phasing/{condition}.{haplotype}.bam", condition= conditions, haplotype = haplotypes),
        #expand(out_folder + "{sample}/phasing/{sample}.{haplotype}.bam", sample = samples, haplotype = haplotypes ),
        #expand(out_folder + "merged/{condition}.merged.bam", condition= conditions)

rule merge_total_bams:
    input:
        expand(out_folder + "{sample}/alignment/{sample}.aligned.bam", sample = samples )
    output:
        expand(out_folder + "merged/{condition}.merged.bam", condition= conditions),
        expand(out_folder + "merged/{condition}.merged.bam.bai", condition= conditions )
    params:
        threads = 8
    run:
        for i in conditions:
            samples_loc = meta_df[meta_df['condition'] == i]['sample']
            files_loc = [out_folder + i + "/alignment/" + i + ".aligned.bam" for i in samples_loc]
            out_file = out_folder + "merged/" + i + ".merged.bam"
            shell("ml samtools; samtools merge --threads {params.threads} -o " + out_file + " " + " ".join(files_loc))
            shell("ml samtools; samtools index {out_file}")

rule phase_bam_merged:
    input: 
        bam = out_folder + "merged/{condition}.merged.bam"
    output:
        bam = out_folder + "phased/{condition}.phased.bam",
        bai = out_folder + "phased/{condition}.phased.bam.bai",
        vcf = out_folder + "phased/{condition}.phased.vcf",
        h1 = out_folder + "phased/{condition}.phased.h1.bam",
        h2 = out_folder + "phased/{condition}.phased.h2.bam", 
        h3 = out_folder + "phased/{condition}.phased.h3.bam"
    params:
        threads = 8
    shell:
        "ml samtools;"
        "{longcallD} call -t{params.threads} {genome} {input.bam} {phase_mode} -b {output.bam} > {output.vcf};"
        "samtools index {output.bam};"
        "samtools view -bh -d HP:1 {output.bam} > {output.h1};"
        "samtools view -bh -d HP:2 {output.bam} > {output.h2};"
        "samtools view -bh -e '![HD]' {output.bam} > {output.h3};"
        "samtools index {output.h1};"
        "samtools index {output.h2};"
        "samtools index {output.h3}"

rule quant_stringtie:
    input:
        gtf = isoform_gtf,
        bam = out_folder + "phased/{condition}.phased.{haplotype}.bam"
    output:
        quant = out_folder + "phased/stringtie/{condition}.phased.{haplotype}/t_data.ctab",
        gene = out_folder + "phased/stringtie/{condition}.phased.{haplotype}/gene_quant.txt"
    shell:
        "{stringtie} -o {out_folder}/phased/stringtie/{wildcards.condition}.phased.{wildcards.haplotype}/quant.txt -A {out_folder}/phased/stringtie/{wildcards.condition}.phased.{wildcards.haplotype}/gene_quant.txt -eB -G {input.gtf} {input.bam}"

rule stringtie_filter:
    input:
        counts = expand(out_folder + "phased/stringtie/{condition}.phased.{haplotype}/t_data.ctab", condition = conditions, haplotype = haplotypes),
        gtf = isoform_gtf
    output:
        counts = out_folder + "phased/stringtie_filter_fpkm.csv",
        gtf = out_folder + "phased/stringtie_filter.gtf"
    params:
        run_code = "phased",
        script = "scripts/stringtie_filter.R",
        prefix = out_folder + "phased/stringtie",
        min_samples = 1, # eventually put in config
        min_reads = 0 # FPKM
    shell:
        "ml R; Rscript {params.script} --inFolder {out_folder}/phased/ --runCode {params.run_code} --gff {input.gtf} --prefix {params.prefix} --min_samples {params.min_samples} --remove_monoexons"
