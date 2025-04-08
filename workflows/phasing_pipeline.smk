## PHASING LONG READ RNA-SEQ
## Jack Humphrey

## Use longcallD to phase BAMs 

## output list of phased SNPs and indels (VCF)
## output phased BAM file
## split BAM by haplotype
import pandas as pd

longcallD = "/sc/arion/projects/ad-omics/data/software/longcallD-v0.0.4_x64-linux/longcallD"
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

haplotypes = ["h1", "h2", "h3"]
rule all:
    input:
        expand(out_folder + "phasing/{condition}.{haplotype}.bam", condition= conditions, haplotype = haplotypes),
        expand(out_folder + "{sample}/phasing/{sample}.{haplotype}.bam", sample = samples, haplotype = haplotypes ),
        expand(out_folder + "merged/{condition}.merged.bam", condition= conditions)

rule phase_bam:
    input: 
        bam = out_folder + "{sample}/alignment/{sample}.aligned.bam"
    output:
        bam = out_folder + "{sample}/phasing/{sample}.phased.bam",
        bai = out_folder + "{sample}/phasing/{sample}.phased.bam.bai",
        vcf = out_folder + "{sample}/phasing/{sample}.phased.vcf"
    params:
        threads = 8
    shell:
        "ml samtools;"
        "{longcallD} call -t{params.threads} {genome} {input.bam} {phase_mode} -b {output.bam} > {output.vcf};"
        "samtools index {output.bam}"

# h1 - haplotype 1
# h2 - haplotype 2
# h3 - no haplotype assigned
rule split_haplotype:
    input:
        bam = out_folder + "{sample}/phasing/{sample}.phased.bam",
    output:
        h1 = out_folder + "{sample}/phasing/{sample}.h1.bam",
        h2 = out_folder + "{sample}/phasing/{sample}.h2.bam",
        h3 = out_folder + "{sample}/phasing/{sample}.h3.bam",
    shell:
        "samtools view -bh -d HP:1 {input.bam} > {output.h1};"
        "samtools view -bh -d HP:2 {input.bam} > {output.h2};"
        "samtools view -bh -e '![HD]' {input.bam} > {output.h3};"
        "samtools index {output.h1};"
        "samtools index {output.h2};"
        "samtools index {output.h3}"

rule merge_total_bams:
    input:
        expand(out_folder + "{sample}/phasing/{sample}.{haplotype}.bam", sample = samples, haplotype = haplotypes )
    output:
        expand(out_folder + "merged/{condition}.merged.bam", condition= conditions),
        expand(out_folder + "merged/{condition}.merged.bam.bai", condition= conditions )
    run:
        for i in conditions:
            samples_loc = meta_df[meta_df['condition'] == i]['sample']    
            files_loc = [out_folder + i + "/alignment/" + i + ".aligned.bam" for i in samples_loc]
            out_file = out_folder + "merged/" + i + ".merged.bam"
            shell("ml samtools; samtools merge -o " + out_file + " " + " ".join(files_loc))
            shell("ml samtools; samtools index {out_file}")

rule merge_phased_bams:
    input: 
        expand(out_folder + "{sample}/phasing/{sample}.{haplotype}.bam", sample = samples, haplotype = haplotypes )
    output:
        expand(out_folder + "phasing/{condition}.{haplotype}.bam", condition= conditions, haplotype = haplotypes ),
        expand(out_folder + "phasing/{condition}.{haplotype}.bam.bai", condition= conditions, haplotype = haplotypes )
    run:
        for i in conditions:
            samples_loc = meta_df[meta_df['condition'] == i]['sample']    
            for hap in haplotypes:
                files_loc = [out_folder + i + "/phasing/" + i + "." + hap + ".bam" for i in samples_loc]
                out_file = out_folder + "phasing/" + i + "." + hap + ".bam"
                shell("ml samtools; samtools merge -o " + out_file + " " + " ".join(files_loc))  
                shell("ml samtools; samtools index {out_file}")
