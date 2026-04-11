# LongCallR - phase long-read RNA-seq for allele-specific tests

# Jack Humphrey 2026
import pandas as pd
import os

shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoquant; ml samtools; ml bedtools;")

longcallr="/sc/arion/projects/ad-omics/data/software/longcallR/target/release/longcallR"
asj_to_bed="/sc/arion/projects/ad-omics/data/software/longcallR/allele_specific/asj_to_bed.py"
#genome=/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v38.primary_assembly/GRCh38.primary_assembly.genome.fa
#gtf=/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v38.primary_assembly/gencode.v38.primary_assembly.annotation.gtf
rediportal="/sc/arion/projects/als-omics/microglia_isoseq/isoseq-pipeline/lcr_test/TABLE1_hg38_v3.txt.gz"

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

isoquant_prefix = out_folder + "isoquant/" + data_code + "/" + data_code + "/" + data_code
miss_prefix = out_folder + "isoquant/filter1_missingness/" + data_code
sqanti_prefix = out_folder + "isoquant/SQANTI/" + data_code + "_isoquant"
filter_prefix = out_folder + "isoquant/filter2_sqanti/" + data_code

junctionFolder = config["junctionFolder"]

# check all files in sample config exist
import yaml
import io
with open(sample_config, 'r') as stream:
    data_loaded = yaml.safe_load(stream)

short_read = data_loaded[1]['illumina bam']
long_read = data_loaded[1]['long read files']
all_files = short_read + long_read
file_check = [os.path.isfile(f) for f in all_files]

if not all(file_check):
    print(" * missing files in sample config!")
    missing_files = [all_files[i] for i in [j for j in range(len(file_check)) if not file_check[j]]]
    print(missing_files)
    exit(1)

rule all:
    input:

rule longcallR:
    output: prefix + prefix + ".phased.bam"
    run:

        shell("longcallr -b $bam -f $genome -p hifi-masseq -t 8 -o ${prefix}/${prefix}")
        shell("ml samtools; samtools index ${prefix}/${prefix}.phased.bam")

rule split_bam:
    input:
        prefix + prefix + ".phased.bam"
    output:
        h1 = prefix + prefix + ".hap1.bam"
        h2 = prefix + prefix + ".hap2.bam"
    run:
        shell("ml samtools;\
            samtools view -h -b -d HP:1 ${prefix}/${prefix}.phased.bam > ${prefix}/${prefix}.hap1.bam\
            samtools index ${prefix}/${prefix}.hap1.bam\
            samtools view -h -b -d HP:2 ${prefix}/${prefix}.phased.bam > ${prefix}/${prefix}.hap2.bam\
            samtools index ${prefix}/${prefix}.hap2.bam")

rule allele_specific_splicing:
$longcallr asj -a $gtf -b ${prefix}/${prefix}.phased.bam -f $genome -o ${prefix}/${prefix} -t 8
$asj_to_bed ${prefix}/${prefix}.asj.tsv 0.05 > ${prefix}/${prefix}.asj.005.bed


rule allele_specific_expression:
$longcallr ase -a $gtf -b ${prefix}/${prefix}.phased.bam -o ${prefix}/${prefix} -t 8

rule allele_specific_editing:
    python longcallR-asediting_v6.py \
    -b ${prefix}/${prefix}.phased.bam \
    -r $rediportal \
    -f $genome \
    -o ${prefix}/${prefix} \
    -a ${gtf} \
    -t 8 \
    --min_coverage 10 \
    --min_editing_rate 0.01

rule assemble_results:
# in R
x <- list.files(pattern = "TSPAN14.*asj.tsv", recursive = TRUE); names(x) <- dirname(x)
d <- map_df(x, read_tsv, .id = "sample", col_types = "cccnnnnnnnllc") %>% arrange(P_value)

x <- list.files(pattern = "TSPAN14.*asediting.tsv", recursive = TRUE); names(x) <- dirname(x)
d2 <- map_df(x, read_tsv, .id = "sample", col_types = "cnccnnnnnnnnnnnl") %>% arrange(P_value)

longcallr=/sc/arion/projects/ad-omics/data/software/longcallR/target/release/longcallR
asj_to_bed=/sc/arion/projects/ad-omics/data/software/longcallR/allele_specific/asj_to_bed.py
genome=/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v38.primary_assembly/GRCh38.primary_assembly.genome.fa
gtf=/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v38.primary_assembly/gencode.v38.primary_assembly.annotation.gtf
#bam=/sc/arion/projects/als-omics/microglia_isoseq/isoseq-pipeline/results/MG-18_1_MFG/alignment/MG-18_1_MFG.aligned.bam
#region=chr19:51225012-51226118
region=chr19:51187019-51283924
region=chr19:48793064-48813029
region=chr6:41147573-41164994
region=chr16:81714393-81989729
region=chr10:80438784-80556752
rediportal=/sc/arion/projects/als-omics/microglia_isoseq/isoseq-pipeline/lcr_test/TABLE1_hg38_v3.txt.gz

bam_list=all_bams.txt

for bam in $(cat $bam_list );do
    prefix=$(basename $bam)
    prefix=${prefix%.aligned.bam}.TSPAN14
    echo ${prefix} 
    mkdir -p ${prefix}
    # phase reads
    $longcallr -b $bam -f $genome -p hifi-masseq -t 8 -o ${prefix}/${prefix} -r $region
    samtools index ${prefix}/${prefix}.phased.bam
    
    # split BAM into two haplotypes
    # Haplotype 1
    samtools view -h -b -d HP:1 ${prefix}/${prefix}.phased.bam > ${prefix}/${prefix}.hap1.bam
    samtools index ${prefix}/${prefix}.hap1.bam

    # Haplotype 2
    samtools view -h -b -d HP:2 ${prefix}/${prefix}.phased.bam > ${prefix}/${prefix}.hap2.bam
    samtools index ${prefix}/${prefix}.hap2.bam
    
    # allele-specific tests
    $longcallr asj -a $gtf -b ${prefix}/${prefix}.phased.bam -f $genome -o ${prefix}/${prefix} -t 8
    $longcallr ase -a $gtf -b ${prefix}/${prefix}.phased.bam -o ${prefix}/${prefix} -t 8
    # make BED file for junctions
    $asj_to_bed ${prefix}/${prefix}.asj.tsv 0.05 > ${prefix}/${prefix}.asj.005.bed

    # experimental - allele-specific editing
    python longcallR-asediting_v6.py \
    -b ${prefix}/${prefix}.phased.bam \
    -r $rediportal \
    -f $genome \
    -o ${prefix}/${prefix} \
    -a ${gtf} \
    -t 8 \
    --min_coverage 10 \
    --min_editing_rate 0.01

done
