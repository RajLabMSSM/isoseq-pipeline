
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
pbmm2_threads = "8"
# file either TSV or XLSX
if ".tsv" in metadata:
    meta_df = pd.read_csv(metadata, sep = '\t')
if ".xlsx" in metadata:
    meta_df = pd.read_excel(metadata)

samples = meta_df['sample']
#sampleFA = meta_df['fasta_path']
#sampleREP = meta_df['cluster_report_path']

print(samples)

metadata_dict = meta_df.set_index("sample").T.to_dict()
#meta_df = meta_df.set_index("sample")

#print(samples)
#rawFolder = config['rawFolder']
out_folder = config['out_folder']
data_code = config['data_code']
ref_genome = config['ref_genome']


mmi = ref_genome + ".mmi"
genome = ref_genome + ".fa"

## barcode and UMI design
barcode_design = config['barcode_design']

mmi = ref_genome + ".mmi"
genome = ref_genome + ".fa"

ref_gtf = config['ref_gtf']

cupcake_cores = 16
# short read junctions
#junctionFolder = "/sc/arion/projects/als-omics/microglia_isoseq/short_read/junctions"
#junctions =  glob.glob(junctionFolder + "/*SJ.out.tab")
#junctionArgs = " ".join(["-c " + f for f in junctions])

#chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "M"]

rule all:
    input:
        expand(out_folder + "pbmm2/" + "{sample}.bam", sample = samples),
        expand(out_folder + "dedup/" + "{sample}.dedup.info.csv", sample = samples),
        expand(out_folder + "lima/" + "{sample}.aligned.bam", sample = samples),
        expand(out_folder + "tag/" + "{sample}.barcodes.txt", sample = samples)

        #out_folder + "cupcake/" + data_code + ".cupcake.collapsed.gff"
        #out_folder + "flnc_bam/all_samples.flnc.aligned.bam",
      #out_folder + "cupcake/" + data_code + ".demux_fl_count.csv"
      #out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt",
      #out_folder + "SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz", 
      #out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",

print(metadata_dict)

# remove primers
rule lima:
    input:
        primers = "ref/primers.fasta"
    output:
        out_folder + "lima/" + "{sample}.5p--3p.bam"
    params:
        prefix = out_folder + "lima/" + "{sample}.bam"
    run:
        ccs = metadata_dict[wildcards.sample]["ccs_bam_path"]
        shell("lima --isoseq --dump-clips {ccs} {input.primers} {params.prefix}")

#
rule detect_barcodes_umi:
    input:
        out_folder + "lima/" + "{sample}.5p--3p.bam"
    output:
        out_folder + "tag/" + "{sample}.bam"
    shell:
        "isoseq3 tag {input} {output} --design {barcode_design}"# T-8U-12B"

rule write_barcodes:
    input:
        out_folder + "lima/" + "{sample}.5p--3p.bam"
    output:
        out_folder + "tag/" + "{sample}.barcodes.txt"
    shell:
        "samtools view {input} | awk '{{x=$(NF - 2);gsub(\"XC:Z:\",\"\",x); print x }}' > {output}

rule refine:
    input:
        flt = out_folder + "tag/" + "{sample}.bam",
        primers = "ref/primers.fasta"
    output:
        out_folder + "refine/" + "{sample}.bam"
    shell:
        "isoseq3 refine {input.flt} {input.primers} {output} --require-polya"


rule dedup:
    input:
        out_folder + "refine/" + "{sample}.bam"
    output:
        bam = out_folder + "dedup/" + "{sample}.bam",
        fasta = out_folder + "dedup/" + "{sample}.fasta"
    shell:
        "isoseq3 dedup {input} {output.bam} --verbose"

rule make_dedup_csv:
    input:
        out_folder + "dedup/" + "{sample}.fasta"
    output:
        out_folder + "dedup/" + "{sample}.dedup.info.csv"
    params:
        script = "/sc/arion/projects/ad-omics/data/software/cDNA_Cupcake/singlecell/make_csv_for_dedup.py"
    shell:
        "cd {out_folder}/dedup/;"
        "cp {wildcards.sample}.fasta dedup.fasta;" 
        "python {params.script} ;"
        "mv dedup.info.csv {wildcards.sample}.dedup.info.csv"

rule align_lima:
    input:
        bam = out_folder + "lima/" + "{sample}.5p--3p.bam"
    output:
        bam = out_folder + "lima/" + "{sample}.aligned.bam"
    run:
        shell("pbmm2 align --sort -j {pbmm2_threads} --sort-threads 4 -m 3G --preset=ISOSEQ \
        --log-level INFO --unmapped {mmi} {input.bam} | \
        samtools calmd -b - {genome} > {output.bam}")

rule align:
    input:
        bam = out_folder + "dedup/" + "{sample}.bam"
    output:
        bam = out_folder + "pbmm2/" + "{sample}.bam"
    run:
        shell("pbmm2 align --sort -j {pbmm2_threads} --sort-threads 4 -m 3G --preset=ISOSEQ \
        --log-level INFO --unmapped {mmi} {input.bam} | \
        samtools calmd -b - {genome} > {output.bam}")
        shell("samtools index {output.bam}")

