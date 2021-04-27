import os
import pandas as pd

out_folder = config["out_folder"]
data_code = config["data_code"]

genome = config["ref_genome"] + ".fa"
mmi = config["ref_genome"] + ".mmi"
genome_code = config["genome_code"]

GTF = config["ref_gtf"]
GTF_code = config["gtf_code"]

pbmm2_threads = "8"
talon_threads = "16" # for testing
genome_code = "hg38"

sqanti_threads = "4"

metadata = config["metadata"]

prefix = os.path.join(out_folder, "run_talon", data_code, data_code )

# read in xlsx metadata to pandas
meta_df = pd.read_excel(metadata)
samples = meta_df['sample']

metadata_dict = meta_df.set_index("sample").T.to_dict()

shell.prefix('export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoseq-pipeline; ')

localrules: write_config

rule all:
    input:
        expand(out_folder + "{sample}/samtools/{sample}.flagstat.txt", sample = samples),
        prefix + ".db",
        prefix + "_talon_abundance_filtered.tsv",
        prefix + "_talon_observedOnly.gtf",
        prefix + "_SQANTI_classification.txt"

# initialise TALON db
rule create_db:
    output:
        prefix + ".db"
    shell:
        "talon_initialize_database --f {GTF} --g {genome_code} --a {GTF_code} --o {prefix} "


# PBMM2 alignment
# add MD flag to reads for TALON
rule align_flnc_bam:
    #output:
     #   bam = out_folder + "pbmm2/{sample}.aligned.md.bam"
    output:
        #tmp_bam = out_folder + "pbmm2/{sample}.aligned.temp.bam"
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    run:
        input_bam = metadata_dict[wildcards.sample]["flnc_bam_path"]

        shell("pbmm2 align --sort -j {pbmm2_threads} --sort-threads 4 -m 3G --preset=ISOSEQ \
                 --log-level INFO --unmapped {mmi} {input_bam} | \
              samtools calmd -b - {genome} > {output.bam}")

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

# label reads
rule talon_label_reads:
    input:
        out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output:
        out_folder + "{sample}/label_reads/{sample}.aligned.md.labelled.sam"
    shell:
        "conda activate talon;"
        "talon_label_reads --f {input} --g {genome} --t 8 --tmpDir={out_folder}/{wildcards.sample}/label_reads/temp/ --deleteTmp "
        "--o {out_folder}/{wildcards.sample}/label_reads/{wildcards.sample};"
        "mv {out_folder}/{wildcards.sample}/label_reads/{wildcards.sample}_labeled.sam {output} ;" 
 
# create TALON config CSV - CURRENTLY DEPRECATED
rule write_config:
    input:
        sams = expand( out_folder + "{sample}/label_reads/{sample}.aligned.md.labelled.sam", sample = samples)
    output:
        config_csv = prefix + "_talon_config.csv"
    run:
        # pandas magic -  name, sample description, platform, sam file        
        sam_paths = [os.path.abspath(i) for i in input.sams]
        config_df = meta_df
        config_df["file"] = sam_paths
        config_df = config_df[ ['sample', 'description', 'platform', 'file'] ]
        config_df.to_csv(output.config_csv, header = False, index = False, sep = ",")

# filter out mono-exonic reads and intra-priming
rule filter_reads:
    input:
       out_folder + "{sample}/label_reads/{sample}.aligned.md.labelled.sam"
    output:
        out_folder + "{sample}/filtered/{sample}.aligned.md.labelled.filtered.sam" 
    params:
        fA_threshold = 0.7 # permissive threshold, can reduce later
    shell:
        "sh scripts/filter_sam.awk {input} {output} {params.fA_threshold}"

# run TALON successively for each file
rule run_talon:
    input:
        sams = expand( out_folder + "{sample}/filtered/{sample}.aligned.md.labelled.filtered.sam", sample = samples),
        talon_db = prefix + ".db"
    output:
        #"test.txt"
        expand(out_folder + "run_talon/{sample}_QC.log", sample = samples)
    run:
        sam_paths = [os.path.abspath(i) for i in input.sams]
        config_df = meta_df
        config_df["file"] = sam_paths
        config_df = config_df[ ['sample', 'description', 'platform', 'file'] ]
        for i in range(config_df.shape[0]):
            sample_loc = config_df["sample"][i] 
            print(" * processing " + sample_loc )
            output_file = out_folder + "run_talon/" + sample_loc + "_talon_read_annot.tsv" 
            # if already created then skip for this sample
            if os.path.exists(output_file):
                pass
            else:
                # create individual config file
                config_loc = config_df.iloc[i:i+1]
                config_loc_file = prefix + "_" + sample_loc + "_config.csv"
                config_loc.to_csv(config_loc_file, header = False, index = False, sep = ",")
                shell("conda activate talon; talon --f {config_loc_file} --db {input.talon_db} \
                  --build {genome_code} --threads {talon_threads} \
                  --o {out_folder}/run_talon/{sample_loc} ")


# filter loosely
rule talon_filter:
    input:
        expand(out_folder + "run_talon/{sample}_QC.log", sample = samples),
        talon_db = prefix + ".db"
    output:
        prefix + "_whitelist.txt"
    params:
        max_frac_a = 0.5,
        min_count = 5,
        min_datasets = 5
    shell:
        "talon_filter_transcripts --db {input.talon_db} "
        " -a {GTF_code} "
        " --maxFracA={params.max_frac_a} --minCount={params.min_count} --minDatasets={params.min_datasets} "
        " --o {output}"

# get abundance counts
rule talon_abundance:
    input:
        expand(out_folder + "run_talon/{sample}_QC.log", sample = samples),
        talon_db = prefix + ".db",
        whitelist = prefix + "_whitelist.txt"
    output:
        prefix + "_talon_abundance_filtered.tsv"
    shell:
        "talon_abundance --db {input.talon_db} "
        " --whitelist {input.whitelist} "
        " -b {genome_code} -a {GTF_code} "
        " --o={prefix}"
        
# create GTF
rule create_GTF:
    input:
        expand(out_folder + "run_talon/{sample}_QC.log", sample = samples),
        talon_db = prefix + ".db",
        whitelist = prefix + "_whitelist.txt"
    output:
        gtf_full = prefix + "_talon_observedOnly.gtf",
        gtf_clean = prefix + "_talon_clean.gtf"
    shell:
        "talon_create_GTF --db={input.talon_db} "
        " --whitelist {input.whitelist} "
        " -b {genome_code} -a {GTF_code} --o={prefix} --observed;"
        " gffread -ET {output.gtf_full} > {output.gtf_clean}"

## FILTER MISSINGNESS
# remove all transcripts with greater than X% missingness
# currently - present in at least 2 samples
# remove monoexonic transcripts to reduce overhead for SQANTI
# input must be simple GFF - long GTF lines break rtracklayer
rule filter_missingness:
    input:
        counts = prefix + "_talon_abundance_filtered.tsv",
        gtf = prefix + "_talon_clean.gtf"
    output:
        counts = prefix + "_filtered_strong.csv",
        gtf = prefix + "_filtered_strong.gtf"
    params:
        script = "scripts/filter_missingness.R",
        input_prefix = prefix,
        output_prefix = prefix,
        min_samples = 5, # eventually put in config
        min_reads = 5
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} --counts {input.counts} --gff {input.gtf} --prefix {params.output_prefix} --min_samples {params.min_samples} --min_reads {params.min_reads} --remove_monoexons"

    
# run SQANTI using filtered GTF
rule SQANTI:
    input:
         gtf = prefix + "_filtered_strong.gtf",
         abundance = prefix + "_filtered_strong.csv"
    output:
         out = prefix + "_SQANTI_classification.txt",
         gtf = prefix  + "_SQANTI_corrected.gtf",
         fasta = prefix + "_SQANTI_corrected.fasta" 
    params:
        sample = data_code + "_SQANTI" ,
        outDir = out_folder + "run_talon/" + data_code + "/",
        nCores = sqanti_threads,
        nChunks = 12,
        software= "/sc/arion/projects/ad-omics/data/software",
        #junctions = junctionArgs,
        #junctions = "\'" + junctionFolder + "/*SJ.out.tab\'" ,
        gtf = GTF,
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
        #" -c {params.junctions} "
        " --cage_peak {params.cage} --polyA_motif_list {params.polya} "
        "--skipORF " # ORF finding is slow, can skip if testing
        #"-c {params.intropolis}"
        " --fl_count {input.abundance}"
        " --gtf {input.gtf} "
        #" --isoAnnotLite --gff3 {params.isoAnnotGFF}"
        " {params.gtf} {genome} "



# filter SQANTI output


# run CPAT to predict ORFs
rule CPAT:
