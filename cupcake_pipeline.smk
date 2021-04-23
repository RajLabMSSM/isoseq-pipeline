
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

cupcake_cores = 32
# short read junctions
#junctionFolder = "/sc/arion/projects/als-omics/microglia_isoseq/short_read/junctions"
#junctions =  glob.glob(junctionFolder + "/*SJ.out.tab")
#junctionArgs = " ".join(["-c " + f for f in junctions])

localrules: create_chain_config

#chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "M"]

rule all:
    input:
      out_folder + "cupcake/" + data_code + ".cupcake.collapsed.gff"
        #out_folder + "flnc_bam/all_samples.flnc.aligned.bam",
      #out_folder + "cupcake/" + data_code + ".demux_fl_count.csv"
      #out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt",
      #out_folder + "SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz", 
      #out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",

# run FASTQC on each BAM separately
rule fastqc:
    input: "{sample}/minimap/{sample}.hq.sam"
    output: "{sample}/qc/{sample}.hq_fastqc.html"
    shell:
        "ml fastqc;"
        "fastqc --outdir={wildcards.sample}/qc/ --format bam {input}"

# prep FLNC report files for demultiplexing later
rule prep_flnc_reports:
    input:
        metadata
    params:
        script = "scripts/prepare_flnc_reports.R"
    output:
        out_folder + "cupcake/" + data_code + ".merged.flnc_report.csv"
    shell:
        "ml R/3.6.0; Rscript {params.script} {input} {output}" 

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

## Cupcake Tools

# collapse redundant reads
# deal with 5' truncated reads - don't use dun-merge-5-shorter or filter_away
# filter away  
# CHECK - does cupcake mind if FASTA and cluster_report are gzipped?
rule cupcake_collapse:
    input:
        bam = out_folder + "merged_bam/" + data_code + ".flnc.aligned.md.bam"
    output:
         gff = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.gff",
         #fasta = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.rep.fa",
         #group = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.group.txt",
         #stat = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.read_stat.txt",
         #count = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.abundance.txt"
         #fasta2 = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.rep.fa",
         #chain_gff = out_folder + "cupcake/chain/cupcake.collapsed.gff",
         #chain_group = out_folder + "cupcake/chain/cupcake.collapsed.group.txt",
         #chain_count = out_folder + "cupcake/chain/cupcake.collapsed.abundance.txt" 
    params:
        prefix = out_folder + "cupcake/" + data_code + ".cupcake",
        prefix_collapsed = out_folder + "cupcake/" + data_code + ".cupcake.collapsed",
        cupcake_dir = "/sc/arion/projects/ad-omics/data/software/cDNA_Cupcake/build/scripts-3.7/"
    shell:
        #"python {params.cupcake_dir}/collapse_isoforms_by_sam.py --input {input.fasta} "
        "collapse_isoforms_by_sam.py "
        "-b {input.bam} -o {params.prefix} --gen_mol_count --cpus {cupcake_cores};"
        # get abundances
        #"python {params.cupcake_dir}/get_abundance_post_collapse.py {params.prefix}.collapsed {input.cluster_report};"
        # filter 5' truncated transcripts
        #"filter_away_subset.py {params.prefix}.collapsed ; "
        # collapse again to create groups file
        #"collapse_isoforms_by_sam.py --input {output.fasta2} "
                #"-s {input.sam_sorted} --dun-merge-5-shorter -o {params.prefix}.collapsed;"
        # get abundance counts
        #"get_abundance_post_collapse.py {params.prefix}.collapsed.collapsed {input.cluster_report};"
        # copy files to chain directory - omit sample name from file name
        #"cp {output.gff} {output.chain_gff};"
        #"cp {output.group} {output.chain_group};"
        #"cp {output.count} {output.chain_count}"

## DEMULTIPLEX
rule demultiplex_abundances:
    input:
        fasta = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.rep.fa",
        flnc_report = out_folder + "cupcake/" + data_code + ".merged.flnc_report.csv",
        read_stat = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.read_stat.txt",
    output:
        out_folder + "cupcake/" + data_code + ".demux_fl_count.csv"
    params:
        script = "/sc/arion/projects/ad-omics/data/software/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py"
    shell:    
        "python {params.script} --mapped_fafq {input.fasta} --read_stat {input.read_stat} --classify_csv {input.flnc_report} -o {output} "

## FILTER MISSINGNESS
# remove all transcripts with greater than X% missingness
# currently - present in at least 2 samples
# remove monoexonic transcripts to reduce overhead for SQANTI
rule filter_missingness:
    input:
        counts = out_folder + "cupcake/" + data_code + ".demux_fl_count.csv",
        gff = out_folder + "cupcake/" + data_code + ".cupcake.collapsed.gff"
    output:
        counts = out_folder + "cupcake_filtered/all_samples.demux_fl_count_filtered.csv",
        gff = out_folder + "cupcake_filtered/all_samples.cupcake.collapsed.filtered.gtf"
    params:
        script = "scripts/filter_missingness.R",
        input_prefix = out_folder + "cupcake/" + data_code + "",
        output_prefix = out_folder + "cupcake_filtered/all_samples",
        min_samples = 2, # eventually put in config
        min_reads = 2
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} -i {params.input_prefix} -o {params.output_prefix} --min_samples {params.min_samples} --min_reads {params.min_reads} --remove_monoexons"

## SQANTI
# classifies each transcript in the GTF
rule SQANTI_all:
    input:
        gff = out_folder + "cupcake_filtered/all_samples.cupcake.collapsed.filtered.gtf",
        abundance = out_folder + "cupcake_filtered/all_samples.demux_fl_count_filtered.csv"
    output:
        #"test.txt"
        fasta = out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_corrected.fasta",
        gtf = out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_corrected.gtf",
        report = out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt"
    params:
        sample = out_folder + ".cupcake.collapsed.filtered",
        outDir = out_folder + "SQANTI3/",
        nCores = 16,
        nChunks = 12,
        software= "/sc/arion/projects/ad-omics/data/software",
        #junctions = junctionArgs,
        #junctions = "\'" + junctionFolder + "/*SJ.out.tab\'" ,
        gtf = ref_gtf,
        #genome = referenceFa + ".fa",
        intropolis = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified",
        cage = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/hg38.cage_peak_phase1and2combined_coord.bed",
        polya = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/human.polyA.list.txt",
        isoAnnotGFF = "/sc/arion/projects/ad-omics/data/references/hg38_reference/RefSeq/Homo_sapiens_GRCh38_RefSeq_78.gff3"
    shell:
        "conda activate SQANTI3.env; module purge;"
        "export PYTHONPATH=$PYTHONPATH:{params.software}/cDNA_Cupcake/sequence;"
        "export PYTHONPATH=$PYTHONPATH:{params.software}/cDNA_Cupcake/;"
        "python {params.software}/SQANTI3/sqanti3_qc.py -t {params.nCores} --aligner_choice=minimap2"
        " --dir {params.outDir} "
        " --out {params.sample} "
        #" -c {params.junctions} "
        " --cage_peak {params.cage} --polyA_motif_list {params.polya} " 
        #"--skipORF " # ORF finding is slow, can skip if testing
        #"-c {params.intropolis}"
        " --fl_count {input.abundance}"
        " --gtf {input.gff} "
        #" --isoAnnotLite --gff3 {params.isoAnnotGFF}"
        " {params.gtf} {genome} "

rule SQANTI_all_filter:
    input:
        classification = out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt",
        fasta = out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_corrected.fasta",
        #sam = "{sample}/cupcake/{sample}.renamed_corrected.sam",
        gtf = out_folder + "SQANTI3/all_samples.cupcake.collapsed.filtered_corrected.gtf"
        #faa = out_folder + "SQANTI3/all_samples.chained_corrected.faa"
    output:
        out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",
        out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite.gtf"
    params:
        python = "/sc/arion/work/$USER/conda/envs/isoseq-pipeline/bin/python",
        software= "/sc/arion/projects/ad-omics/data/software"
    shell:
        "mkdir -p all_samples/SQANTI3_filtered; "
        "conda activate SQANTI3.env; module purge;"
        "export PYTHONPATH=$PYTHONPATH:{params.software}/cDNA_Cupcake/sequence;"
        "export PYTHONPATH=$PYTHONPATH:{params.software}/cDNA_Cupcake/;"
        "python {params.software}/SQANTI3/sqanti3_RulesFilter.py "
        #" --faa {input.faa} " #--sam {input.sam} "
        " {input.classification} {input.fasta} {input.gtf} "
        " -r 6 -a 0.6 " # intra-priming
        " --filter_mono_exonic " # remove all mono-exons
        "; "   
        " mv all_samples/SQANTI3/*png all_samples/SQANTI3_filtered/;"
        " mv all_samples/SQANTI3/*filtered*lite* all_samples/SQANTI3_filtered/ "


#### IsoAnnot

rule isoaAnotLite:
    input:
        classification = out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",
        gtf = out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite.gtf",
        junctions = out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_junctions.txt",
        isoAnnotGFF = "/sc/arion/projects/ad-omics/data/references/hg38_reference/RefSeq/Homo_sapiens_GRCh38_RefSeq_78.gff3"
    output:
        out_folder + "IsoAnnot/all_samples.filtered_lite_tappAS_annot_from_SQANTI3.gff3"
    params:
        software =  "/sc/arion/projects/ad-omics/data/software",
        prefix = out_folder + "IsoAnnot/all_samples.filtered_lite"
    shell:
        "mkdir -p all_samples/IsoAnnot/;"
        "conda activate SQANTI3.env; module purge;"
        "python {params.software}/SQANTI3/utilities/IsoAnnotLite_SQ3.py {input.gtf} {input.classification} {input.junctions} -gff3 {input.isoAnnotGFF} -o {params.prefix}"
#### MISC

rule collapseAnnotation:
    input: ref_gtf 
    output: ref_gtf + ".genes"
    params: script = "scripts/collapse_annotation.py"
    shell: "/sc/arion/work/humphj04/conda/envs/isoseq-pipeline/bin/python {params.script} {input} {output}"

rule rnaseqc:
    input:
        geneGTF = ref_gtf + ".genes",
        bam =  "{sample}/minimap/{sample}_sorted.bam"
    params:
        out =  "{sample}/qc/"
    output:
         "{sample}/qc/{sample}.metrics.tsv"
    shell:
        "ml rnaseqc;"
        "rnaseqc {input.geneGTF} {input.bam} {params.out} "
        " --sample={wildcards.sample} "
        " --unpaired --coverage --verbose --mapping-quality 0 --base-mismatch=1000 --detection-threshold=1"

# multiqc, version 1.8.dev0 works with rnaseqc outputs
rule multiQC:
    input:
        expand("{sample}/qc/{sample}.metrics.tsv", sample = samples),
        expand("{sample}/qc/{sample}.flagstat.txt", sample = samples),
        expand("{sample}/qc/{sample}.idxstat.txt", sample = samples),
        expand("{sample}/qc/{sample}.hq_fastqc.html", sample = samples)
    output:
         "multiqc/multiqc_report.html"
    shell:
        "export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8;"
        "multiqc -f --outdir multiqc/ ." 

# sort and tabix index final GFF
rule indexGFF:
    input:
        out_folder + "SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite.gtf"
    output:
        gff = out_folder + "SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz",
        index = out_folder + "SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz.tbi"
    params:
        gff3sort = "/sc/arion/projects/ad-omics/data/software/gff3sort/gff3sort.pl"
    shell:
        "{params.gff3sort} {input} | bgzip > {output.gff};"
        "ml tabix;"
        "tabix {output.gff} "



