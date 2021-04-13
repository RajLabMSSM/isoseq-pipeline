
cupcake_path = "/sc/arion/projects/ad-omics/data/software/cDNA_Cupcake"

shell.prefix('export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoseq-pipeline; module purge')
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
rawFolder = config['rawFolder']
dataCode = config['dataCode']

referenceFa = config['referenceFa']
referenceGTF = config['referenceGTF']

primers = "reference/NEB_primers_01_2019.fa"

# short read junctions
junctionFolder = "/sc/arion/projects/als-omics/microglia_isoseq/short_read/junctions"
junctions =  glob.glob(junctionFolder + "/*SJ.out.tab")
junctionArgs = " ".join(["-c " + f for f in junctions])

localrules: create_chain_config

chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "M"]

rule all:
    input:
      dataCode + "/flnc_bam/all_samples.flnc.aligned.bam",
      dataCode + "/cupcake/all_samples.demux_fl_count.csv"
      #dataCode + "/SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt",
      #dataCode + "/SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz", 
      #dataCode + "/SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",
      
#  "test.txt",
         #dataCode + "/SQANTI3/all_samples.cupcake.collapsed_classification.txt",
        #expand(dataCode + "/isoseq3-cluster/all_samples.chr{chr}.hq.fasta", chr = chromosomes + ["U"] ),
        #dataCode + "/TAMA/all_samples_merge.txt",
       #"TAMA_merge/tama_merge_config.txt",
       #expand("{sample}/TAMA/{sample}.bed", sample = samples), 
    #   expand( "{sample}/cupcake/{sample}.cupcake.abundance.txt", sample = samples),
#       expand( "{sample}/SQANTI3/{sample}.{method}_classification.txt", sample = samples, method = ["stringtie","cupcake"]),
        #expand( "{sample}/qc/{sample}.metrics.tsv", sample = dataCode + "")
        #dataCode + "/SQANTI3/all_samples.chained_classification.txt",
    #   expand( "{sample}/stringtie/{sample}.stringtie.gtf", sample = samples)
#        dataCode + "/SQANTI3_filtered/all_samples.chained_classification.filtered_lite_classification.txt",
#        dataCode + "/SQANTI3_filtered/all_samples.chained_classification.filtered.sorted.gff.gz.tbi",
    #   expand(fastqFolder + "{sample}.classification.txt", sample = samples),
         #"multiqc/multiqc_report.html",

# CCS - call circular consensus sequences from the SMRTcell movies
# todo: add chunking for paralellisation
#Input Filter Options:
#  --min-passes   INT    Minimum number of full-length subreads required to generate CCS for a ZMW. [3]
#  --min-snr      FLOAT  Minimum SNR of subreads to use for generating CCS [2.5]

#Draft Filter Options:
#  --min-length   INT    Minimum draft length before polishing. [10]
#  --max-length   INT    Maximum draft length before polishing. [50000]

#Output Filter Options:
#  --min-rq       FLOAT  Minimum predicted accuracy in [0, 1]. [0.99]   
rule CCS:
    input: 
        rawFolder + "{sample}/{sample}.subreadset.xml"
    output: 
        bam =  "{sample}/isoseq3-ccs/{sample}.ccs.bam",
        report =  "{sample}/isoseq3-ccs/{sample}.ccs.report.txt"
    shell:
        "ccs -j 0 {input} {output.bam} --report-file {output.report}"

# primer removal and demultiplexing
rule isoseq_lima:
    input:
         "{sample}/isoseq3-ccs/{sample}.ccs.bam",
    output:
         "{sample}/isoseq3-lima/{sample}.fl.bam"
    shell:
        "lima --isoseq --different --min-passes 1 --split-bam-named --dump-clips --dump-removed -j 0 {input} {primers} {output}"

# trimming of polya tails and removal of concatemers
rule isoseq_refine:
    input:
         "{sample}/isoseq3-lima/{sample}.fl.bam"
    output:
         "{sample}/isoseq3-refine/{sample}.flnc.bam"
    shell:
        "isoseq3 refine --require-polya {input} {primers} {output}"

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
        dataCode + "/flnc_bam/all_samples.merged.flnc_report.csv"
    shell:
        "ml R/3.6.0; Rscript {params.script} {input} {output}" 

# merge flnc.bam files
rule merge_flnc_bams:
    input:
        metaDF["flnc_bam_path"]
    output:
        dataCode + "/flnc_bam/all_samples.flnc.bam"
    params:
        bams = " -in ".join(metaDF["flnc_bam_path"])
    shell:
        "ml bamtools;bamtools merge -in {params.bams} -out {output}"

# align to reference genome using pbmm2
rule align_flnc_bam:
    input:
        bam = dataCode + "/flnc_bam/all_samples.flnc.bam",
        mmi = referenceFa + ".mmi"
    output:
        bam = dataCode + "/flnc_bam/all_samples.flnc.aligned.bam"
    shell:
        "pbmm2 align --sort -j 32 --sort-threads 4 -m 3G --preset=ISOSEQ --log-level INFO --unmapped {input.mmi} {input.bam} {output.bam}"

# get bam statistics
rule samtools:
    input:  dataCode + "/flnc_bam/all_samples.flnc.aligned.bam"
    output:
        bam =  dataCode + "/minimap/all_samples.sorted.bam",
        bai =  dataCode + "/minimap/all_samples.sorted.bam.bai",
        flagstat =  dataCode + "/qc/all_samples.flagstat.txt",
        idxstat =  dataCode + "/qc/all_samples.idxstat.txt"
    shell:
        "samtools view -bh {input} | samtools sort > {output.bam}; "
        "samtools index {output.bam};"
        "samtools flagstat {output.bam} > {output.flagstat};"
        "samtools idxstats {output.bam} > {output.idxstat} "

# make FASTA version of BAM for Cupcake
rule make_fasta:
    input: 
        dataCode + "/flnc_bam/all_samples.flnc.aligned.bam"
    output:
        dataCode + "/flnc_bam/all_samples.flnc.aligned.fasta"
    shell:
        "ml samtools;"
        "samtools fasta {input} > {output}"

# make SAM version of BAM for Cupcake
rule make_sam:
    input:
        dataCode + "/flnc_bam/all_samples.flnc.aligned.bam"
    output:
        dataCode + "/flnc_bam/all_samples.flnc.aligned.sam"
    shell:
        "ml samtools;"
        "samtools view -h {input} > {output}"
## Cupcake Tools

# collapse redundant reads
# deal with 5' truncated reads - don't use dun-merge-5-shorter or filter_away
# filter away  
# CHECK - does cupcake mind if FASTA and cluster_report are gzipped?
rule cupcake_collapse:
    input:
        fasta = dataCode + "/flnc_bam/all_samples.flnc.aligned.fasta",
        #fasta =   dataCode + "/prepend/all_samples.merged.hq.prepend.fasta",
        sam_sorted = dataCode + "/flnc_bam/all_samples.flnc.aligned.sam"
    output:
         gff = dataCode + "/cupcake/all_samples.cupcake.collapsed.gff",
         fasta = dataCode + "/cupcake/all_samples.cupcake.collapsed.rep.fa",
         #group = dataCode + "/cupcake/all_samples.cupcake.collapsed.group.txt",
         stat = dataCode + "/cupcake/all_samples.cupcake.collapsed.read_stat.txt",
         #count = dataCode + "/cupcake/all_samples.cupcake.collapsed.abundance.txt"
         #fasta2 = dataCode + "/cupcake/all_samples.cupcake.collapsed.rep.fa",
         #chain_gff = dataCode + "/cupcake/chain/cupcake.collapsed.gff",
         #chain_group = dataCode + "/cupcake/chain/cupcake.collapsed.group.txt",
         #chain_count = dataCode + "/cupcake/chain/cupcake.collapsed.abundance.txt" 

    params:
        prefix = dataCode + "/cupcake/all_samples.cupcake",
        prefix_collapsed = dataCode + "/cupcake/all_samples.cupcake.collapsed",
        cupcake_dir = "/sc/arion/projects/ad-omics/data/software/cDNA_Cupcake/build/scripts-3.7/"
    shell:
        #"python {params.cupcake_dir}/collapse_isoforms_by_sam.py --input {input.fasta} "
        "collapse_isoforms_by_sam.py --input {input.fasta} "
        "-s {input.sam_sorted} -o {params.prefix};"
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
        fasta = dataCode + "/cupcake/all_samples.cupcake.collapsed.rep.fa",
        flnc_report = dataCode + "/flnc_bam/all_samples.merged.flnc_report.csv",
        read_stat = dataCode + "/cupcake/all_samples.cupcake.collapsed.read_stat.txt",
    output:
        dataCode + "/cupcake/all_samples.demux_fl_count.csv"
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
        counts = dataCode + "/cupcake/all_samples.demux_fl_count.csv",
        gff = dataCode + "/cupcake/all_samples.cupcake.collapsed.gff"
    output:
        counts = dataCode + "/cupcake_filtered/all_samples.demux_fl_count_filtered.csv",
        gff = dataCode + "/cupcake_filtered/all_samples.cupcake.collapsed.filtered.gtf"
    params:
        script = "scripts/filter_missingness.R",
        input_prefix = dataCode + "/cupcake/all_samples",
        output_prefix = dataCode + "/cupcake_filtered/all_samples",
        min_samples = 2, # eventually put in config
        min_reads = 2
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} -i {params.input_prefix} -o {params.output_prefix} --min_samples {params.min_samples} --min_reads {params.min_reads} --remove_monoexons"

## SQANTI
# classifies each transcript in the GTF
rule SQANTI_all:
    input:
        gff = dataCode + "/cupcake_filtered/all_samples.cupcake.collapsed.filtered.gtf",
        abundance = dataCode + "/cupcake_filtered/all_samples.demux_fl_count_filtered.csv"
    output:
        #"test.txt"
        fasta = dataCode + "/SQANTI3/all_samples.cupcake.collapsed.filtered_corrected.fasta",
        gtf = dataCode + "/SQANTI3/all_samples.cupcake.collapsed.filtered_corrected.gtf",
        report = dataCode + "/SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt"
    params:
        sample = dataCode + ".cupcake.collapsed.filtered",
        outDir = dataCode + "/SQANTI3/",
        nCores = 16,
        nChunks = 12,
        software= "/sc/arion/projects/ad-omics/data/software",
        #junctions = junctionArgs,
        junctions = "\'" + junctionFolder + "/*SJ.out.tab\'" ,
        gtf = referenceGTF,
        genome = referenceFa + ".fa",
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
        " -c {params.junctions} "
        " --cage_peak {params.cage} --polyA_motif_list {params.polya} " 
        #"--skipORF " # ORF finding is slow, can skip if testing
        #"-c {params.intropolis}"
        " --fl_count {input.abundance}"
        " --gtf {input.gff} "
        #" --isoAnnotLite --gff3 {params.isoAnnotGFF}"
        " {params.gtf} {params.genome} "

rule SQANTI_all_filter:
    input:
        classification = dataCode + "/SQANTI3/all_samples.cupcake.collapsed.filtered_classification.txt",
        fasta = dataCode + "/SQANTI3/all_samples.cupcake.collapsed.filtered_corrected.fasta",
        #sam = "{sample}/cupcake/{sample}.renamed_corrected.sam",
        gtf = dataCode + "/SQANTI3/all_samples.cupcake.collapsed.filtered_corrected.gtf"
        #faa = dataCode + "/SQANTI3/all_samples.chained_corrected.faa"
    output:
        dataCode + "/SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",
        dataCode + "/SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite.gtf"
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
        classification = dataCode + "/SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_classification.txt",
        gtf = dataCode + "/SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite.gtf",
        junctions = dataCode + "/SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite_junctions.txt",
        isoAnnotGFF = "/sc/arion/projects/ad-omics/data/references/hg38_reference/RefSeq/Homo_sapiens_GRCh38_RefSeq_78.gff3"
    output:
        dataCode + "/IsoAnnot/all_samples.filtered_lite_tappAS_annot_from_SQANTI3.gff3"
    params:
        software =  "/sc/arion/projects/ad-omics/data/software",
        prefix = dataCode + "/IsoAnnot/all_samples.filtered_lite"
    shell:
        "mkdir -p all_samples/IsoAnnot/;"
        "conda activate SQANTI3.env; module purge;"
        "python {params.software}/SQANTI3/utilities/IsoAnnotLite_SQ3.py {input.gtf} {input.classification} {input.junctions} -gff3 {input.isoAnnotGFF} -o {params.prefix}"
#### MISC

rule collapseAnnotation:
    input: referenceGTF 
    output: referenceGTF + ".genes"
    params: script = "scripts/collapse_annotation.py"
    shell: "/sc/arion/work/humphj04/conda/envs/isoseq-pipeline/bin/python {params.script} {input} {output}"

rule rnaseqc:
    input:
        geneGTF = referenceGTF + ".genes",
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
        dataCode + "/SQANTI3_filtered/all_samples.cupcake.collapsed.filtered_classification.filtered_lite.gtf"
    output:
        gff = dataCode + "/SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz",
        index = dataCode + "/SQANTI3_filtered/all_samples.filtered.sorted.gtf.gz.tbi"
    params:
        gff3sort = "/sc/arion/projects/ad-omics/data/software/gff3sort/gff3sort.pl"
    shell:
        "{params.gff3sort} {input} | bgzip > {output.gff};"
        "ml tabix;"
        "tabix {output.gff} "

