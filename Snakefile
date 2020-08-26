
cupcake_path = "/sc/hydra/projects/ad-omics/data/software/cDNA_Cupcake"

shell.prefix('export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; ml purge; conda activate isoseq-pipeline;')
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
junctionFolder = "/sc/hydra/projects/ad-omics/microglia_isoseq/short_read_junctions"
junctions =  glob.glob(junctionFolder + "/*SJ.out.tab")
junctionArgs = " ".join(["-c " + f for f in junctions])

localrules: create_chain_config

chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "M"]

rule all:
    input:
        expand("all_samples/isoseq3-cluster/all_samples.chr{chr}.hq.fasta", chr = chromosomes + ["U"] ),
        #"all_samples/TAMA/all_samples_merge.txt",
       #"TAMA_merge/tama_merge_config.txt",
       #expand("{sample}/TAMA/{sample}.bed", sample = samples), 
    #   expand( "{sample}/cupcake/{sample}.cupcake.abundance.txt", sample = samples),
#       expand( "{sample}/SQANTI2/{sample}.{method}_classification.txt", sample = samples, method = ["stringtie","cupcake"]),
        expand( "{sample}/qc/{sample}.metrics.tsv", sample = "all_samples")
        #"all_samples/SQANTI2/all_samples.chained_classification.txt",
    #   expand( "{sample}/stringtie/{sample}.stringtie.gtf", sample = samples)
#        "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered_lite_classification.txt",
#        "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered.sorted.gff.gz.tbi",
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

# merge flnc.bam files
rule merge_flnc_bams:
    input:
        metaDF["flnc_bam_path"]
    output:
        "all_samples/flnc_bam/all_samples.flnc.bam"
    params:
        bams = " -in ".join(metaDF["flnc_bam_path"])
    shell:
        "ml bamtools;bamtools merge -in {params.bams} -out {output}"

# align to reference genome using pbmm2
rule align_flnc_bam:
    input:
        bam = "all_samples/flnc_bam/all_samples.flnc.bam",
        mmi = referenceFa + ".mmi"
    output:
        bam = "all_samples/flnc_bam/all_samples.flnc.aligned.bam"
    shell:
        "pbmm2 align --sort -j 32 --sort-threads 4 -m 3G --preset=ISOSEQ --log-level INFO --unmapped {input.mmi} {input.bam} {output.bam}"

# split bams by chromosome
rule split_flnc_bam:
    input:
        bam = "all_samples/flnc_bam/all_samples.flnc.aligned.bam"
    output:
        expand("all_samples/flnc_bam/all_samples.flnc.aligned.chr{chr}.bam", chr = chromosomes )
    shell:
        "ml bamtools;"
        "bamtools split -in {input.bam} -reference  -refPrefix \"\" "

# get out unmapped reads for clustering
rule split_unmapped_flnc_bam:
    input:
        bam = "all_samples/flnc_bam/all_samples.flnc.aligned.bam"
    output:
        bam = "all_samples/flnc_bam/all_samples.flnc.aligned.chrU.bam"
    shell:
        "ml samtools;"
        "samtools view -bh -f 4 {input.bam} > {output.bam}" 


# cluster and polish together
rule isoseq3_cluster:
    input:
        "all_samples/flnc_bam/all_samples.flnc.aligned.chr{chr}.bam"
        #"all_samples/flnc_bam/all_samples.flnc.bam"
        #"{sample}/isoseq3-refine/{sample}.flnc.bam"
    params:
        fasta_gz = "all_samples/isoseq3-cluster/all_samples.chr{chr}.hq.fasta.gz",
        bam = "all_samples/isoseq3-cluster/all_samples.chr{chr}.bam"
        #fasta_gz =  "{sample}/isoseq3-cluster/{sample}.hq.fasta.gz"
        
    output:
        fasta = "all_samples/isoseq3-cluster/all_samples.chr{chr}.hq.fasta",
        report = "all_samples/isoseq3-cluster/all_samples.chr{chr}.cluster_report.csv"
        #fasta =  "{sample}/isoseq3-cluster/{sample}.hq.fasta",
        #report =  "{sample}/isoseq3-cluster/{sample}.polished.cluster_report.csv"
    shell:
        "isoseq3 cluster --verbose --use-qvs -j 0 {input} {params.bam};"
        "gunzip {params.fasta_gz}"

# prepend chrom name to clustered FASTA lines before merging
rule prepend_chr_name:
    input:
         fasta = "all_samples/isoseq3-cluster/all_samples.chr{chr}.hq.fasta"
    output:
        fasta =  "all_samples/isoseq3-cluster/all_samples.chr{chr}.hq.prepend.fasta"
    shell:
        "ml bioawk;"
        "bioawk -c fastx '{{ print $name\"/chr{wildcards.chr}\"; print $seq}}' {input.fasta} > {output.fasta}"

# merge FASTA together for full alignment
rule merge_fasta:
    input:
        expand("all_samples/isoseq3-cluster/all_samples.chr{chr}.hq.prepend.fasta", chr = chromosomes + ["U"] )
    output:
        "all_samples/isoseq3-cluster/all_samples.merged.hq.prepend.fasta"    
    shell:
        "cat {input} > {output}"

# currently isoseq3 and clustering is outsourced to Nancy
rule symlinkFiles:
    input: 
        fasta = sampleFA,
        report = sampleREP
    output:
        expand("{sample}/isoseq3-cluster/{sample}.hq.fasta", sample = samples),
        expand("{sample}/isoseq3-cluster/{sample}.cluster_report.csv", sample = samples)
    run:
        for s in samples:
            fa_in = metaDF.loc[s]['fasta_path']
            rep_in = metaDF.loc[s]['cluster_report_path']
            
            fa_out = s + "/isoseq3-cluster/" + s + ".hq.fasta"
            rep_out = s + "/isoseq3-cluster/" + s + ".cluster_report.csv"

            os.makedirs(s + "/isoseq3-cluster/", exist_ok = True)

            if not os.path.exists(fa_out):
                os.symlink(fa_in, fa_out)   
            if not os.path.exists(rep_out):
                os.symlink(rep_in, rep_out)

#### MINIMAP 

rule minimapIndex:
    input: referenceFa + ".fa"
    output: referenceFa + ".mmi"
    shell:
        "minimap2 -d {output} {input}" 

rule minimap:
    input: 
        fastq = "all_samples/isoseq3-cluster/all_samples.merged.hq.prepend.fasta",
        #fastq =  "{sample}/isoseq3-cluster/{sample}.hq.fasta",
        ref = referenceFa + ".fa",
        index = referenceFa + ".mmi"
    params: 
        "-ax splice -t 4 -uf --secondary=no -C5"
        #config['minimapParams']
    output: 
        #sam = "all_samples/minimap/all_samples.hq.sam",
        #sam_sorted = "all_samples/minimap/all_samples.hq.sorted.sam"
        sam =  "{sample}/minimap/{sample}.hq.sam",
        sam_sorted =  "{sample}/minimap/{sample}.hq.sorted.sam"
    shell:
        "minimap2 {params} {input.index} {input.fastq} > {output.sam};"
        "sort -k 3,3 -k 4,4n {output.sam} > {output.sam_sorted}"

rule samtools:
    input:  "{sample}/minimap/{sample}.hq.sam"
    output:
        bam =  "{sample}/minimap/{sample}_sorted.bam",
        bai =  "{sample}/minimap/{sample}_sorted.bam.bai",
        flagstat =  "{sample}/qc/{sample}.flagstat.txt",
        idxstat =  "{sample}/qc/{sample}.idxstat.txt"
    shell:
        "samtools view -bh {input} | samtools sort > {output.bam}; "
        "samtools index {output.bam};"
        "samtools flagstat {output.bam} > {output.flagstat};"
        "samtools idxstats {output.bam} > {output.idxstat} "

rule fastqc:
    input: "{sample}/minimap/{sample}.hq.sam"
    output: "{sample}/qc/{sample}.hq_fastqc.html"
    shell:
        "ml fastqc;"
        "fastqc --outdir={wildcards.sample}/qc/ --format sam {input}"

## TAMA tools

rule TAMA_collapse:
    input:
        sam_sorted = "{sample}/minimap/{sample}.hq.sorted.sam",
        genome = referenceFa + ".fa"
    output:
        bed = "{sample}/TAMA/{sample}.bed",
        txt = "{sample}/TAMA/{sample}_read.txt"
    params:
        script = "/sc/hydra/projects/ad-omics/data/software/tama/tama_collapse.py"
    shell:
        'conda activate py2bio;'
        'python {params.script} -s {input.sam_sorted} -f {input.genome} -p {wildcards.sample}/TAMA/{wildcards.sample} -x no_cap -rm low_mem'


rule create_TAMA_merge_config:
    input:
        expand( "{sample}/TAMA/{sample}.bed", sample = samples)
    output:
        config = "TAMA_merge/tama_merge_config.txt"
    run:
        tamaMergeRows = []
        for s in samples:
            entry = s + "/TAMA/" + s + ".bed\tno_cap\t1,1,1\t" + s
            tamaMergeRows.append(entry)

        with open(output.config, 'w') as filehandle:
                for listitem in tamaMergeRows:
                    filehandle.write('%s\n' % listitem)

rule TAMA_merge:
    input:
        config = "TAMA_merge/tama_merge_config.txt"
    output:
        "all_samples/TAMA/all_samples_merge.txt"
    params:
        script = "/sc/hydra/projects/ad-omics/data/software/tama/tama_merge.py"
    shell:
        'conda activate py2bio;'
        'python {params.script} -f {input.config} -p all_samples'



## Cupcake Tools

# collapse redundant reads
# deal with 5' truncated reads - don't use dun-merge-5-shorter or filter_away
# filter away  
# CHECK - does cupcake mind if FASTA and cluster_report are gzipped?
rule cupcake_collapse:
    input:
        fasta =   "{sample}/isoseq3-cluster/{sample}.hq.fasta",
        sam_sorted =  "{sample}/minimap/{sample}.hq.sorted.sam",
        cluster_report =  "{sample}/isoseq3-cluster/{sample}.cluster_report.csv"
    output:
         gff = "{sample}/cupcake/{sample}.cupcake.collapsed.gff",
         fasta = "{sample}/cupcake/{sample}.cupcake.collapsed.rep.fa",
         group = "{sample}/cupcake/{sample}.cupcake.collapsed.group.txt",
         stat = "{sample}/cupcake/{sample}.cupcake.collapsed.read_stat.txt",
         count = "{sample}/cupcake/{sample}.cupcake.collapsed.abundance.txt"
         #fasta2 = "{sample}/cupcake/{sample}.cupcake.collapsed.rep.fa",
         #chain_gff = "{sample}/cupcake/chain/cupcake.collapsed.gff",
         #chain_group = "{sample}/cupcake/chain/cupcake.collapsed.group.txt",
         #chain_count = "{sample}/cupcake/chain/cupcake.collapsed.abundance.txt" 

    params:
        prefix = "{sample}/cupcake/{sample}.cupcake",
        prefix_collapsed = "{sample}/cupcake/{sample}.cupcake.collapsed"
    shell:
        "collapse_isoforms_by_sam.py --input {input.fasta} "
        "-s {input.sam_sorted} -o {params.prefix};"
        # get abundances
        "get_abundance_post_collapse.py {params.prefix}.collapsed {input.cluster_report};"
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


rule create_chain_config:
    output:
        config = "all_samples/chain.config.txt"
    run:
        chainSampleRows = []
        for i in samples:
            l = "".join(["SAMPLE=", i,";", i, "/cupcake/chain/"])
            chainSampleRows.append(l)
        chainFileRows = ["GROUP_FILENAME=cupcake.collapsed.group.txt", \
                         "COUNT_FILENAME=cupcake.collapsed.abundance.txt", \
                         "GFF_FILENAME=cupcake.collapsed.gff" ]
        allRows = chainSampleRows + [""] + chainFileRows    
        
        with open(output.config, 'w') as filehandle:
                for listitem in allRows:
                    filehandle.write('%s\n' % listitem)


#### SQANTI - per-sample if chaining won't work

rule SQANTI:
    input:
        #fasta = "{sample}/isoseq3-cluster/{sample}.hq.fasta"
        gff = "{sample}/{method}/{sample}.{method}.collapsed.filtered.gff"
    output:
        report = "{sample}/SQANTI2/{sample}.{method}_classification.txt"
    params:
        sample = "{sample}.{method}",
        outDir = "{sample}/SQANTI2/",
        python = "/sc/arion/work/$USER/conda/envs/isoseq-pipeline/bin/python",
        sqantiPath= "/sc/arion/projects/ad-omics/data/software/SQANTI2",
        nCores = 12,
        #abundance = "{sample}/cupcake/{sample}.{method}.abundance.txt",
        gtf = referenceGTF,
        genome = referenceFa + ".fa",
        junctions = "\'" + junctionFolder + "/*SJ.out.tab\'" ,
        intropolis = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/SQANTI2/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified",
        cage = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/SQANTI2/hg38.cage_peak_phase1and2combined_coord.bed",
        polya = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/SQANTI2/human.polyA.list.txt"
    shell:
        #"export PATH=/sc/hydra/projects/ad-omics/data/software/UCSC/:$PATH;"
        #"module unload gcc;ml R/3.6.0; "
        "export PYTHONPATH=$PYTHONPATH:/hpc/users/humphj04/pipelines/cDNA_Cupcake/sequence/;"
        "{params.python} {params.sqantiPath}/sqanti_qc2.py -t {params.nCores} --aligner_choice=minimap2"
        " --dir {params.outDir} "
        " --out {params.sample} "
        " -c {params.junctions} "
        " --cage_peak {params.cage} --polyA_motif_list {params.polya} -c {params.intropolis}"
        #" --fl_count {params.abundance}"
        " --gtf {input.gff} " 
        " {params.gtf} {params.genome} "

# sqanti filter - only works with Cupcake output for now
rule SQUANTI_filter:
    input:
        classification = "{sample}/SQANTI2/{sample}.{method}_classification.txt",
        fasta = "{sample}/SQANTI2/{sample}.cupcake.collapsed_corrected.fasta",
        #sam = "{sample}/cupcake/{sample}.renamed_corrected.sam",
        gtf = "{sample}/SQANTI2/{sample}.cupcake.collapsed_corrected.gtf",
        faa = "{sample}/SQANTI2/{sample}.cupcake.collapsed_corrected.faa"
    output:
        "{sample}/SQANTI2/{sample}_classification.filtered_lite_classification.txt"
    params:
        python = "/sc/hydra/work/$USER/conda/envs/isoseq-pipeline/bin/python",
        sqantiPath= "/sc/hydra/projects/ad-omics/data/software/SQANTI2"
    shell:
        "{params.python} {params.sqantiPath}/sqanti_filter2.py "
        " --faa {input.faa} " #--sam {input.sam} "
        " {input.classification} {input.fasta} {input.gtf} "    

rule chain_samples:
    input:
        chain_gff = expand("{sample}/cupcake/chain/cupcake.collapsed.gff", sample = samples),
        chain_config = "all_samples/chain.config.txt"
    params:
        n_cores = 4
    output:
        gff = "all_samples/all_samples.chained.gff",
    shell:
        "chain_samples.py {input.chain_config} count_fl --cpus {params.n_cores}; "
        "if [ ! -d all_samples/ ]; then mkdir all_samples; fi ;"
        "mv all_samples.chained* all_samples/;"
        "rm tmp* "

# collapse again after chaining? experimental
rule collapse_post_chain:
    input:
        gff = "all_samples/all_samples.chained.gff"
    output:
        "collapse_isoforms_by_sam.py --input {input.fasta} "
        "-s {input.sam_sorted} -o {params.prefix};"


## SQANTI

rule SQANTI_all:
    input:
        gff = "all_samples/all_samples.chained.gff"
    output:
        fasta = "all_samples/SQANTI2/all_samples.chained_corrected.fasta",
        gtf = "all_samples/SQANTI2/all_samples.chained_corrected.gtf",
        report = "all_samples/SQANTI2/all_samples.chained_classification.txt"
    params:
        sample = "all_samples.chained",
        outDir = "all_samples/SQANTI2/",
        nCores = 32,
        nChunks = 16,
        python = "/sc/hydra/work/$USER/conda/envs/isoseq-pipeline/bin/python",
        sqantiPath= "/sc/hydra/projects/ad-omics/data/software/SQANTI2",
        junctions = "\'" + junctionFolder + "/*SJ.out.tab\'" ,
        abundance = "all_samples/all_samples.chained_count.txt",
        gtf = referenceGTF,
        genome = referenceFa + ".fa",
        intropolis = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/SQANTI2/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified",
        cage = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/SQANTI2/hg38.cage_peak_phase1and2combined_coord.bed",
        polya = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/SQANTI2/human.polyA.list.txt"
    shell:
        "export PYTHONPATH=$PYTHONPATH:/hpc/users/humphj04/pipelines/cDNA_Cupcake/sequence/;"
        "{params.python} {params.sqantiPath}/sqanti_qc2.py -t {params.nCores} --aligner_choice=minimap2"
        " --dir {params.outDir} "
        " --out {params.sample} "
        " -c {params.junctions} "
        " --cage_peak {params.cage} --polyA_motif_list {params.polya} " 
        "--skipORF " # skipping ORF finding for now as it's very slow
        #"-c {params.intropolis}"
        " --fl_count {params.abundance}"
        " --gtf {input.gff} "
        " {params.gtf} {params.genome} "

rule SQUANTI_all_filter:
    input:
        classification = "all_samples/SQANTI2/all_samples.chained_classification.txt",
        fasta = "all_samples/SQANTI2/all_samples.chained_corrected.fasta",
        #sam = "{sample}/cupcake/{sample}.renamed_corrected.sam",
        gtf = "all_samples/SQANTI2/all_samples.chained_corrected.gtf"
        #faa = "all_samples/SQANTI2/all_samples.chained_corrected.faa"
    output:
        "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered_lite_classification.txt",
        "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered_lite.gtf"
    params:
        python = "/sc/hydra/work/$USER/conda/envs/isoseq-pipeline/bin/python",
        sqantiPath= "/sc/hydra/projects/ad-omics/data/software/SQANTI2"
    shell:
        " mkdir -p all_samples/SQANTI2_filtered; "
        "{params.python} {params.sqantiPath}/sqanti_filter2.py "
        #" --faa {input.faa} " #--sam {input.sam} "
        " {input.classification} {input.fasta} {input.gtf} "
        " -r 6 -a 0.6 " # intra-priming
        " --filter_mono_exonic " # remove all mono-exons
        #" --skipGTF --skipFaFq " # don't write out new FA and GTF - very slow currently
        "; "   
        " mv all_samples/SQANTI2/*png all_samples/SQANTI2_filtered/;"
        " mv all_samples/SQANTI2/*filtered*lite* all_samples/SQANTI2_filtered/ "

#### STRINGTIE2

# assemble minimap-aligned reads into transcripts - alternative to cupcake_collapse
# some monoexon transcripts are given strand of "." - remove them in bioawk
rule stringtie:
    input:
        bam =  "{sample}/minimap/{sample}_sorted.bam"
    params:
        gtf = referenceGTF,
        stringtiePath = "/sc/hydra/projects/ad-omics/data/software/stringtie"
    output:
        gtf = "{sample}/stringtie/{sample}.stringtie.collapsed.gtf",
        gff = "{sample}/stringtie/{sample}.stringtie.collapsed.gff"
    shell:
        "ml bioawk; "
        "{params.stringtiePath}/stringtie --rf -G {params.gtf} -L -o {output.gtf} {input.bam};"
        "gffread -E {output.gtf} -o- | awk \'$7 != \".\"\' > {output.gff}"
 
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
        "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered_lite.gtf"
    output:
        gff = "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered.sorted.gff.gz",
        index = "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered.sorted.gff.gz.tbi"
    params:
        gff3sort = "/sc/hydra/projects/ad-omics/data/software/gff3sort/gff3sort.pl"
    shell:
        "{params.gff3sort} {input} | bgzip > {output.gff};"
        "ml tabix;"
        "tabix {output.gff} "


