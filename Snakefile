
shell.prefix('export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; ml purge;')
#shell.prefix('export PS1="";source activate isoseq-pipeline;ml R/3.6.0;')
import pandas as pd
import glob
metadata = config['metadata']

metaDF = pd.read_csv(metadata, sep = '\t')
samples = metaDF['sample']

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

rule all:
    input:
        "TAMA_merge/tama_merge_config.txt",
        expand("{sample}/TAMA/{sample}.bed", sample = samples), 
    #   expand( "{sample}/cupcake/{sample}.cupcake.abundance.txt", sample = samples),
#       expand( "{sample}/SQANTI2/{sample}.{method}_classification.txt", sample = samples, method = ["stringtie","cupcake"]),
        expand( "{sample}/qc/{sample}.metrics.tsv", sample = samples),
        "all_samples/SQANTI2/all_samples.chained_classification.txt",
#   expand( "{sample}/stringtie/{sample}.stringtie.gtf", sample = samples)
        "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered_lite_classification.txt",
        "all_samples/SQANTI2_filtered/all_samples.chained_classification.filtered.sorted.gff.gz.tbi"
    #   expand(fastqFolder + "{sample}.classification.txt", sample = samples),
        # "multiqc/multiqc_report.html",
        #expand( "rnaseqc/{samp}.metrics.tsv", samp = samples),
        #outFolder + "multiqc/multiqc_report.html",
        #expand(outFolder + "rnaseqc/{samp}.metrics.tsv", samp = samples),
        #reference + ".mmi",
        #expand( "sorted/{samp}_sorted.bam", samp = samples)
#       #expand("{outFolder}{samples}_sorted.bam", samples = samples, outFolder = outFolder)


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
# polish
rule isoseq_cluster:
    input:
        "{sample}/isoseq3-refine/{sample}.flnc.bam"
    params:
        fasta_gz =  "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta.gz"
    output:
        fasta =  "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta",
        report =  "{sample}/isoseq3-cluster/{sample}.polished.cluster_report.csv"
    shell:
        "isoseq3 cluster --verbose --use-qvs -j 0 {input} {output.fasta};"
        "gunzip {params.fasta_gz}"

#### MINIMAP 

rule minimapIndex:
    input: referenceFa + ".fa"
    output: referenceFa + ".mmi"
    shell:
        "minimap2 -d {output} {input}" 
rule minimap:
    input: 
        fastq =  "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta",
        ref = referenceFa + ".fa",
        index = referenceFa + ".mmi"
    params: 
        "-ax splice -t 4 -uf --secondary=no -C5"
        #config['minimapParams']
    output: 
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
        'python {params.script} -s {input.sam_sorted} -f {input.genome} -p {wildcards.sample}/TAMA/{wildcards.sample} -x no_cap'

#annotation_capped.bed        capped  1,1,1   caplib
#annotation_nocap.bed        no_cap  2,1,1   nocaplib

rule create_TAMA_merge_config:
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


## Cupcake Tools

# collapse redundant reads
# deal with 5' truncated reads
# filter away   
rule cupcake_collapse:
    input:
        fasta =   "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta",
        sam_sorted =  "{sample}/minimap/{sample}.hq.sorted.sam",
        cluster_report =  "{sample}/isoseq3-cluster/{sample}.polished.cluster_report.csv"
    output:
         gff = "{sample}/cupcake/{sample}.cupcake.collapsed.filtered.gff",
         fasta = "{sample}/cupcake/{sample}.cupcake.collapsed.rep.fa",
         group = "{sample}/cupcake/{sample}.cupcake.collapsed.group.txt",
         stat = "{sample}/cupcake/{sample}.cupcake.collapsed.read_stat.txt",
         count = "{sample}/cupcake/{sample}.cupcake.collapsed.abundance.txt",
         fasta2 = "{sample}/cupcake/{sample}.cupcake.collapsed.filtered.rep.fa",
         chain_gff = "{sample}/cupcake/chain/cupcake.collapsed.filtered.gff",
         chain_group = "{sample}/cupcake/chain/cupcake.collapsed.group.txt",
         chain_count = "{sample}/cupcake/chain/cupcake.collapsed.abundance.txt" 

    params:
        prefix = "{sample}/cupcake/{sample}.cupcake",
        prefix_collapsed = "{sample}/cupcake/{sample}.cupcake.collapsed"
    shell:
        "collapse_isoforms_by_sam.py --input {input.fasta} "
        "-s {input.sam_sorted} --dun-merge-5-shorter -o {params.prefix};"
        # get abundances
        "get_abundance_post_collapse.py {params.prefix}.collapsed {input.cluster_report};"
        # filter 5' truncated transcripts
        "filter_away_subset.py {params.prefix}.collapsed ; "
        # collapse again to create groups file
        #"collapse_isoforms_by_sam.py --input {output.fasta2} "
                #"-s {input.sam_sorted} --dun-merge-5-shorter -o {params.prefix}.collapsed.filtered;"
        # get abundance counts
        #"get_abundance_post_collapse.py {params.prefix}.collapsed.filtered.collapsed {input.cluster_report};"
        # copy files to chain directory - omit sample name from file name
        "cp {output.gff} {output.chain_gff};"
        "cp {output.group} {output.chain_group};"
        "cp {output.count} {output.chain_count}"

chainFileRows = ["GROUP_FILENAME=cupcake.collapsed.group.txt", "COUNT_FILENAME=cupcake.collapsed.abundance.txt", "GFF_FILENAME=cupcake.collapsed.filtered.gff"]

rule create_chain_config:
    output:
        config = "chain.config.txt"
    run:
        chainSampleRows = []
        for i in samples:
            l = "SAMPLE=" + i + ";" + i + "/cupcake/chain/"
            chainSampleRows.append(l)
        allRows = chainSampleRows + [""] + chainFileRows    
        
        with open(output.config, 'w') as filehandle:
                for listitem in allRows:
                    filehandle.write('%s\n' % listitem)

rule chain_samples:
    input:
        chain_gff = expand("{sample}/cupcake/chain/cupcake.collapsed.filtered.gff", sample = samples),
        chain_config = "chain.config.txt"
    output:
        gff = "all_samples/all_samples.chained.gff"
    shell:
        "chain_samples.py {input.chain_config} count_fl; "
        "if [ ! -d all_samples/ ]; then mkdir all_samples; fi ;"
        "mv all_samples.chained* all_samples/;"
        "rm tmp* "


#### SQANTI

rule SQANTI:
    input:
        #fasta = "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta"
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



rule SQANTI_all:
    input:
        gff = "all_samples/all_samples.chained.gff"
    output:
        report = "all_samples/SQANTI2/all_samples.chained_classification.txt"
    params:
        sample = "all_samples.chained",
        outDir = "all_samples/SQANTI2/",
        python = "/sc/hydra/work/$USER/conda/envs/isoseq-pipeline/bin/python",
        sqantiPath= "/sc/hydra/projects/ad-omics/data/software/SQANTI2",
        nCores = 12,
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
        python = "/sc/arion/work/$USER/conda/envs/isoseq-pipeline/bin/python",
        sqantiPath= "/sc/arion/projects/ad-omics/data/software/SQANTI2"
    shell:
        "{params.python} {params.sqantiPath}/sqanti_filter2.py "
        #" --faa {input.faa} " #--sam {input.sam} "
        " {input.classification} {input.fasta} {input.gtf} ;"   
        " mv *png all_samples/SQANTI_filtered/"


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
rule multiQC:
    input:
        "{sample}/qc/{sample}.metrics.tsv",
        "{sample}/qc/{sample}.flagstat.txt",
        "{sample}/qc/{sample}.idxstat.txt"
    output:
         "multiqc/multiqc_report.html"
    shell:
        "export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8;"
        "multiqc -f --outdir {outFolder}multiqc/ {outFolder}" 

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

