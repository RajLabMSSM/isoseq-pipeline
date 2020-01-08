
import pandas as pd

metadata = config['metadata']

metaDF = pd.read_csv(metadata, sep = '\t')
samples = metaDF['sample']

#print(samples)
rawFolder = config['rawFolder']
dataCode = config['dataCode']

referenceFa = config['referenceFa']
referenceGTF = config['referenceGTF']

primers = "reference/NEB_primers_01_2019.fa"

rule all:
	input:
	#	expand( "{sample}/cupcake/{sample}.cupcake.abundance.txt", sample = samples),
		expand( "{sample}/SQANTI2/{sample}.{method}_classification.txt", sample = samples, method = ["stringtie","cupcake"]),
		expand( "{sample}/qc/{sample}.metrics.tsv", sample = samples),
		expand( "{sample}/stringtie/{sample}.stringtie.gtf", sample = samples)
	#	expand(fastqFolder + "{sample}.classification.txt", sample = samples),
		# "multiqc/multiqc_report.html",
		#expand( "rnaseqc/{samp}.metrics.tsv", samp = samples),
		#outFolder + "multiqc/multiqc_report.html",
		#expand(outFolder + "rnaseqc/{samp}.metrics.tsv", samp = samples),
		#reference + ".mmi",
		#expand( "sorted/{samp}_sorted.bam", samp = samples)
#		#expand("{outFolder}{samples}_sorted.bam", samples = samples, outFolder = outFolder)


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

# collapse redundant reads
rule cupcake_collapse:
	input:
		fasta =   "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta",
		sam_sorted =  "{sample}/minimap/{sample}.hq.sorted.sam",
		cluster_report =  "{sample}/isoseq3-cluster/{sample}.polished.cluster_report.csv"
	output:
		 "{sample}/cupcake/{sample}.cupcake.collapsed.gff",
		 "{sample}/cupcake/{sample}.cupcake.collapsed.rep.fa",
		 "{sample}/cupcake/{sample}.cupcake.collapsed.group.txt",
		 "{sample}/cupcake/{sample}.cupcake.collapsed.read_stat.txt",
		 "{sample}/cupcake/{sample}.cupcake.collapsed.abundance.txt"
	params:
		prefix = "{sample}/cupcake/{sample}.cupcake"
	shell:
		"collapse_isoforms_by_sam.py --input {input.fasta} "
   		"-s {input.sam_sorted} --dun-merge-5-shorter -o {params.prefix};"
		"get_abundance_post_collapse.py {params.prefix}.collapsed {input.cluster_report}"

# assemble minimap-aligned reads into transcripts - alternative to cupcake_collapse
rule stringtie:
	input:
		bam =  "{sample}/minimap/{sample}_sorted.bam"
	params:
		gtf = referenceGTF,
		stringtiePath = "/sc/orga/projects/ad-omics/data/software/stringtie"
	output:
		gtf = "{sample}/stringtie/{sample}.stringtie.collapsed.gtf",
		gff = "{sample}/stringtie/{sample}.stringtie.collapsed.gff"
	shell:
		"{params.stringtiePath}/stringtie -G {params.gtf} -L -o {output.gtf} {input.bam};"
		"gffread -E {output.gtf} -o {output.gff}"
 

rule collapseAnnotation:
	input: referenceGTF 
	output: referenceGTF + ".genes"
	params: script = "scripts/collapse_annotation.py"
	shell: "/sc/orga/work/humphj04/conda/envs/isoseq-pipeline/bin/python {params.script} {input} {output}"

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
rule SQANTI:
	input:
		#fasta = "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta"
		gff = "{sample}/{method}/{sample}.{method}.collapsed.gff"
	output:
		report = "{sample}/SQANTI2/{sample}.{method}_classification.txt"
	params:
		sample = "{sample}.{method}",
		outDir = "{sample}/SQANTI2/",
		python = "/sc/orga/work/$USER/conda/envs/isoseq-pipeline/bin/python",
		sqantiPath= "/sc/orga/projects/ad-omics/data/software/SQANTI2",
		nCores = 4,
		#abundance = "{sample}/cupcake/{sample}.{method}.abundance.txt",
		gtf = referenceGTF,
		genome = referenceFa + ".fa",
		intropolis = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified",
		cage = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/hg38.cage_peak_phase1and2combined_coord.bed",
		polya = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/human.polyA.list.txt"
	shell:
		#"export PATH=/sc/orga/projects/ad-omics/data/software/UCSC/:$PATH;"
		#"module unload gcc;ml R/3.6.0; "
		"export PYTHONPATH=$PYTHONPATH:/hpc/users/humphj04/pipelines/cDNA_Cupcake/sequence/;"
		"{params.python} {params.sqantiPath}/sqanti_qc2.py -t {params.nCores} --aligner_choice=minimap2"
		" --dir {params.outDir} "
		" --out {params.sample} "
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
		python = "/sc/orga/work/$USER/conda/envs/isoseq-pipeline/bin/python",
                sqantiPath= "/sc/orga/projects/ad-omics/data/software/SQANTI2"
	shell:
		"{params.python} {params.sqantiPath}/sqanti_filter2.py "
		" --faa {input.faa} " #--sam {input.sam} "
		" {input.classification} {input.fasta} {input.gtf} " 	


