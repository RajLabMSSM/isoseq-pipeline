
import pandas as pd

metadata = config['metadata']

metaDF = pd.read_csv(metadata, sep = '\t')
samples = metaDF['sample']

#print(samples)

dataCode = config['dataCode']

referenceFa = config['referenceFa']
referenceGTF = config['referenceGTF']

primers = "reference/NEB_primers_01_2019.fa"

rule all:
	input:
		expand( "{sample}/cupcake/{sample}.hq.collapsed.abundance.txt", sample = samples)
	#	expand(fastqFolder + "{sample}.classification.txt", sample = samples),
		# "multiqc/multiqc_report.html",
		#expand( "rnaseqc/{samp}.metrics.tsv", samp = samples),
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
		"lima --isoseq --different --min-passes 1 --split-bam-named --dump-clips --dump-removed -j 0 {input.file} {primers} {output}"

# trimming of polya tails and removal of concatemers
rule isoseq_refine:
	input:
		 "{sample}/isoseq3-lima/{sample}.fl.bam"
	output:
		 "{sample}/isoseq3-refine/{sample}.flnc.bam"
	shell:
		"isoseq3 refine --require-polya {input} {primers} {output}"

rule isoseq_cluster:
	input:
		 "{sample}/isoseq3-refine/{sample}.flnc.bam"
	output:
		bam =  "{sample}/isoseq3-cluster/{sample}.polished.bam",
		fastq =  "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta.gz",
		report =  "{sample}/isoseq3-cluster/{sample}.polished.cluster_report.csv"
	shell:
		"isoseq3 cluster --verbose --use-qvs -j 0 {input} {output.bam}"

rule minimapIndex:
	input: referenceFa + ".fa"
	output: referenceFa + ".mmi"
	shell:
		"minimap2 -d {output} {input}" 
rule minimap:
	input: 
		fastq =  "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta.gz",
		ref = referenceFa + ".fa",
		index = referenceFa + ".mmi"
	params: 
		"-ax splice -t 30 -uf --secondary=no -C5"
		#config['minimapParams']
	output: 
		sam =  "{sample}/minimap/{sample}.hq.sam",
		sam_sorted =  "{sample}/minimap/{sample}.hq.sorted.sam"
	shell:
		"minimap2 {params} {input.index} {input.fastq} > {output.sam};"
		"sort -k 3,3 -k 4,4n {output.sam} > {output.sorted_sam}"
rule samtools:
	input:  "{sample}/minimap/{sample}.hq.sam"
	output:
		bam =  "{sample}/minimap/{samples}_sorted.bam",
		bai =  "{sample}/minimap/{samples}_sorted.bam.bai",
		flagstat =  "{sample}/qc/{samples}.flagstat.txt",
                idxstat =  "{sample}/qc/{samples}.idxstat.txt"
	shell:
		"samtools view -bh {input} | samtools sort > {output.bam}; "
		"samtools index {output.bam}"
                "samtools flagstat {output.bam} > {output.flagstat};"
                "samtools idxstats {output.bam} > {output.idxstat} "

rule cupcake_collapse:
	input:
		fastq =   "{sample}/isoseq3-cluster/{sample}.polished.hq.fasta.gz"
		sam_sorted =  "{sample}/minimap/{sample}.hq.sorted.sam"
		cluster_report =  "{sample}/isoseq3-cluster/{sample}.polished.cluster_report.csv"
	output:
		 "{sample}/cupcake/{sample}.hq.collapsed.gff",
		 "{sample}/cupcake/{sample}.hq.collapsed.rep.fq",
		 "{sample}/cupcake/{sample}.hq.collapsed.group.txt",
		 "{sample}/cupcake/{sample}.hq.collapsed.read_stat.txt",
		 "{sample}/cupcake/{sample}.hq.collapsed.abundance.txt"
	params:
		prefix = "{sample}.hq"
	shell:
		"collapse_isoforms_by_sam.py --input {input.fastq} --fq "
   		"-s {input.sam_sorted} --dun-merge-5-shorter -o {params.prefix};"
		"get_abundance_post_collapse.py {params.prefix}.collapsed {input.cluster_report}"

rule collapseAnnotation:
	input: referenceGTF 
	output: referenceGTF + ".genes.gtf"
	params: script = "scripts/collapse_annotation.py"
	shell: "/sc/orga/work/humphj04/conda/envs/isoseq-pipeline/bin/python {params.script} {input} {output}"

rule rnaseqc:
	input:
		geneGTF = referenceGTF + ".genes.gtf",
		bam =  "sorted/{samples}_sorted.bam"
	params:
		out =  "{sample}/qc/"
	output:
		 "{sample}/qc/{samples}.metrics.tsv"
	shell:
		"ml rnaseqc;"
		"rnaseqc {input.geneGTF} {input.bam} {params.out} "
		" --sample={wildcards.samples} "
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
		fastq = "{sample}/cupcake/{sample}.hq.collapsed.rep.fq"
	output:
		report = fastqFolder + "{sample}.classification.txt"
	params:
		#python = "/sc/orga/work/$USER/conda/envs/isoseq-pipeline/bin/python",
		sqantiPath= "/sc/orga/projects/ad-omics/data/software/SQANTI2",
		nCores = 4,
		gtf = referenceGTF,
		genome = referenceFa + ".fa",
		intropolis = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz",
		cage = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/hg38.cage_peak_phase1and2combined_coord.bed.gz",
		polya = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/human.polyA.list.txt"
	shell:
	`	"ml R/3.6.0; "
		"export PYTHONPATH=$PYTHONPATH:/hpc/users/humphj04/pipelines/cDNA_Cupcake/sequence/;"
		"python {params.sqantiPath}/sqanti_qc2.py -t {params.nCores} --aligner_choice=minimap2"
		" --cage_peak {params.cage} --polyA_motif_list {params.polya} -c {params.intropolis}"
		" {input.fastq} {params.gtf} {params.genome} "
