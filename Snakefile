
import pandas as pd

metadata = config['metadata']

metaDF = pd.read_csv(metadata, sep = '\t')
samples = metaDF['sample']

#print(samples)

dataCode = config['dataCode']
outFolder = "results/" + dataCode + "/"

fastqFolder = config['fastqFolder']

referenceFa = config['referenceFa']
referenceGTF = config['referenceGTF']

rule all:
	input:
		expand(outFolder + "SQANTI2/{sample}/{sample}.classification.txt", sample = samples),
		#outFolder + "multiqc/multiqc_report.html",
		#expand(outFolder + "rnaseqc/{samp}.metrics.tsv", samp = samples),
		#reference + ".mmi",
		#expand(outFolder + "sorted/{samp}_sorted.bam", samp = samples)
#		#expand("{outFolder}{samples}_sorted.bam", samples = samples, outFolder = outFolder)
	
rule minimapIndex:
	input: referenceFa + ".fa"
	output: referenceFa + ".mmi"
	shell:
		"minimap2 -d {output} {input}" 
rule minimap:
	input: 
		fastq = fastqFolder + "{sample}.fastq",
		ref = referenceFa + ".fa",
		index = referenceFa + ".mmi"
	params: config['minimapParams']
	output: 
		sam = outFolder + "aligned/{sample}.sam",
	shell:
		"minimap2 {params} {input.index} {input.fastq} > {output.sam};"

rule samtools:
	input: outFolder + "aligned/{samples}.sam"
	output:
		bam = outFolder + "sorted/{samples}_sorted.bam",
		bai = outFolder + "sorted/{samples}_sorted.bam.bai"
	shell:
		"samtools view -bh {input} | samtools sort > {output.bam}; "
		"samtools index {output.bam}"

rule collapseAnnotation:
	input: referenceGTF 
	output: referenceGTF + ".genes.gtf"
	params: script = "scripts/collapse_annotation.py"
	shell: "/sc/orga/work/humphj04/conda/envs/isoseq-pipeline/bin/python {params.script} {input} {output}"

rule rnaseqc:
	input:
		geneGTF = referenceGTF + ".genes.gtf",
		bam = outFolder + "sorted/{samples}_sorted.bam"
	params:
		out = outFolder + "rnaseqc/"
	output:
		outFolder + "rnaseqc/{samples}.metrics.tsv"
	shell:
		"ml rnaseqc;"
		"rnaseqc {input.geneGTF} {input.bam} {params.out} "
		" --sample={wildcards.samples} "
		" --unpaired --coverage --verbose --mapping-quality 0 --base-mismatch=1000 --detection-threshold=1"

rule samtoolsQC:
	input:
		bam = outFolder + "sorted/{samples}_sorted.bam"
	output:
		flagstat = outFolder + "qc/{samples}.flagstat.txt",
		idxstat = outFolder + "qc/{samples}.idxstat.txt"
	shell:	
		"samtools flagstat {input.bam} > {output.flagstat};"
		"samtools idxstats {input.bam} > {output.idxstat} "

rule multiQC:
	input:
                expand(outFolder + "rnaseqc/{samp}.metrics.tsv", samp = samples),
                expand(outFolder +"qc/{samp}.flagstat.txt", samp = samples),
                expand(outFolder +"qc/{samp}.idxstat.txt", samp = samples)
	output:
		outFolder + "multiqc/multiqc_report.html"
	shell:
		"export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8;"
		"multiqc -f --outdir {outFolder}multiqc/ {outFolder}" 

rule CollapseIsoforms:
	input:
		fastq = fastqFolder + "{sample}.fastq",
		sam = outFolder + "aligned/{sample}.sam"
	output:
		outFolder + "collapsed/{sample}.collapsed.gff"
	params:
		out = outFolder + "collapsed/{sample}",
		samSorted = outFolder + "aligned/{sample}.sorted.sam",
		cupcakePath = "/hpc/users/humphj04/pipelines/cDNA_Cupcake/build/scripts-3.6"
	shell:
		"sort -k 3,3 -k 4,4n {input.sam} > {params.samSorted}; " 
		"{params.cupcakePath}/collapse_isoforms_by_sam.py "
		" --input {input.fastq} --fq "
   		"-s {params.samSorted} --dun-merge-5-shorter -o {params.out}"

rule SQANTI:
	input:
		fastq = outFolder + "collapsed/{sample}.collapsed.rep.fq"
	output:
		report = outFolder + "SQANTI2/{sample}/{sample}.classification.txt"
	params:
		sample = "{sample}",
		outDir = outFolder + "SQANTI2/{sample}/",
		python = "/sc/orga/work/$USER/conda/envs/isoseq-pipeline/bin/python",
		sqantiPath= "/sc/orga/projects/ad-omics/data/software/SQANTI2",
		nCores = 4,
		gtf = referenceGTF + ".genes.gtf",
		genome = referenceFa + ".fa",
		intropolis = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified",
		cage = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/hg38.cage_peak_phase1and2combined_coord.bed",
		polya = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQANTI2/human.polyA.list.txt"
	shell:
		"export PATH=/sc/orga/projects/ad-omics/data/software/UCSC/:$PATH;"
		"ml R/3.6.0; "
		"export PYTHONPATH=$PYTHONPATH:/hpc/users/humphj04/pipelines/cDNA_Cupcake/sequence/;"
		"{params.python} {params.sqantiPath}/sqanti_qc2.py -t {params.nCores} "
		" --dir {params.outDir} "
		" --out {params.outDir}/{params.sample} "
		" --cage_peak {params.cage} --polyA_motif_list {params.polya} -c {params.intropolis}"
		" {input.fastq} {params.gtf} {params.genome} "
