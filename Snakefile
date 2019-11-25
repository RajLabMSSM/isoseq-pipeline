
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
		expand(fastqFolder + "{sample}_transcripts.classification.txt", sample = samples),
		outFolder + "multiqc/multiqc_report.html",
		expand(outFolder + "rnaseqc/{samp}.metrics.tsv", samp = samples),
		#reference + ".mmi",
		expand(outFolder + "sorted/{samp}_sorted.bam", samp = samples)
#		#expand("{outFolder}{samples}_sorted.bam", samples = samples, outFolder = outFolder)
	
rule minimapIndex:
	input: referenceFa + ".fa"
	output: referenceFa + ".mmi"
	shell:
		"minimap2 -d {output} {input}" 
rule minimap:
	input: 
		fastq = fastqFolder + "{sample}.fastq.gz",
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

rule SQUANTI:
	input:
		fastq = fastqFolder + "{sample}_transcripts.fastq"
	output:
		report = fastqFolder + "{sample}_transcripts.classification.txt"
	params:
		nCores = 4,
		gtf = referenceGTF + ".genes.gtf",
		genome = referenceFa + ".fa",
		intropolis = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQUANTI2/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz",
		cage = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQUANTI2/hg38.cage_peak_phase1and2combined_coord.bed.gz",
		polya = "/sc/orga/projects/ad-omics/data/references/hg38_reference/SQUANTI2/human.polyA.list.txt"
	shell:
		"export PYTHONPATH=$PYTHONPATH:/hpc/users/humphj04/pipelines/cDNA_Cupcake/sequence;"
		"python squanti_qc2.py -t {params.nCores} "
		" --cage_peak {params.cage} --polyA_motif_list {params.polya} -c {params.intropolis}"
		" {input.fastq} {params.gtf} {params.genome} "
