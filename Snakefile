
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
		outFolder + "multiqc/multiqc_report.html",
		expand(outFolder + "rnaseqc/{samp}.metrics.tsv", samp = samples),
		#reference + ".mmi",
		expand(outFolder + "sorted/{samp}_transcripts_sorted.bam", samp = samples)
#		#expand("{outFolder}{samples}_transcripts_sorted.bam", samples = samples, outFolder = outFolder)
	
rule minimapIndex:
	input: referenceFa + ".fa"
	output: referenceFa + ".mmi"
	shell:
		"minimap2 -d {output} {input}" 
rule minimap:
	input: 
		fastq = fastqFolder + "{sample}_transcripts.fastq",
		ref = referenceFa + ".fa",
		index = referenceFa + ".mmi"
	params: config['minimapParams']
	output: 
		sam = outFolder + "aligned/{sample}_transcripts.sam",
	shell:
		"minimap2 {params} {input.index} {input.fastq} > {output.sam};"

rule samtools:
	input: outFolder + "aligned/{samples}_transcripts.sam"
	output:
		bam = outFolder + "sorted/{samples}_transcripts_sorted.bam",
		bai = outFolder + "sorted/{samples}_transcripts_sorted.bam.bai"
	shell:
		"samtools view -bh {input} | samtools sort > {output.bam}; "
		"samtools index {output.bam}"

rule collapseAnnotation:
	input: referenceGTF 
	output: referenceGTF + ".genes.gtf"
	params: script = "scripts/collapse_annotation.py"
	shell: "python {params.script} {input} {output}"

rule rnaseqc:
	input:
		geneGTF = referenceGTF + ".genes.gtf",
		bam = outFolder + "sorted/{samples}_transcripts_sorted.bam"
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
		bam = outFolder + "sorted/{samples}_transcripts_sorted.bam"
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
