
import pandas as pd

metadata = config['metadata']

metaDF = pd.read_csv(metadata, sep = '\t')
samples = metaDF['sample']

print(samples)

referenceFa = config['referenceFa']
referenceGTF = config['referenceGTF']

outFolder = config['outFolder']
inFolder = config['inFolder']
rule all:
	input:
		"multiqc/multiqc_report.html"
		#expand("qc/{samp}.metrics.tsv", samp = samples),
		#expand("qc/{samp}.flagstat.txt", samp = samples),
		#expand("qc/{samp}.idxstat.txt", samp = samples),  
		#reference + ".mmi",
		#expand("sorted/{samp}_hq_transcripts_sorted.bam", samp = samples)
#		expand("{outFolder}{samples}_hq_transcripts_sorted.bam", samples = samples, outFolder = outFolder)
	
rule minimapIndex:
	input: referenceFa + ".fa"
	output: referenceFa + ".mmi"
	shell:
		"minimap2 -d {output} {input}" 
rule minimap:
	input: 
		fastq = "fastq/{samples}_hq_transcripts.fastq",
		ref = referenceFa + ".fa",
		index = referenceFa + ".mmi"
	params: " -ax splice -uf -C5 "	
	output: sam = "aligned/{samples}_hq_transcripts.sam"
	shell:
		"minimap2 {params} {input.ref} {input.fastq} > {output.sam}"

rule samtools:
	input: "aligned/{samples}_hq_transcripts.sam"
	output:
		bam = "sorted/{samples}_hq_transcripts_sorted.bam",
		bai = "sorted/{samples}_hq_transcripts_sorted.bam.bai"
	shell:
		"samtools view -bh {input} | samtools sort > {output.bam}; "
		"samtools index {output.bam}"

rule collapseAnnotation:
	input: referenceGTF 
	output: referenceGTF + ".genes.gtf"
	params: script = "scripts/collapse_annotation.py"
	shell: "python {params.script} {input} {output}"

rule QC:
	input:
		geneGTF = referenceGTF + ".genes.gtf",
		bam = "sorted/{samples}_hq_transcripts_sorted.bam"
	output:
		rnaseqc = "qc/{samples}.metrics.tsv",
		flagstat = "qc/{samples}.flagstat.txt",
		idxstat = "qc/{samples}.idxstat.txt"
	shell:	
		"samtools flagstat {input.bam} > {output.flagstat};"
		"samtools idxstats {input.bam} > {output.idxstat};"
		"rnaseqc {input.geneGTF} {input.bam} -s {wildcards.samples} --coverage qc/"

rule multiQC:
	input:
                expand("qc/{samp}.metrics.tsv", samp = samples),
                expand("qc/{samp}.flagstat.txt", samp = samples),
                expand("qc/{samp}.idxstat.txt", samp = samples)
	output:
		"multiqc/multiqc_report.html"
	shell:
		"multiqc -f --outdir multiqc/ ." 
