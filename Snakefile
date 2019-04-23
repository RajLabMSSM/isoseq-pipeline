
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
		expand("qc/{samp}.metrics.tsv", samp = samples),  
		#reference + ".mmi",
		expand("sorted/{samp}_hq_transcripts_sorted.bam", samp = samples)
#		expand("{outFolder}{samples}_hq_transcripts_sorted.bam", samples = samples, outFolder = outFolder)
	
rule minimapIndex:
	input: referenceFa + ".fa"
	output: referenceFa + ".mmi"
	shell:
		"minimap2 -d {output} {input}" 
rule minimap:
	input: 
		fastq = "fastq/{sample}_hq_transcripts.fastq",
		ref = referenceFa + ".fa",
		index = referenceFa + ".mmi"
	params: " -ax splice -uf -C5 "	
	output: sam = "aligned/{sample}_hq_transcripts.sam"
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

rule rnaseqc:
	input:
		geneGTF = referenceGTF + ".genes.gtf",
		bam = "sorted/{samples}_hq_transcripts_sorted.bam"
	output:
		"qc/{sample}.metrics.tsv"
	shell:
		"rnaseqc {input.geneGTF} {input.bam} --coverage --unpaired {out}"
