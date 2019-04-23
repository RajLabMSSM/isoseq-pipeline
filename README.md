# Pipeline for working with Isoseq3 output files
**Jack Humphrey 2019**

Dependencies:
snakemake
minimap2
samtools
rnaseqc
bx-python
multiqc

Conda recipe coming soon

conda create -n isoseq-pipeline python=3.6 snakemake samtools rnaseqc minimap2
conda activate isoseq-pipeline
pip install bx-python
pip install multiqc
