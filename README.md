# Pipeline for working with Isoseq3 output files
**Jack Humphrey 2019**

## Dependencies:
- snakemake
- minimap2
- samtools
- rnaseqc
- bx-python
- multiqc

## Conda recipe

```
conda create -n isoseq-pipeline python=3.6 snakemake samtools rnaseqc minimap2
conda activate isoseq-pipeline
pip install bx-python
pip install multiqc
```

## Resources

- How isoseq works:
https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.1.md

- What to do with the output of Isoseq:
https://github.com/Magdoll/cDNA_Cupcake/wiki/Best-practice-for-aligning-Iso-Seq-to-reference-genome:-minimap2,-GMAP,-STAR,-BLAT
https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/What-to-do-after-Iso-Seq-Cluster%3F

- cDNA_Cupcake - collapse long reads into unique transcripts
https://github.com/Magdoll/cDNA_Cupcake

- SQANTI2 - identify isoforms from aligned long reads:
https://github.com/Magdoll/SQANTI2

- Cogent - reconstruct coding genome from long reads without a reference genome
https://github.com/Magdoll/Cogent

- bioawk - awk + built in parsing of bio data (SAM, GTF, BED, etc)
https://github.com/lh3/bioawk

- Lists of non-polyadenylated genes for verify polyA annealing data
https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-2-r16


