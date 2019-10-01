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
conda create -n isoseq-pipeline python=3.6 snakemake samtools minimap2
conda activate isoseq-pipeline
pip install multiqc
pip install bx-python
```

## Symlink reference files

```
mkdir reference
ln -s /hpc/users/humphj04/GENCODE/gencode.v30lift37.annotation.gtf reference/gencode.v30lift37.annotation.gtf
ln -s /sc/orga/projects/ad-omics/ricardo/Data/1000G_phase1/human_g1k_v37.fasta reference/human_g1k_v37.fa
```

## Running on test data

```
snakemake --configfile config-files/test_config.yaml
```

## Outline of pipeline

1. Align reads to reference genome using minimap2

2. Index BAM files, run flagstat and idxstat for QC

2. Get more QC with RNASeqQC

3. Pool QC metrics together with multiQC

Future:

* Collapse reads with TAMA / SQUANTI2 / CupCake

* Infer transcript function with IsoAnnot - unreleased yet

* Create Kallisto reference for short reads


## Resources

* How isoseq works:  
    - https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.1.md

* What to do with the output of Isoseq:
    * https://github.com/Magdoll/cDNA_Cupcake/wiki/Best-practice-for-aligning-Iso-Seq-to-reference-genome:-minimap2,-GMAP,-STAR,-BLAT
    * https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/What-to-do-after-Iso-Seq-Cluster%3F


## Tools

* [TAMA - Transcriptome Annotation by Modular Algorithms](https://github.com/GenomeRIK/tama/wiki)
- collapse aligned reads to find unique transcripts
- merge multiple transcriptomes together

* [cDNA_Cupcake - collapse long reads into unique transcripts](https://github.com/Magdoll/cDNA_Cupcake)

*  [SQANTI2 - identify isoforms from aligned long reads](https://github.com/Magdoll/SQANTI2)

* [isoAnnot - database of isoform functions - not available yet ]()

* [tappAS -Your application to understand the functional implications of alternative splicing](http://tappas.org/)

* [Cogent - reconstruct coding genome from long reads without a reference genome](https://github.com/Magdoll/Cogent)


## Papers

[Kuo et al - Illuminating the dark side of the human transcriptome](https://www.biorxiv.org/content/biorxiv/early/2019/09/24/780015.full.pdf)

Analyses Universal Human Reference RNA Iso-seq sample, comparing 5 different pipelines.
Discusses multiple sources of error/contamination that can be misinterpreted as novel transcripts:
* genomic fragments
* internal priming of pre-mRNA
* 5' degradation leading to novel shortened isoforms
* wobble of alignment between read and reference due to read errors
* chimeric reads caused by polishing of long reads by short reads from homologous but different transcripts

![Fig 1 Kuo et al](https://www.biorxiv.org/content/biorxiv/early/2019/09/24/780015/F1.large.jpg?width=800&height=600&carousel=1)




## Misc

* bioawk - awk + built in parsing of bio data (SAM, GTF, BED, etc)
    * https://github.com/lh3/bioawk

* Lists of non-polyadenylated genes for verify polyA annealing data
    * https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-2-r16


