# Pipeline for working with Isoseq3 output files
**Jack Humphrey 2019**

## Dependencies:
- snakemake
- minimap2
- samtools
- rnaseqc
- bx-python
- multiqc
- bcbiogff
- gffread
- biopython
- [cDNA_cupcake](https://github.com/Magdoll/cDNA_Cupcake)
- [SQUANTI2](https://github.com/Magdoll/SQANTI2)
- [gtfToGenePred](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
## Conda recipe

```
conda create -c bioconda -n isoseq-pipeline python=3.7 snakemake samtools=1.9 minimap2
conda install -n isoseq-pipeline psutil biopython
conda install -n isoseq-pipeline -c bioconda isoseq3=3.2 pbccs=4.0
conda install -n isoseq-pipeline -c bioconda bcbiogff gffread lima pbcoretools bamtools pysam ucsc-gtftogenepred openssl=1.0 pbbam

conda activate isoseq-pipeline
pip install multiqc
pip install bx-python

# install cupcake_cDNA - Liz Tseng's code
git clone git@github.com:Magdoll/cDNA_Cupcake.git
cd cDNA_Cupcake
python setup.py build
python setup.py install --prefix=<where your conda environment is installed>
cd ..

# clone SQANTI2
git clone git@github.com:Magdoll/SQANTI2.git

# if doesn't install via conda then download UCSC tool and put in PATH
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred scripts/
chmod +x scripts/gtfToGenePred
echo "export PATH="$PWD/scripts/:\$PATH >> ~/.bashrc"
```

## Symlink reference files

```
mkdir reference
# if you want hg19:
ln -s /hpc/users/humphj04/GENCODE/gencode.v30lift37.annotation.gtf reference/gencode.v30lift37.annotation.gtf
ln -s /sc/orga/projects/ad-omics/ricardo/Data/1000G_phase1/human_g1k_v37.fasta reference/human_g1k_v37.fa

# hg38 - SQANTI only supports this
ln -s /sc/orga/projects/ad-omics/data/references/hg38_reference/hg38.fa reference/hg38.fa
cp /sc/orga/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.annotation.gtf.gz reference/gencode.v30.annotation.gtf.gz
gunzip reference/gencode.v30.annotation.gtf.gz
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

* Installing isoseq on conda:
   - https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda

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


* [IsoAnnotLite - annotate novel isoforms from PacBio reads](http://tappas.org/what-if-i-come-from-pacbio/)

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


