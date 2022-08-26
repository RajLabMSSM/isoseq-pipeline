import pandas as pd
import os

# stringtie pipeline
R_VERSION = "R/4.0.3"
shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate snakemake; ml R/4.0.3;")

ref_fasta = config["ref_genome"] + ".fa"
ref_gtf = config["ref_gtf"]
metadata = config["metadata"]
data_code = config["data_code"]
out_folder = config["out_folder"]

# read in metadata
meta_df = pd.read_excel(metadata)
samples = meta_df['sample']

metadata_dict = meta_df.set_index("sample").T.to_dict()
#isoquant = "/sc/arion/projects/ad-omics/data/software/IsoQuant/isoquant.py"

#stringtie = "/sc/arion/projects/ad-omics/data/software/stringtie-2.2.1.Linux_x86_64/stringtie"
stringtie = "/sc/arion/projects/ad-omics/data/software/stringtie-2.2.1.compiled/stringtie"

salmon = "/sc/arion/projects/ad-omics/data/software/salmon-1.9.0_linux_x86_64/bin/salmon"
gffread = "/hpc/packages/minerva-common/cufflinks/2.2.1/bin/gffread"
# stringtie specific params
stringtie_threads = "4"
merge_threads = "4"
# sqanti
sqanti_threads = "8"

prefix = out_folder + "stringtie/" + data_code
miss_prefix = out_folder + "stringtie/filter1/" + data_code
sqanti_prefix = out_folder + "stringtie/SQANTI/" +  data_code
filter_prefix = out_folder + "stringtie/filter/" + data_code
cpat_prefix = out_folder + "stringtie/CPAT/" + data_code
cpat_folder = "/sc/arion/projects/ad-omics/data/references/CPAT"

td_prefix = out_folder + "stringtie/TransDecode/" + data_code
td_out = out_folder + "stringtie/TransDecode/" 
td_folder = "/sc/arion/projects/ad-omics/data/software/TransDecoder-v5.5.0/"

suppa_prefix = out_folder + "stringtie/SUPPA/" + data_code

junctionFolder = "/sc/arion/projects/als-omics/microglia_isoseq/short_read_junctions/junctions/"
rule all:
    input:
        filter_prefix + "_filter_sqanti.cds.sorted.gtf.gz",
        miss_prefix + "_filter_fpkm.csv",
        miss_prefix + "_filter.gtf",
        prefix + "_all_samples_merged_stringtie.gtf",
        expand( out_folder + "{sample}/stringtie/sample_{sample}/t_data.ctab", sample = samples)
        #prefix + "_extended_annotations.gtf",
        #sqanti_prefix + "_classification.txt",
        #filter_prefix + "_filter_sqanti.sorted.gtf.gz",
        #expand( "{sample}.sorted.gtf.gz.tbi", sample = [filter_prefix + "_filter_sqanti.cds.gtf", sqanti_prefix + "_corrected.gtf.cds.gff"]), #td_prefix + ".transdecoder.genome.gff3" ] ),
        #cpat_prefix + ".ORF_seqs.fa",
        #expand(suppa_prefix + ".events_{event_type}_strict.ioe", event_type = ["SE", "MX","RI","AF", "AL", "A3", "A5"]),
        #suppa_prefix + ".all_suppa_events.ioe",
        #suppa_prefix + "_events.psi.gz",
        #td_prefix + ".transdecoder.genome.gff3"

rule run_stringtie:
    input:
        gtf = ref_gtf,
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output:
        gtf = out_folder + "{sample}/stringtie/{sample}.stringtie.gtf"
    shell:
        "{stringtie} -p {stringtie_threads} -o {output.gtf} -L -G {input.gtf} {input.bam}"

rule merge_stringtie:
    input:
        ref = ref_gtf,
        gtf = expand( out_folder + "{sample}/stringtie/{sample}.stringtie.gtf", sample = samples)
    output:
        prefix + "_all_samples_merged_stringtie.gtf"
    shell:
      "{stringtie} -p {merge_threads} --merge -o {output} -i -G {input.ref} {input.gtf}"

rule quant_stringtie:
    input:
        gtf = prefix + "_all_samples_merged_stringtie.gtf",
        bam = out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam"
    output:
        quant = out_folder + "{sample}/stringtie/sample_{sample}/t_data.ctab"
    shell:
        "{stringtie} -o {out_folder}/{wildcards.sample}/stringtie/sample_{wildcards.sample}/quant.txt -eB -G {input.gtf} {input.bam}"

# remove monoexons
# keep all annotated found at least once
# keep novel tx if found twice
rule stringtie_filter:
    input:
        counts = expand(out_folder + "{sample}/stringtie/sample_{sample}/t_data.ctab", sample = samples),
        gtf = prefix + "_all_samples_merged_stringtie.gtf"
    output:
        counts = miss_prefix + "_filter_fpkm.csv",
        gtf = miss_prefix + "_filter.gtf"
    params:
        script = "scripts/stringtie_filter.R",
        prefix = miss_prefix,
        min_samples = 2, # eventually put in config
        min_reads = 0 # FPKM
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} --inFolder {out_folder} --gff {input.gtf} --prefix {params.prefix} --min_samples {params.min_samples} --remove_monoexons"

# run SQANTI using filtered GTF
rule SQANTI:
    input:
         gtf = miss_prefix + "_filter.gtf",
         abundance = miss_prefix + "_filter_fpkm.csv"
    output:
         out = sqanti_prefix + "_classification.txt",
         gtf = sqanti_prefix  + "_corrected.gtf",
         fasta = sqanti_prefix + "_corrected.fasta",
         gff = sqanti_prefix + "_corrected.gtf.cds.gff"
    params:
        sample = data_code,
        outDir = out_folder + "stringtie/SQANTI/",
        nCores = sqanti_threads,
        nChunks = 8,
        software= "/sc/arion/projects/ad-omics/data/software",
        #junctions = junctionArgs,
        junctions = "\'" + junctionFolder + "/*SJ.out.tab\'" ,
        gtf = ref_gtf,
        #genome = referenceFa + ".fa",
        intropolis = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified",
        cage = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/hg38.cage_peak_phase1and2combined_coord.bed",
        polya = "/sc/arion/projects/ad-omics/data/references/hg38_reference/SQANTI3/human.polyA.list.txt",
        isoAnnotGFF = "/sc/arion/projects/ad-omics/data/references/hg38_reference/RefSeq/Homo_sapiens_GRCh38_RefSeq_78.gff3"
    shell:
        "conda activate SQANTI3.env; module purge;"
        "ml R/4.0.3;"
        "export PYTHONPATH=$PYTHONPATH:{params.software}/cDNA_Cupcake/sequence;"
        "export PYTHONPATH=$PYTHONPATH:{params.software}/cDNA_Cupcake/;"
        "python {params.software}/SQANTI3/sqanti3_qc.py -t {params.nCores} "
        " --dir {params.outDir} "
        " --out {params.sample} "
        " -c {params.junctions} "
        " --cage_peak {params.cage} --polyA_motif_list {params.polya} "
        #"--skipORF " # ORF finding is slow, can skip if testing
        #"-c {params.intropolis}"
        #" --fl_count {input.abundance}"
        " --gtf {input.gtf} "
        #" --isoAnnotLite --gff3 {params.isoAnnotGFF}"
        " {params.gtf} {ref_fasta} "

## filter SQANTI
rule filter_sqanti:
    input:
        #counts = miss_prefix + "_miss_counts.csv",
        tpm = miss_prefix + "_filter_fpkm.csv",
        gff = sqanti_prefix + "_corrected.gtf.cds.gff",
        sqanti = sqanti_prefix + "_classification.txt",
        fasta = sqanti_prefix + "_corrected.fasta"
    output:
        counts = filter_prefix + "_filter_sqanti_counts.csv",
        #tpm = filter_prefix + "_filter_sqanti_fkpm.csv",
        gff = filter_prefix + "_filter_sqanti.cds.gtf",
        sqanti = filter_prefix + "_filter_sqanti_classification.tsv",
        fasta =  filter_prefix + "_filter_sqanti.fasta"
    params:
        script = "scripts/filter_sqanti.R"
    shell:
        "ml R/4.0.3; Rscript {params.script} --counts {input.tpm} --input {miss_prefix} --output {filter_prefix} --sqanti {input.sqanti} --fasta {input.fasta} --gff {input.gff}"


# sort and tabix index final GFF
rule indexGFF:
    input:
        gtf = filter_prefix + "_filter_sqanti.cds.gtf",
    output:
        gtf = filter_prefix + "_filter_sqanti.cds.sorted.gtf.gz",
        index = filter_prefix + "_filter_sqanti.cds.sorted.gtf.gz.tbi"
    params:
        gff3sort = "/sc/arion/projects/ad-omics/data/software/gff3sort/gff3sort.pl"
    shell:
        "ml tabix;"
        "{params.gff3sort} {input.gtf} | bgzip > {output.gtf};"
        "tabix {output.gtf} "


# later
rule create_fasta:
    input:
        prefix + "_all_samples_merged_stringtie.gtf"
    output:
        prefix + "_all_samples_merged_stringtie.fasta"
    shell:
        "{gffread} -w {output} -g {genome} {input}"

rule salmon_index:
    input:
        prefix + "_all_samples_merged_stringtie.fasta"
    output:
        prefix + "_all_samples_merged_stringtie.salmon"
    shell:
        "{salmon} index -t {input} -i {output}"
#
