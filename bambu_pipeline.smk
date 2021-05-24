import pandas as pd
import os

# bambu pipeline

# uses r_test environment where bambu is installed on R 4.0.1
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

# bambu specific params
matching_strategy = "precise"
model_strategy = "default_ccs"
bambu_cores = "8"

sqanti_threads = "8"

prefix = out_folder + "bambu/" + data_code
miss_prefix = out_folder + "bambu/missingness/" + data_code
sqanti_prefix = out_folder + "bambu/SQANTI/" +  data_code + "_bambu"  
filter_prefix = out_folder + "bambu/filter/" + data_code
cpat_prefix = out_folder + "bambu/CPAT/" + data_code
cpat_folder = "/sc/arion/projects/ad-omics/data/references/CPAT"

td_prefix = out_folder + "bambu/TransDecode/" + data_code
td_out = out_folder + "bambu/TransDecode/" 
td_folder = "/sc/arion/projects/ad-omics/data/software/TransDecoder-v5.5.0/"

suppa_prefix = out_folder + "bambu/SUPPA/" + data_code

junctionFolder = "/sc/arion/projects/als-omics/microglia_isoseq/short_read/junctions"

rule all:
    input:
        prefix + "_extended_annotations.gtf",
        sqanti_prefix + "_classification.txt",
        #filter_prefix + "_miss_sqanti.sorted.gtf.gz",
        expand( "{sample}.sorted.gtf.gz.tbi", sample = [filter_prefix + "_miss_sqanti.cds.gtf", sqanti_prefix + "_corrected.gtf.cds.gff"]), #td_prefix + ".transdecoder.genome.gff3" ] ),
        cpat_prefix + ".ORF_seqs.fa",
        expand(suppa_prefix + ".events_{event_type}_strict.ioe", event_type = ["SE", "MX","RI","AF", "AL", "A3", "A5"]),
        suppa_prefix + ".all_suppa_events.ioe",
        suppa_prefix + "_events.psi.gz",
        td_prefix + ".transdecoder.genome.gff3"

rule create_annotation:
    input: 
        gtf = ref_gtf
    output:
        rdata = out_folder + "bambu/bambu_annotation.RData"
    params:
        script = "scripts/bambu_annotation.R"
    shell:
        "Rscript {params.script} -i {input.gtf} -o {output.rdata}"

rule run_bambu:
    input:
        bams = expand( out_folder + "{sample}/pbmm2/{sample}.aligned.md.bam", sample = samples),
        anno = out_folder + "bambu/bambu_annotation.RData"
    params:
        script = "scripts/bambu_run.R"
    output:
        prefix + "_extended_annotations.gtf",
        prefix + "_counts_transcript.txt",
        prefix + "_bambu.RData"
    shell:
        "Rscript {params.script} --cores {bambu_cores} --fasta {ref_fasta} --anno {input.anno} --prefix {prefix} {input.bams}"

## FILTER MISSINGNESS
# remove all transcripts with greater than X% missingness
# currently - present in at least 2 samples
# remove monoexonic transcripts to reduce overhead for SQANTI
# input must be simple GFF - long GTF lines break rtracklayer
rule filter_missingness:
    input:
        counts = prefix + "_counts_transcript.txt",
        gtf = prefix + "_extended_annotations.gtf"
    output:
        counts = miss_prefix + "_miss_counts.csv",
        tpm = miss_prefix + "_miss_tpm.csv",
        gtf = miss_prefix + "_miss.gtf"
    params:
        script = "scripts/filter_missingness.R",
        prefix = miss_prefix,
        min_samples = 5, # eventually put in config
        min_reads = 1 # 1 TPM
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} --matrix {input.counts} --gff {input.gtf} --prefix {params.prefix} --min_samples {params.min_samples} --min_reads {params.min_reads} --remove_monoexons --tpm"

# run SQANTI using filtered GTF
rule SQANTI:
    input:
         gtf = miss_prefix + "_miss.gtf",
         abundance = miss_prefix + "_miss_counts.csv"
    output:
         out = sqanti_prefix + "_classification.txt",
         gtf = sqanti_prefix  + "_corrected.gtf",
         fasta = sqanti_prefix + "_corrected.fasta",
         gff = sqanti_prefix + "_corrected.gtf.cds.gff"
    params:
        sample = data_code + "_bambu" ,
        outDir = out_folder + "bambu/SQANTI/",
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
        counts = miss_prefix + "_miss_counts.csv",
        tpm = miss_prefix + "_miss_tpm.csv",
        gff = sqanti_prefix + "_corrected.gtf.cds.gff",
        sqanti = sqanti_prefix + "_classification.txt",
        fasta = sqanti_prefix + "_corrected.fasta"
    output:
        counts = filter_prefix + "_miss_sqanti_counts.csv",
        tpm = filter_prefix + "_miss_sqanti_tpm.csv",
        gff = filter_prefix + "_miss_sqanti.cds.gtf",
        sqanti = filter_prefix + "_miss_sqanti_classification.tsv",
        fasta =  filter_prefix + "_miss_sqanti.fasta"
    params:
        script = "scripts/filter_sqanti.R"
    shell:
        "Rscript {params.script} --input {miss_prefix} --output {filter_prefix} --sqanti {input.sqanti} --fasta {input.fasta} --gff {input.gff}" 

# sort and tabix index final GFF
rule indexGFF:
    input:
        "{gff}"
    output:
        gtf = "{gff}.sorted.gtf.gz",
        index = "{gff}.sorted.gtf.gz.tbi"
    params:
        gff3sort = "/sc/arion/projects/ad-omics/data/software/gff3sort/gff3sort.pl"
    shell:
        "ml tabix;"
        "{params.gff3sort} {input} | bgzip > {output.gtf};"
        "tabix {output.gtf} "

# run CPAT to predict ORFs
rule cpat:
    input:
        fasta = filter_prefix + "_miss_sqanti.fasta"
    output:
        cpat_prefix + ".ORF_seqs.fa"
    shell:
        "conda activate isoseq-pipeline;"
        "cpat.py -x {cpat_folder}/Human_Hexamer.tsv -d {cpat_folder}/Human_logitModel.RData --top-orf=5 -g {input.fasta} -o {cpat_prefix}"

# also TransDecoder
rule TransDecoder:
    input:
        fasta = filter_prefix + "_miss_sqanti.fasta",
        gtf = filter_prefix + "_miss_sqanti.cds.gtf"
    output:
        td_gff = td_out + "longest_orfs.gff3",
        gff = td_prefix + "_miss_sqanti.gff3",
        cds_gff = td_prefix + ".transdecoder.genome.gff3"
    shell:
        "conda activate isoseq-pipeline;"
        "{td_folder}/util/gtf_to_alignment_gff3.pl {input.gtf} > {output.gff};"
        "{td_folder}/TransDecoder.LongOrfs -S -t {input.fasta};"
        "mv *.transdecoder_dir/* {td_out} ; "
        "{td_folder}/util/cdna_alignment_orf_to_genome_orf.pl "
        "   {output.td_gff} {output.gff} {input.fasta} > {output.cds_gff};"
        "rm -r *.transdecoder_dir*; rm pipeliner*;"

# get AS events from GTF 
rule SUPPA_events:
    input:
        filter_prefix + "_miss_sqanti.cds.gtf"
    output:
        events = expand(suppa_prefix + ".events_{event_type}_strict.ioe", event_type = ["SE", "MX","RI","AF", "AL", "A3", "A5"]),
        total = suppa_prefix + ".all_suppa_events.ioe"
    params:
        prefix = suppa_prefix + ".events"
    shell:
        "conda activate isoseq-pipeline;"
        "suppa.py generateEvents -i {input} -o {params.prefix} -e SE SS MX RI FL -f ioe --pool-genes;"
        " awk 'FNR==1 && NR!=1 {{ while (/^<header>/) getline; }} 1 {{print}}' {output.events} > {output.total}"

# convert TPM
rule SUPPA_convert_TPM:
    input:
        filter_prefix + "_miss_sqanti_tpm.csv"
    output:
        suppa_prefix + "_suppa_tpm.tsv"
    params:
        script = "scripts/convert_tpm_for_suppa.R"
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} -i {input} -o {output}"

# quantify PSI 
rule SUPPA_quantify:
    input:
        events = suppa_prefix + ".all_suppa_events.ioe",
        tpm = suppa_prefix + "_suppa_tpm.tsv"
    output:
        suppa_prefix + "_events.psi.gz"
    params:
        prefix = suppa_prefix + "_events"
    shell:
        "conda activate isoseq-pipeline;"
        "suppa.py psiPerEvent -i {input.events} -e {input.tpm} -o {params.prefix};"
        "gzip {suppa_prefix}_events.psi"
        

