# qc_pipeline.smk
# Long-read QC + alignment (minimap2 for ONT, pbmm2 for PacBio FLNC), selected by
# `prep`. Included by Snakefile (owns shell.prefix, config, samples, out_folder,
# groups, mmi).

minimap2_bin = "/sc/arion/projects/ad-omics/data/software/minimap2-2.31/minimap2"
genome = referenceFa
rrna   = "/sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/gencode.v30.rRNA.interval.list"

# canonical aligned BAM: both aligners write here, so downstream never branches on platform
output_bam = out_folder + "{sample}/alignment/{sample}.aligned.bam"

# minimap2 preset per prep (config can override)
if prep == "nanopore_cdna":
    minimap_preset = "splice:hq"; minimap_flags = "-uf"
if prep == "nanopore_direct":
    minimap_preset = "splice";    minimap_flags = "-uf -k14"   # direct RNA
if prep == "pacbio":
    minimap_preset = "splice:hq"; minimap_flags = "-uf"
if config.get('minimap_preset'):
    minimap_preset = config['minimap_preset']
if config.get('minimap_flags'):
    minimap_flags = config['minimap_flags']

# raw pre-alignment reads for NanoStat: ONT fastq, or PacBio unaligned FLNC bam
def raw_reads(wildcards):
    if prep == "pacbio":
        return pacbio_flnc(wildcards.sample)
    return metadata_dict[wildcards.sample]["long_read_fastq"]


# reference
if not config.get('ref_genome_index'):
    rule create_index:
        input:  genome
        output: mmi
        shell:
            "{minimap2_bin} -x {minimap_preset} -d {output} {input}"


# alignment: one aligner active per run (gated by prep), both -> canonical BAM
if prep in ("nanopore_cdna", "nanopore_direct"):

    rule nanopore_alignment:
        input:
            index = mmi,
            junc  = junctions_bed
        output:
            sam = temp(out_folder + "{sample}/alignment/{sample}.aligned.sam")
        threads: 8
        run:
            fastq_file = metadata_dict[wildcards.sample]["long_read_fastq"]
            shell(
                "{minimap2_bin} "
                "-ax {minimap_preset} {minimap_flags} "
                "--secondary=no -G 500k "
                "--junc-bed {input.junc} "
                "-t {threads} "
                "{input.index} {fastq_file} "
                "> {output.sam}"
            )

    rule sam_to_bam:
        input:
            out_folder + "{sample}/alignment/{sample}.aligned.sam"
        output:
            bam = out_folder + "{sample}/alignment/{sample}.aligned.bam",
            bai = out_folder + "{sample}/alignment/{sample}.aligned.bam.bai"
        threads: 8
        shell:
            "samtools sort -@ {threads} -m 2G -o {output.bam} {input}; "
            "samtools index {output.bam}"


if prep == "pacbio":

    # flnc resolved by pacbio_flnc() (metadata or CCS); minimap2 mmi reused as pbmm2 index
    rule pacbio_alignment:
        input:
            index = mmi,
            flnc  = lambda wc: pacbio_flnc(wc.sample)
        output:
            bam = out_folder + "{sample}/alignment/{sample}.aligned.bam",
            bai = out_folder + "{sample}/alignment/{sample}.aligned.bam.bai"
        threads: 8
        shell:
            "pbmm2 align --sort -j {threads} --sort-threads 4 -m 3G "
            "--preset=ISOSEQ --log-level INFO --unmapped "
            "{input.index} {input.flnc} | "
            "samtools calmd -b - {genome} > {output.bam}; "
            "samtools index {output.bam}"


# QC metrics
rule samtools_qc:
    input:
        bam = output_bam,
        bai = out_folder + "{sample}/alignment/{sample}.aligned.bam.bai"
    output:
        flag = out_folder + "{sample}/qc/{sample}.flagstat.txt",
        idx  = out_folder + "{sample}/qc/{sample}.idxstat.txt"
    shell:
        "samtools flagstat {input.bam} > {output.flag}; "
        "samtools idxstats {input.bam} > {output.idx}"


rule rnaseqc:
    input:
        bam = output_bam,
        bai = out_folder + "{sample}/alignment/{sample}.aligned.bam.bai",
        gtf = genes_gtf
    output:
        out_folder + "{sample}/qc/{sample}.metrics.tsv"
    params:
        outdir = out_folder + "{sample}/qc/"
    shell:
        "ml rnaseqc; "
        "rnaseqc {input.gtf} {input.bam} {params.outdir} "
        "--sample={wildcards.sample} "
        "--unpaired --coverage --verbose "
        "--mapping-quality 0 --base-mismatch=1000 --detection-threshold=1"


# FastQC disabled for long reads: its length binning is misleading on ONT/PacBio;
# NanoStat (below) replaces the read-count/length/quality metrics.
# rule fastqc:
#     input:
#         fastq = raw_reads
#     output:
#         html = out_folder + "{sample}/qc/{sample}.raw_fastqc.html",
#         zip  = out_folder + "{sample}/qc/{sample}.raw_fastqc.zip"
#     params:
#         outdir = out_folder + "{sample}/qc/"
#     threads: 8
#     shell:
#         "ml fastqc; "
#         "fastqc --threads {threads} --outdir={params.outdir} {input.fastq}; "
#         "base=$(basename {input.fastq}); "
#         "base=${{base%.gz}}; base=${{base%.fastq}}; base=${{base%.fq}}; "
#         "mv {params.outdir}${{base}}_fastqc.html {output.html}; "
#         "mv {params.outdir}${{base}}_fastqc.zip  {output.zip}"


# NanoStat: long-read count/length/quality (raw + aligned), collated for plotting
rule nanostat_raw:
    input:
        reads = raw_reads
    output:
        out_folder + "{sample}/qc/{sample}.nanostat.raw.tsv"
    threads: 4
    run:
        flag = "--ubam" if prep == "pacbio" else "--fastq"   # pacbio raw = unaligned FLNC bam
        shell("NanoStat %s {input.reads} --tsv -t {threads} > {output}" % flag)


rule nanostat_aligned:
    input:
        bam = output_bam,
        bai = out_folder + "{sample}/alignment/{sample}.aligned.bam.bai"
    output:
        out_folder + "{sample}/qc/{sample}.nanostat.aligned.tsv"
    threads: 4
    shell:
        "NanoStat --bam {input.bam} --tsv -t {threads} > {output}"


rule collate_nanostat:
    input:
        expand(out_folder + "{sample}/qc/{sample}.nanostat.raw.tsv",     sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.nanostat.aligned.tsv", sample=samples)
    output:
        out_folder + data_code + "_nanostat_collated.tsv"
    params:
        script = "scripts/collate_nanostat.py"
    shell:
        "python {params.script} --inFolder {out_folder} -o {output}"


rule picard:
    input:
        bam    = output_bam,
        reflat = reflat_file
    output:
        out_folder + "{sample}/qc/{sample}.RNASeqMetrics"
    shell:
        "ml picard; "
        "java -jar $PICARD CollectRnaSeqMetrics "
        "I={input.bam} O={output} "
        "REF_FLAT={input.reflat} "
        "STRAND=FIRST_READ_TRANSCRIPTION_STRAND "
        "RIBOSOMAL_INTERVALS={rrna} "
        "VALIDATION_STRINGENCY=LENIENT"


rule read_lengths:
    input:
        output_bam
    output:
        mapped   = out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt",
        unmapped = out_folder + "{sample}/qc/{sample}.unmapped.readlengths.txt"
    shell:
        "ml samtools; ml bioawk; "
        "samtools view -F 4 {input} | "
        "bioawk -c sam '{{print length($seq)}}' > {output.mapped}; "
        "samtools view -f 4 {input} | "
        "bioawk -c sam '{{print length($seq)}}' > {output.unmapped}"


# collate_read_lengths.R searches --inFolder recursively (pass out_folder, not the launch dir)
rule collate_lengths:
    input:
        expand(out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt", sample=samples)
    output:
        out_folder + data_code + "_read_lengths_collated.tsv.gz"
    params:
        script = "scripts/collate_read_lengths.R"
    shell:
        "conda activate bambu_env; "
        "Rscript {params.script} --inFolder {out_folder} -o {output}"


# multiqc has no NanoStat module; nanostat is collated separately (above)
rule multiQC:
    input:
        expand(out_folder + "{sample}/qc/{sample}.RNASeqMetrics", sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.metrics.tsv",   sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.flagstat.txt",  sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.idxstat.txt",   sample=samples)
    output:
        out_folder + "multiqc/multiqc_report.html"
    params:
        multiqc = "/sc/arion/projects/als-omics/conda/envs/snakemake/bin/multiqc"
    shell:
        "export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8; "
        "{params.multiqc} -f --outdir {out_folder}multiqc/ {out_folder}"
