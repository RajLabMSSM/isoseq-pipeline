# Alignment and QC. one aligner is active per run, selected by prep

minimap2_bin = "/sc/arion/projects/ad-omics/data/software/minimap2-2.31/minimap2"
genome = referenceFa
rrna   = "/sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/gencode.v30.rRNA.interval.list"

output_cram = out_folder + "{sample}/alignment/{sample}.aligned.cram"
qc_bam      = out_folder + "{sample}/alignment/{sample}.qc.bam"

# minimap2 preset per library prep, overridable from config
if prep == "nanopore_cdna":
    minimap_preset = "splice:hq"; minimap_flags = "-ub"
if prep == "nanopore_direct":
    minimap_preset = "splice";    minimap_flags = "-uf -k14"
if prep == "pacbio":
    minimap_preset = "splice:hq"; minimap_flags = "-uf"
if config.get('minimap_preset'):
    minimap_preset = config['minimap_preset']
if config.get('minimap_flags'):
    minimap_flags = config['minimap_flags']

combined_junc_bed = out_folder + data_code + ".combined_junctions.bed"


# raw pre-alignment reads for NanoStat, ont fastq or pacbio unaligned flnc bam
def raw_reads(wildcards):
    if prep == "pacbio":
        return pacbio_flnc(wildcards.sample)
    return metadata_dict[wildcards.sample]["long_read_fastq"]


def junc_bed_for_alignment():
    return combined_junc_bed if use_short_read_junctions else junctions_bed


# pychopper output for cdna, otherwise the raw long reads
def align_fastq(wildcards):
    if prep == "nanopore_cdna" and use_pychopper:
        return out_folder + "{s}/pychopper/{s}.full_length.fastq".format(s=wildcards.sample)
    return metadata_dict[wildcards.sample]["long_read_fastq"]


#####
## Alignment
#####

# annotation bed12 plus empirical STAR introns, combined into one --junc-bed
if use_short_read_junctions:

    rule build_junc_bed:
        input:
            anno = junctions_bed,
            sj   = junction_folder
        output:
            combined_junc_bed
        params:
            script    = "scripts/build_junc_bed.py",
            min_uniq  = short_read_junc_min_uniq,
            canonical = "--canonical_only" if short_read_junc_canonical_only else ""
        shell:
            """
            python {params.script} \
            --annotation {input.anno} \
            --sj_folder {input.sj} \
            --min_uniq {params.min_uniq} {params.canonical} \
            -o {output}
            """


if not config.get('ref_genome_index'):

    rule create_index:
        input:
            genome
        output:
            mmi
        shell:
            "{minimap2_bin} -x {minimap_preset} -d {output} {input}"


# orient and trim full-length cdna reads, and rescue fused reads
if prep == "nanopore_cdna" and use_pychopper:

    rule pychopper:
        input:
            fastq = lambda wc: metadata_dict[wc.sample]["long_read_fastq"]
        output:
            fl = temp(out_folder + "{sample}/pychopper/{sample}.full_length.fastq")
        params:
            pdir  = out_folder + "{sample}/pychopper",
            pre   = out_folder + "{sample}/pychopper/{sample}",
            kit   = pychopper_kit,
            qflag = ("-q " + str(pychopper_q)) if pychopper_q is not None else ""
        threads: 8
        shell:
            """
            conda activate pychopper_env
            mkdir -p {params.pdir}
            zcat -f {input.fastq} > {params.pre}.input.fastq
            pychopper -k {params.kit} -m edlib -t {threads} {params.qflag} \
            -r {params.pre}.report.pdf \
            -u {params.pre}.unclassified.fastq \
            -w {params.pre}.rescued.fastq \
            -S {params.pre}.stats.tsv \
            {params.pre}.input.fastq {params.pre}.primary.fastq
            cat {params.pre}.primary.fastq {params.pre}.rescued.fastq > {output.fl}
            rm -f {params.pre}.input.fastq {params.pre}.primary.fastq {params.pre}.unclassified.fastq {params.pre}.rescued.fastq
            """


if prep in ("nanopore_cdna", "nanopore_direct"):

    rule nanopore_alignment:
        input:
            index = mmi,
            junc  = junc_bed_for_alignment(),
            fastq = align_fastq
        output:
            sam = temp(out_folder + "{sample}/alignment/{sample}.aligned.sam")
        threads: 8
        shell:
            """
            {minimap2_bin} \
            -ax {minimap_preset} {minimap_flags} \
            --secondary=no -G 500k \
            --junc-bed {input.junc} \
            -t {threads} \
            {input.index} {input.fastq} \
            > {output.sam}
            """

    rule sam_to_cram:
        input:
            out_folder + "{sample}/alignment/{sample}.aligned.sam"
        output:
            cram = output_cram,
            crai = output_cram + ".crai"
        threads: 8
        shell:
            """
            samtools sort -@ {threads} -m 2G -O bam {input} | \
            samtools view -@ {threads} -C -T {genome} -o {output.cram} -
            samtools index {output.cram}
            """


if prep == "pacbio":

    # the minimap2 .mmi is reused as the pbmm2 index
    rule pacbio_alignment:
        input:
            index = mmi,
            flnc  = lambda wc: pacbio_flnc(wc.sample)
        output:
            cram = output_cram,
            crai = output_cram + ".crai"
        threads: 8
        shell:
            """
            pbmm2 align --sort -j {threads} --sort-threads 4 -m 3G \
            --preset=ISOSEQ --log-level INFO --unmapped \
            {input.index} {input.flnc} | \
            samtools calmd -b - {genome} | \
            samtools view -C -T {genome} -o {output.cram} -
            samtools index {output.cram}
            """


# alignments are stored as cram, so QC tools that cannot read cram get a temporary bam
rule cram_to_qc_bam:
    input:
        cram = output_cram,
        crai = output_cram + ".crai"
    output:
        bam = temp(qc_bam),
        bai = temp(qc_bam + ".bai")
    threads: 4
    shell:
        """
        samtools view -@ {threads} -b -T {genome} -o {output.bam} {input.cram}
        samtools index {output.bam}
        """

#####
## QC metrics
#####

rule samtools_qc:
    input:
        bam = qc_bam,
        bai = qc_bam + ".bai"
    output:
        flag = out_folder + "{sample}/qc/{sample}.flagstat.txt",
        idx  = out_folder + "{sample}/qc/{sample}.idxstat.txt"
    shell:
        """
        samtools flagstat {input.bam} > {output.flag}
        samtools idxstats {input.bam} > {output.idx}
        """


rule rnaseqc:
    input:
        bam = qc_bam,
        bai = qc_bam + ".bai",
        gtf = genes_gtf
    output:
        out_folder + "{sample}/qc/{sample}.metrics.tsv"
    params:
        outdir = out_folder + "{sample}/qc/"
    shell:
        """
        ml rnaseqc
        rnaseqc {input.gtf} {input.bam} {params.outdir} \
        --sample={wildcards.sample} \
        --unpaired --coverage --verbose \
        --mapping-quality 0 --base-mismatch=1000 --detection-threshold=1
        """


# fastqc is deliberately not run, its length binning is misleading on long reads.
# NanoStat covers the same read-level metrics
rule nanostat_raw:
    input:
        reads = raw_reads
    output:
        out_folder + "{sample}/qc/{sample}.nanostat.raw.tsv"
    threads: 4
    run:
        flag = "--ubam" if prep == "pacbio" else "--fastq"
        shell("NanoStat %s {input.reads} --tsv -t {threads} > {output}" % flag)


rule nanostat_aligned:
    input:
        bam = qc_bam,
        bai = qc_bam + ".bai"
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
        bam    = qc_bam,
        reflat = reflat_file
    output:
        out_folder + "{sample}/qc/{sample}.RNASeqMetrics"
    shell:
        """
        ml picard
        java -jar $PICARD CollectRnaSeqMetrics \
        I={input.bam} O={output} \
        REF_FLAT={input.reflat} \
        STRAND=FIRST_READ_TRANSCRIPTION_STRAND \
        RIBOSOMAL_INTERVALS={rrna} \
        VALIDATION_STRINGENCY=LENIENT
        """


rule read_lengths:
    input:
        qc_bam
    output:
        mapped   = out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt",
        unmapped = out_folder + "{sample}/qc/{sample}.unmapped.readlengths.txt"
    shell:
        """
        ml samtools
        ml bioawk
        samtools view -F 4 {input} | bioawk -c sam '{{print length($seq)}}' > {output.mapped}
        samtools view -f 4 {input} | bioawk -c sam '{{print length($seq)}}' > {output.unmapped}
        """


# the script searches --inFolder recursively, so it takes out_folder rather than the launch dir
rule collate_lengths:
    input:
        expand(out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt", sample=samples)
    output:
        out_folder + data_code + "_read_lengths_collated.tsv.gz"
    params:
        script = "scripts/collate_read_lengths.R"
    shell:
        """
        conda activate bambu_env
        Rscript {params.script} --inFolder {out_folder} -o {output}
        """


# multiqc has no NanoStat module, so the collated table is written as custom content
# (any *_mqc.tsv is picked up generically). this adds read-level stats to the
# alignment-level metrics the modules already cover
rule nanostat_mqc:
    input:
        out_folder + data_code + "_nanostat_collated.tsv"
    output:
        raw     = out_folder + "nanostat_raw_mqc.tsv",
        aligned = out_folder + "nanostat_aligned_mqc.tsv"
    params:
        script = "scripts/nanostat_mqc.py",
        prefix = out_folder + "nanostat"
    shell:
        "python {params.script} -i {input} -o {params.prefix}"


rule multiQC:
    input:
        expand(out_folder + "{sample}/qc/{sample}.RNASeqMetrics", sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.metrics.tsv",   sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.flagstat.txt",  sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.idxstat.txt",   sample=samples),
        nanostat_raw     = out_folder + "nanostat_raw_mqc.tsv",
        nanostat_aligned = out_folder + "nanostat_aligned_mqc.tsv"
    output:
        out_folder + "multiqc/multiqc_report.html"
    params:
        multiqc = "/sc/arion/projects/als-omics/conda/envs/snakemake/bin/multiqc"
    shell:
        """
        export LC_ALL=en_US.UTF-8
        export LANG=en_US.UTF-8
        {params.multiqc} -f --outdir {out_folder}multiqc/ {out_folder}
        """
