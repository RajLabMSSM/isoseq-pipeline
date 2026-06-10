
# qc_pipeline_rawfastq.smk
# ONT / PacBio long-read QC + alignment (minimap2 2.30).
# IDENTICAL to qc_pipeline.smk EXCEPT the `fastqc` rule (and the matching
# multiQC input) run FastQC on the RAW input FASTQ instead of the aligned BAM,
# so the per-base quality / length distributions reflect the original reads
# rather than post-alignment SEQ records (rev-comped, supplementary-clipped).
# Included by Snakefile, which owns shell.prefix, config parsing, samples,
# out_folder, groups, mmi. This file adds only QC-specific constants + rules.
# Brooke Friedman

# ── QC-specific constants ─────────────────────────────────────────────────────
R_VERSION    = "R/4.0.3"
minimap2_bin = "/sc/arion/projects/ad-omics/data/software/minimap2-2.30/minimap2"

genome = referenceFa                       # genome fasta (with .fa), from Snakefile
# rRNA interval list for Picard
rrna   = "/sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/gencode.v30.rRNA.interval.list"

# ── output BAM path (per prep type) ──────────────────────────────────────────
if prep in ("nanopore_cdna", "nanopore_direct"):
    output_bam = out_folder + "{sample}/alignment/{sample}.aligned.bam"
if prep == "pacbio":
    output_bam = out_folder + "{sample}/pbmm2/{sample}.aligned.bam"

# ── minimap2 preset selection (config can override) ───────────────────────────
#   nanopore_cdna   -> splice:hq  (ONT R10.4.1 HAC)
#   nanopore_direct -> splice -k14 (direct RNA)
#   pacbio          -> splice:hq
if prep == "nanopore_cdna":
    minimap_preset = "splice:hq"; minimap_flags = "-uf"
if prep == "nanopore_direct":
    minimap_preset = "splice";    minimap_flags = "-uf -k14"
if prep == "pacbio":
    minimap_preset = "splice:hq"; minimap_flags = "-uf"
if config.get('minimap_preset'):
    minimap_preset = config['minimap_preset']
if config.get('minimap_flags'):
    minimap_flags = config['minimap_flags']

# ── helper: resolve the raw long-read FASTQ for a sample (for FastQC input) ────
def raw_fastq(wildcards):
    return metadata_dict[wildcards.sample]["long_read_fastq"]

# ══════════════════════════════════════════════════════════════════════════════
# REFERENCE PREPARATION  (run once per reference genome)
# ══════════════════════════════════════════════════════════════════════════════
if not config.get('ref_genome_index'):
    rule create_index:
        """Build minimap2 .mmi index from reference FASTA."""
        input:  genome
        output: mmi
        shell:
            "{minimap2_bin} -x {minimap_preset} -d {output} {input}"

# ══════════════════════════════════════════════════════════════════════════════
# ALIGNMENT  (nanopore)
# ══════════════════════════════════════════════════════════════════════════════
rule nanopore_alignment:
    """
    Align ONT cDNA reads to the genome with minimap2 2.30.
    --junc-bed guides splice junction detection; --secondary=no keeps primaries;
    -G 500k allows large introns. SAM is temporary -> sorted BAM in sam_to_bam.
    """
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
    """Convert SAM -> coordinate-sorted BAM -> index."""
    input:
        out_folder + "{sample}/alignment/{sample}.aligned.sam"
    output:
        bam = out_folder + "{sample}/alignment/{sample}.aligned.bam",
        bai = out_folder + "{sample}/alignment/{sample}.aligned.bam.bai"
    threads: 8
    shell:
        "samtools sort -@ {threads} -m 2G -o {output.bam} {input}; "
        "samtools index {output.bam}"


# ══════════════════════════════════════════════════════════════════════════════
# PACBIO ALIGNMENT  (kept for compatibility; only used when prep == pacbio)
# ══════════════════════════════════════════════════════════════════════════════
rule align_flnc_bam:
    """PacBio alignment with pbmm2."""
    input:  mmi
    output: bam = out_folder + "{sample}/pbmm2/{sample}.aligned.bam"
    run:
        input_bam = metadata_dict[wildcards.sample]["flnc_bam_path"]
        shell(
            "pbmm2 align --sort -j 8 --sort-threads 4 -m 3G --preset=ISOSEQ "
            "--log-level INFO --unmapped {mmi} {input_bam} | "
            "samtools calmd -b - {genome} > {output.bam}"
        )


# ══════════════════════════════════════════════════════════════════════════════
# QC METRICS
# ══════════════════════════════════════════════════════════════════════════════
rule samtools_qc:
    """flagstat (mapping summary) + idxstats (per-chrom counts)."""
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
    """RNA-SeQC for long reads (--unpaired; raised mismatch limit)."""
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


rule fastqc:
    """
    FastQC on the RAW long-read FASTQ (pre-alignment input).

    FastQC names its outputs after the input file's basename, which is NOT
    guaranteed to equal {sample} (e.g. CSHL fastqs are combined_PAKxxxxx_*.fastq.gz).
    So we run FastQC, then rename the html/zip to deterministic {sample}.raw_*
    names that the multiQC rule and rule all can predict.
    """
    input:
        fastq = raw_fastq
    output:
        html = out_folder + "{sample}/qc/{sample}.raw_fastqc.html",
        zip  = out_folder + "{sample}/qc/{sample}.raw_fastqc.zip"
    params:
        outdir = out_folder + "{sample}/qc/"
    threads: 8
    shell:
        "ml fastqc; "
        "fastqc --threads {threads} --outdir={params.outdir} {input.fastq}; "
        # FastQC strips a recognized extension off the basename; reproduce that
        # to find what it wrote, then rename to the deterministic {sample} name.
        "base=$(basename {input.fastq}); "
        "base=${{base%.gz}}; base=${{base%.fastq}}; base=${{base%.fq}}; "
        "mv {params.outdir}${{base}}_fastqc.html {output.html}; "
        "mv {params.outdir}${{base}}_fastqc.zip  {output.zip}"


rule picard:
    """Picard CollectRnaSeqMetrics."""
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
    """Extract mapped/unmapped read length distributions from BAM."""
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


rule collate_lengths:
    """
    Collate per-sample read length files into one TSV for plotting.
    NOTE: collate_read_lengths.R searches --inFolder recursively for the per-sample
    *.readlengths.txt files. We pass out_folder so it finds them regardless of the
    directory Snakemake is launched from (the old bug: it searched the launch dir
    and silently wrote an empty file).
    """
    input:
        expand(out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt", sample=samples)
    output:
        out_folder + data_code + "_read_lengths_collated.tsv.gz"
    params:
        script = "scripts/collate_read_lengths.R"
    shell:
        "ml {R_VERSION}; "
        "Rscript {params.script} --inFolder {out_folder} -o {output}"


rule multiQC:
    """Aggregate all QC outputs into a single MultiQC report."""
    input:
        expand(out_folder + "{sample}/qc/{sample}.RNASeqMetrics",     sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.metrics.tsv",       sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.flagstat.txt",      sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.idxstat.txt",       sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.raw_fastqc.zip",    sample=samples)
    output:
        out_folder + "multiqc/multiqc_report.html"
    params:
        multiqc = "/sc/arion/projects/als-omics/conda/envs/snakemake/bin/multiqc"
    shell:
        "export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8; "
        "{params.multiqc} -f --outdir {out_folder}multiqc/ {out_folder}"
