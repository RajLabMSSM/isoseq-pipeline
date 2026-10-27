minimap2_bin      = "/sc/arion/projects/ad-omics/data/software/minimap2-2.31/minimap2"
output_cram       = out_folder + "{sample}/alignment/{sample}.aligned.cram"
qc_bam            = out_folder + "{sample}/alignment/{sample}.qc.bam"
rrna              = out_folder + data_code + ".rRNA.interval_list"
combined_junc_bed = out_folder + data_code + ".combined_junctions.bed"
pychopper_stats   = expand(out_folder + "{sample}/pychopper/{sample}.stats.tsv", sample=samples)


def raw_fastq(wildcards):
    return metadata_dict[wildcards.sample]["long_read_fastq"]


# minimap2 --junc-bed: the GENCODE junction bed plus the cohort's STAR introns (>= short_read_junc_min_uniq
# unique reads, canonical motif), so ONT reads are aligned to splice sites the short reads confirm.
# exclude_samples are dropped here (tdpkd: TDP-5 aberrant differentiation, CTRL-4 low-depth outlier).
# The published v9 build (2026-09-21) did not drop them here: its tdpkd junc-bed had 508,354 short-read
# junctions, vs 500,324 with the two dropped.
rule build_junc_bed:
    input:
        anno = config['junctions_bed'],
        sj   = junction_folder
    output:
        combined_junc_bed
    params:
        exclude = ",".join(sorted(excluded_samples))
    shell:
        """
        python scripts/build_junc_bed.py \
        --annotation {input.anno} \
        --sj_folder {input.sj} \
        --min_uniq {config[short_read_junc_min_uniq]} --canonical_only \
        --exclude_samples '{params.exclude}' \
        -o {output}
        """


rule pychopper:
    input:
        fastq = raw_fastq
    output:
        fl    = temp(out_folder + "{sample}/pychopper/{sample}.full_length.fastq"),
        stats = out_folder + "{sample}/pychopper/{sample}.stats.tsv"
    params:
        pdir = out_folder + "{sample}/pychopper",
        pre  = out_folder + "{sample}/pychopper/{sample}"
    threads: 8
    shell:
        """
        conda activate pychopper_env
        mkdir -p {params.pdir}
        zcat -f {input.fastq} > {params.pre}.input.fastq
        pychopper -k {config[pychopper_kit]} -m edlib -t {threads} \
        -r {params.pre}.report.pdf \
        -u {params.pre}.unclassified.fastq \
        -w {params.pre}.rescued.fastq \
        -S {params.pre}.stats.tsv \
        {params.pre}.input.fastq {params.pre}.primary.fastq
        cat {params.pre}.primary.fastq {params.pre}.rescued.fastq > {output.fl}
        rm -f {params.pre}.input.fastq {params.pre}.primary.fastq {params.pre}.unclassified.fastq {params.pre}.rescued.fastq
        """


rule nanopore_alignment:
    input:
        index = config['ref_genome_index'],
        junc  = combined_junc_bed,
        fastq = out_folder + "{sample}/pychopper/{sample}.full_length.fastq"
    output:
        sam = temp(out_folder + "{sample}/alignment/{sample}.aligned.sam")
    threads: 8
    shell:
        """
        {minimap2_bin} \
        -ax {config[minimap_preset]} {config[minimap_flags]} \
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


# temporary bam for the QC tools that cannot read cram
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


rule nanostat_raw:
    input:
        reads = raw_fastq
    output:
        out_folder + "{sample}/qc/{sample}.nanostat.raw.tsv"
    threads: 4
    shell:
        "NanoStat --fastq {input.reads} --tsv -t {threads} > {output}"


rule nanostat_aligned:
    input:
        bam = qc_bam,
        bai = qc_bam + ".bai"
    output:
        out_folder + "{sample}/qc/{sample}.nanostat.aligned.tsv"
    threads: 4
    shell:
        "NanoStat --bam {input.bam} --tsv -t {threads} > {output}"


rule rrna_intervals:
    input:
        gtf    = ref_gtf,
        genome = genome
    output:
        rrna
    params:
        bed  = out_folder + data_code + ".rRNA.bed",
        dict = out_folder + data_code + ".genome.dict"
    shell:
        """
        ml picard
        awk -F'\\t' -v OFS='\\t' '$3 == "gene" && $9 ~ /gene_type "(rRNA|rRNA_pseudogene|Mt_rRNA)"/ {{ match($9, /gene_id "[^"]+"/); print $1, $4 - 1, $5, substr($9, RSTART + 9, RLENGTH - 10), 0, $7 }}' {input.gtf} > {params.bed}
        java -jar $PICARD CreateSequenceDictionary R={input.genome} O={params.dict}
        java -jar $PICARD BedToIntervalList I={params.bed} O={output} SD={params.dict}
        rm -f {params.bed} {params.dict}
        """


rule picard:
    input:
        bam    = qc_bam,
        reflat = config['reflat_file'],
        rrna   = rrna
    output:
        out_folder + "{sample}/qc/{sample}.RNASeqMetrics"
    shell:
        """
        ml picard
        java -jar $PICARD CollectRnaSeqMetrics \
        I={input.bam} O={output} \
        REF_FLAT={input.reflat} \
        STRAND=FIRST_READ_TRANSCRIPTION_STRAND \
        RIBOSOMAL_INTERVALS={input.rrna} \
        VALIDATION_STRINGENCY=LENIENT
        """


rule read_lengths:
    input:
        qc_bam
    output:
        out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt"
    shell:
        """
        ml samtools
        ml bioawk
        samtools view -F 4 {input} | bioawk -c sam '{{print length($seq)}}' > {output}
        """


rule longread_qc_tables:
    input:
        lengths = expand(out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt", sample=samples),
        raw     = expand(out_folder + "{sample}/qc/{sample}.nanostat.raw.tsv",       sample=samples),
        aligned = expand(out_folder + "{sample}/qc/{sample}.nanostat.aligned.tsv",   sample=samples),
        pyc     = pychopper_stats
    output:
        hist    = out_folder + data_code + "_read_length_hist.tsv.gz",
        summary = out_folder + data_code + "_longread_qc_summary.tsv"
    shell:
        """
        python scripts/longread_qc_tables.py \
        --metadata {config[metadata]} \
        --read-lengths {input.lengths} \
        --nanostat-raw {input.raw} \
        --nanostat-aligned {input.aligned} \
        --pychopper {input.pyc} \
        --hist-out {output.hist} \
        --summary-out {output.summary}
        """


# nanostat is excluded: multiqc names every NanoStat --tsv report "dataset"
rule multiQC:
    input:
        expand(out_folder + "{sample}/qc/{sample}.RNASeqMetrics", sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.metrics.tsv",   sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.flagstat.txt",  sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.idxstat.txt",   sample=samples),
        pychopper_stats
    output:
        out_folder + "multiqc/multiqc_report.html"
    params:
        multiqc = "/sc/arion/projects/als-omics/conda/envs/multiqc/bin/multiqc"
    shell:
        """
        export LC_ALL=en_US.UTF-8
        export LANG=en_US.UTF-8
        {params.multiqc} -f -e nanostat --outdir {out_folder}multiqc/ {out_folder}
        """
