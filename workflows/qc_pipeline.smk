# tool path and the per-sample / per-cohort files shared by the qc rules
minimap2_bin      = "/sc/arion/projects/ad-omics/data/software/minimap2-2.31/minimap2"
output_cram       = out_folder + "{sample}/alignment/{sample}.aligned.cram"
qc_bam            = out_folder + "{sample}/alignment/{sample}.qc.bam"
rrna              = out_folder + data_code + ".rRNA.interval_list"
combined_junc_bed = out_folder + data_code + ".combined_junctions.bed"
pychopper_stats   = expand(out_folder + "{sample}/pychopper/{sample}.stats.tsv", sample=samples)
# the panels of the multiqc page, written by longread_qc_tables.py into the folder multiqc scans
longread_mqc_panels = expand(out_folder + data_code + "_{panel}_mqc.tsv",
                             panel=["longread_yields", "pychopper_classes", "pychopper_reads", "pychopper_lengths",
                                    "minimap_reads", "minimap_lengths"])


# raw long-read fastq of one sample, from the sample sheet
def raw_fastq(wildcards):
    return metadata_dict[wildcards.sample]["long_read_fastq"]


# minimap2 --junc-bed: the GENCODE junction bed plus the cohort's STAR introns (>= short_read_junc_min_uniq
# unique reads, canonical motif), so ONT reads are aligned to splice sites the short reads confirm.
# exclude_samples are dropped here (tdpkd: TDP-5 aberrant differentiation, CTRL-4 low-depth outlier).
rule build_junc_bed:
    input:
        anno = config['junctions_bed'],
        sj   = junction_folder,
        # the script is an input, so editing it re-runs this step
        script = "scripts/build_junc_bed.py"
    output:
        combined_junc_bed
    params:
        exclude = ",".join(sorted(excluded_samples))
    shell:
        """
        python {input.script} \
        --annotation {input.anno} \
        --sj_folder {input.sj} \
        --min_uniq {config[short_read_junc_min_uniq]} --canonical_only \
        --exclude_samples '{params.exclude}' \
        -o {output}
        """


# pychopper finds the cDNA primers of the kit (-k) at both read ends (-m edlib alignment), orients
# each read to the transcript strand and trims the primers. kept = reads with both primers (primary)
# plus reads rescued from fused molecules (-w); unclassified reads where there is no valid primer arrangement (-u) are deleted. 
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
        mkdir -p {params.pdir}                                            # make output directory
        zcat -f {input.fastq} > {params.pre}.input.fastq                  # decompress input FASTQ
        pychopper -k {config[pychopper_kit]} -m edlib -t {threads} \
        -r {params.pre}.report.pdf \
        -u {params.pre}.unclassified.fastq \
        -w {params.pre}.rescued.fastq \
        -S {params.pre}.stats.tsv \
        {params.pre}.input.fastq {params.pre}.primary.fastq
        cat {params.pre}.primary.fastq {params.pre}.rescued.fastq > {output.fl}  # combines full-length and rescued FASTQs into one
        rm -f {params.pre}.input.fastq {params.pre}.primary.fastq {params.pre}.unclassified.fastq {params.pre}.rescued.fastq
        """


# spliced alignment with minimap2 (Li 2018): the preset comes from the config. splice already looks for
# GT-AG splice signals on both strands (-ub). k is fixed when the .mmi index is built, it's k=15 which is default. --secondary=no writes no secondary alignments.
# -G 500k widens the intron search to 500 kb, above the default 200 kb; it is not a hard limit, so some
# introns come out longer. --junc-bed gives a bonus at the splice sites in the combined bed.
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
        -ax {config[minimap_preset]} \
        --secondary=no -G 500k \
        --junc-bed {input.junc} \
        -t {threads} \
        {input.index} {input.fastq} \
        > {output.sam}
        """


# sort and store the alignments as cram against the genome fasta
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


# read and per-chromosome mapping counts
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


# rnaseqc gene / exon metrics. --unpaired counts single reads, which rnaseqc otherwise skips.
# --mapping-quality 0 and --base-mismatch=1000 keep every aligned long read (nanopore reads carry many
# mismatches, the short-read defaults would discard most of them); --detection-threshold=1 calls a
# gene detected from 1 read; --coverage also writes per-transcript coverage
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


# read count and length statistics of the raw reads (before pychopper) - summary stats only
rule nanostat_raw:
    input:
        reads = raw_fastq
    output:
        out_folder + "{sample}/qc/{sample}.nanostat.raw.tsv"
    threads: 4
    shell:
        "NanoStat --fastq {input.reads} --tsv -t {threads} > {output}"


# the same statistics over the reads pychopper kept (the fastq that goes to minimap2), the middle stage of the
# read counts: raw -> after pychopper -> after minimap2
rule nanostat_pychopper:
    input:
        reads = out_folder + "{sample}/pychopper/{sample}.full_length.fastq"
    output:
        out_folder + "{sample}/qc/{sample}.nanostat.pychopper.tsv"
    threads: 4
    shell:
        "NanoStat --fastq {input.reads} --tsv -t {threads} > {output}"


# the same statistics over the primary mapped reads (--no_supplementary, so one record per read) - summary stats only
rule nanostat_aligned:
    input:
        bam = qc_bam,
        bai = qc_bam + ".bai"
    output:
        out_folder + "{sample}/qc/{sample}.nanostat.aligned.tsv"
    threads: 4
    shell:
        "NanoStat --bam {input.bam} --no_supplementary --tsv -t {threads} > {output}"


# rRNA interval list for picard: the awk keeps the gene lines whose gene_type is rRNA,
# rRNA_pseudogene or Mt_rRNA and prints them as bed (0-based start, gene_id as the name, strand);
# picard turns the bed into an interval list against the genome's sequence dictionary
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


# picard rna-seq metrics (coding / utr / intronic / rRNA shares). pychopper has oriented every read
# to the transcript strand, so the read itself is on the transcription strand
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


# length of every primary mapped read, for the read-length histogram (-F 2308 drops unmapped, secondary and
# supplementary records, so each read is counted once at its full trimmed length)
rule read_lengths:
    input:
        qc_bam
    output:
        out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt"
    shell:
        """
        ml samtools
        ml bioawk
        samtools view -F 2308 {input} | bioawk -c sam '{{print length($seq)}}' > {output}
        """


# per-cohort qc tables: read-length histogram, the per-sample summary and every NanoStat metric at the three stages
rule longread_qc_tables:
    input:
        lengths = expand(out_folder + "{sample}/qc/{sample}.mapped.readlengths.txt", sample=samples),
        raw     = expand(out_folder + "{sample}/qc/{sample}.nanostat.raw.tsv",       sample=samples),
        chopped = expand(out_folder + "{sample}/qc/{sample}.nanostat.pychopper.tsv", sample=samples),
        aligned = expand(out_folder + "{sample}/qc/{sample}.nanostat.aligned.tsv",   sample=samples),
        pyc     = pychopper_stats,
        # the sample sheet and the script are inputs, so editing either redoes the qc tables
        metadata = config['metadata'],
        script  = "scripts/longread_qc_tables.py"
    output:
        hist    = out_folder + data_code + "_read_length_hist.tsv.gz",
        summary = out_folder + data_code + "_longread_qc_summary.tsv",
        all_nanostat = out_folder + data_code + "_nanostat_all_metrics.tsv",
        panels  = longread_mqc_panels
    params:
        mqc_prefix = out_folder + data_code
    shell:
        """
        python {input.script} \
        --metadata {input.metadata} \
        --read-lengths {input.lengths} \
        --nanostat-raw {input.raw} \
        --nanostat-pychopper {input.chopped} \
        --nanostat-aligned {input.aligned} \
        --pychopper {input.pyc} \
        --hist-out {output.hist} \
        --summary-out {output.summary} \
        --nanostat-out {output.all_nanostat} \
        --mqc-prefix {params.mqc_prefix}
        """


# nanostat is excluded: multiqc names every NanoStat --tsv report "dataset". the *_mqc.tsv panels add the read flow plots
# and the yields; multiqc_config.yaml hides the parts of multiqc's own output that those replace
rule multiQC:
    input:
        expand(out_folder + "{sample}/qc/{sample}.RNASeqMetrics", sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.metrics.tsv",   sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.flagstat.txt",  sample=samples),
        expand(out_folder + "{sample}/qc/{sample}.idxstat.txt",   sample=samples),
        pychopper_stats,
        longread_mqc_panels,
        settings = "scripts/multiqc_config.yaml"
    output:
        out_folder + "multiqc/multiqc_report.html"
    params:
        multiqc = "/sc/arion/projects/als-omics/conda/envs/multiqc/bin/multiqc"
    shell:
        """
        export LC_ALL=en_US.UTF-8
        export LANG=en_US.UTF-8
        {params.multiqc} -f -e nanostat -c {input.settings} --outdir {out_folder}multiqc/ {out_folder}
        """
