# stringtie_pipeline.smk
# StringTie3 (Shinder et al. 2025) transcript assembly, run PER GROUP = Join & Call
# (Jetzinger et al. 2025).
#
# StringTie's native multi-sample flow (-L per sample + `stringtie --merge`) is
# Call & Join. To do Join & Call we merge each group's aligned BAMs into ONE pooled
# BAM and assemble it once, so rare isoforms get pooled read support and aren't lost
# to per-sample coverage thresholds.
#
# Mode is config-driven:
#   stringtie_mode: "long_read"  -> stringtie3 -L         (ONT/PacBio long reads)
#   stringtie_mode: "hybrid"     -> stringtie3 --mix      (pools short BAMs too;
#                                   requires short_read_bam_path in metadata)
#
# Included by Snakefile (no shell.prefix / config / rule all here).

stringtie         = "/sc/arion/projects/ad-omics/data/software/stringtie-3.0.0.Linux_x86_64/stringtie"
stringtie_threads = 4
stringtie_mode    = config.get("stringtie_mode", "long_read")   # "long_read" | "hybrid"


def group_short_bams(group):
    """Short-read BAMs for a group (hybrid mode), skipping NA/missing entries."""
    paths = []
    for s in sample_groups[group]:
        p = metadata_dict[s].get("short_read_bam_path")
        if p is not None and str(p) not in ("nan", "NA", ""):
            paths.append(p)
    return paths


rule stringtie_pool_long_bams:
    """Pool a group's long-read aligned BAMs into one BAM (the Join & Call unit)."""
    input:
        bams = lambda wc: group_bams(wc.group)
    output:
        bam = out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.bam",
        bai = out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.bam.bai"
    threads: 4
    shell:
        "conda activate isoseq-pipeline; "
        "samtools merge -f -@ {threads} {output.bam} {input.bams}; "
        "samtools index {output.bam}"


rule stringtie_pool_short_bams:
    """Pool a group's short-read BAMs (hybrid mode only)."""
    input:
        bams = lambda wc: group_short_bams(wc.group)
    output:
        bam = out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.short.bam",
        bai = out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.short.bam.bai"
    threads: 4
    shell:
        "conda activate isoseq-pipeline; "
        "samtools merge -f -@ {threads} {output.bam} {input.bams}; "
        "samtools index {output.bam}"


def stringtie_run_inputs(wildcards):
    d = {
        "gtf":      ref_gtf,
        "long_bam": out_folder + "stringtie3/%s/%s_%s.pooled.bam" % (wildcards.group, data_code, wildcards.group),
    }
    if stringtie_mode == "hybrid":
        d["short_bam"] = out_folder + "stringtie3/%s/%s_%s.pooled.short.bam" % (wildcards.group, data_code, wildcards.group)
    return d


rule stringtie_run:
    """
    Single assembly on the pooled BAM (Join & Call). -L enforces long-read mode
    (and -s 1.5 -g 0); --mix expects the long-read BAM as the 2nd input.
    """
    input:
        unpack(stringtie_run_inputs)
    output:
        gtf = out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf"
    threads: stringtie_threads
    run:
        if stringtie_mode == "hybrid":
            shell("{stringtie} -p {threads} --mix -G {input.gtf} -o {output.gtf} "
                  "{input.short_bam} {input.long_bam}")
        else:
            shell("{stringtie} -p {threads} -L -G {input.gtf} -o {output.gtf} {input.long_bam}")
