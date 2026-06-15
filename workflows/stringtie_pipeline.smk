# stringtie_pipeline.smk
# StringTie3 (Shinder et al. 2025) assembly, per group (Join & Call): assemble the
# one pooled BAM. Mode via config stringtie_mode:
#   long_read -> -L ; hybrid -> --mix (needs short_read_bam_path in metadata)

stringtie         = "/sc/arion/projects/ad-omics/data/software/stringtie-3.0.3.Linux_x86_64/stringtie"
stringtie_threads = 4
stringtie_mode    = config.get("stringtie_mode", "long_read")   # "long_read" | "hybrid"


# a group's short-read BAMs (hybrid mode), skipping NA/missing
def group_short_bams(group):
    paths = []
    for s in sample_groups[group]:
        p = metadata_dict[s].get("short_read_bam_path")
        if p is not None and str(p) not in ("nan", "NA", ""):
            paths.append(p)
    return paths


# long-read pooling is the shared pool_group_bams rule (Snakefile); this is short-read only
rule stringtie_pool_short_bams:
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
        "long_bam": out_folder + "pooled/%s_%s.pooled.bam" % (data_code, wildcards.group),
    }
    if stringtie_mode == "hybrid":
        d["short_bam"] = out_folder + "stringtie3/%s/%s_%s.pooled.short.bam" % (wildcards.group, data_code, wildcards.group)
    return d


# -L = long-read mode; --mix takes the long bam as the 2nd input
rule run_stringtie:
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
