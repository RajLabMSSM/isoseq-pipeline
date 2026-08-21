# StringTie3 (Shinder et al. 2025) assembly, run once per group on the pooled bam

stringtie         = "/sc/arion/projects/ad-omics/data/software/stringtie-3.0.3.Linux_x86_64/stringtie"
stringtie_threads = 4
stringtie_mode    = config.get("stringtie_mode", "long_read")   # long_read uses -L, hybrid uses --mix

# -t turns off coverage-based trimming of transcript ends. only has an effect in --mix mode
stringtie_mix_end_trim_flag = "-t " if config.get("stringtie_mix_disable_end_trim", True) else ""


# short-read bams for a group, skipping samples with no path in the metadata
def group_short_bams(group):
    paths = []
    for s in sample_groups[group]:
        p = metadata_dict[s].get("short_read_bam_path")
        if p is not None and str(p) not in ("nan", "NA", ""):
            paths.append(p)
    return paths


# long reads are pooled by pool_group_bams in the Snakefile, this is the short-read equivalent
rule stringtie_pool_short_bams:
    input:
        bams = lambda wc: group_short_bams(wc.group)
    output:
        bam = temp(out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.short.bam"),
        bai = temp(out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.short.bam.bai")
    params:
        ref = ref_genome + ".fa"
    threads: 4
    shell:
        """
        conda activate isoseq-pipeline
        samtools merge -f -@ {threads} --reference {params.ref} {output.bam} {input.bams}
        samtools index {output.bam}
        """


def stringtie_run_inputs(wildcards):
    d = {
        "gtf":      ref_gtf,
        "long_bam": out_folder + "pooled/%s_%s.pooled.bam" % (data_code, wildcards.group),
        "long_bai": out_folder + "pooled/%s_%s.pooled.bam.bai" % (data_code, wildcards.group),
    }
    if stringtie_mode == "hybrid":
        d["short_bam"] = out_folder + "stringtie3/%s/%s_%s.pooled.short.bam" % (wildcards.group, data_code, wildcards.group)
    return d


# --mix expects the short-read bam first and the long-read bam second
rule run_stringtie:
    input:
        unpack(stringtie_run_inputs)
    output:
        gtf = out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf"
    threads: stringtie_threads
    run:
        if stringtie_mode == "hybrid":
            shell("{stringtie} -p {threads} --mix {stringtie_mix_end_trim_flag}-G {input.gtf} -o {output.gtf} "
                  "{input.short_bam} {input.long_bam}")
        else:
            shell("{stringtie} -p {threads} -L -G {input.gtf} -o {output.gtf} {input.long_bam}")
