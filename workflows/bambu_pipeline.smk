# bambu_pipeline.smk
# Bambu transcript discovery, run PER GROUP = Join & Call (Jetzinger et al. 2025).
# Bambu (Chen et al. 2023) discovers jointly across all BAMs handed to one bambu()
# call, so passing a group's replicates pools their evidence -> better recall of
# rare, low-expression novel isoforms. Sensitivity is tunable via NDR (higher -> 1.0
# = more novel); see scripts/bambu_run.R.
#
# Included by Snakefile (no shell.prefix / config / rule all here). Inherits
# out_folder, data_code, ref_gtf, referenceFa, groups, group_bams().

BAMBU_R_VERSION = "R/4.2.2"          # Bambu v3 lives on the R/4.2.2 module
bambu_cores     = 8

# Reference annotation prepared once, shared across all groups.
rule bambu_create_annotation:
    input:
        gtf = ref_gtf
    output:
        rdata = out_folder + "bambu/bambu_annotation.RData"
    params:
        script = "scripts/bambu_annotation.R"
    shell:
        "conda activate snakemake; ml {BAMBU_R_VERSION}; "
        "Rscript {params.script} -i {input.gtf} -o {output.rdata}"


rule bambu_run:
    """
    Join & Call within a group: all of the group's aligned BAMs go to ONE bambu()
    call so discovery pools evidence across replicates. To raise sensitivity for
    rare novel isoforms, add `--NDR <0-1>` (higher = more novel) — see bambu_run.R.
    """
    input:
        bams = lambda wc: group_bams(wc.group),
        anno = out_folder + "bambu/bambu_annotation.RData"
    output:
        gtf    = out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf",
        counts = out_folder + "bambu/{group}/" + data_code + "_{group}_counts_transcript.txt",
        rdata  = out_folder + "bambu/{group}/" + data_code + "_{group}_bambu.RData"
    params:
        prefix = out_folder + "bambu/{group}/" + data_code + "_{group}",
        script = "scripts/bambu_run.R"
    threads: bambu_cores
    shell:
        "conda activate snakemake; ml {BAMBU_R_VERSION}; "
        "Rscript {params.script} --cores {threads} --fasta {referenceFa} "
        "--anno {input.anno} --prefix {params.prefix} {input.bams}"
        # rare-isoform tuning (optional): append  --NDR 0.5  (or higher, up to 1.0)
