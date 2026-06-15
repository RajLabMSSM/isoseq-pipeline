# consensus_pipeline.smk
# Cross-tool consensus per group: gffcompare matches the three assembler GTFs by
# intron chain; keep transcripts present in >= min_consensus_tools of them. First
# of two precision filters (consensus here, SQANTI next). Query order = tracking
# column order: q1 bambu, q2 isoquant, q3 stringtie3.
# merge_consensus then unions the per-group consensus GTFs into ONE per-cohort GTF
# (dedup by intron chain) — that single GTF is what goes to SQANTI.

gffcompare = "/sc/arion/projects/ad-omics/data/software/gffcompare-0.12.10.Linux_x86_64/gffcompare"
min_consensus_tools = config.get("min_consensus_tools", 2)


# 5'/3' ends ignored (good for ONT end-fuzziness); --multi-exon-only excludes single-exon
rule consensus_isoforms:
    input:
        bambu     = out_folder + "bambu/{group}/" + data_code + "_{group}_extended_annotations.gtf",
        isoquant  = out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf",
        stringtie = out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf",
        ref       = ref_gtf
    output:
        combined = out_folder + "consensus/{group}/" + data_code + "_{group}_gffcmp.combined.gtf",
        tracking = out_folder + "consensus/{group}/" + data_code + "_{group}_gffcmp.tracking",
        gtf      = out_folder + "consensus/{group}/" + data_code + "_{group}_consensus.gtf"
    params:
        prefix    = out_folder + "consensus/{group}/" + data_code + "_{group}_gffcmp",
        min_tools = min_consensus_tools,
        script    = "scripts/consensus_from_gffcompare.py",
        python    = "/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/python"
    shell:
        "conda activate isoseq-pipeline; "
        # gffcompare writes {prefix}.combined.gtf, {prefix}.tracking, .loci, .stats
        "{gffcompare} -r {input.ref} -o {params.prefix} "
        "{input.bambu} {input.isoquant} {input.stringtie}; "
        "{params.python} {params.script} "
        "--tracking {output.tracking} --combined {output.combined} "
        "--min-tools {params.min_tools} -o {output.gtf}"


# union the per-group consensus GTFs into one per-cohort reference (gffcompare
# .combined collapses identical intron chains); this single GTF feeds SQANTI
rule merge_consensus:
    input:
        gtfs = expand(out_folder + "consensus/{group}/" + data_code + "_{group}_consensus.gtf", group=groups),
        ref  = ref_gtf
    output:
        gtf = out_folder + "consensus/" + data_code + "_merged_consensus.gtf"
    params:
        prefix = out_folder + "consensus/" + data_code + "_merged_gffcmp"
    shell:
        "conda activate isoseq-pipeline; "
        "{gffcompare} -r {input.ref} -o {params.prefix} {input.gtfs}; "
        "cp {params.prefix}.combined.gtf {output.gtf}"
