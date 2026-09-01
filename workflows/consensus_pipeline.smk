# End-aware consensus across the three assemblers, then the final reference + novel gtf
#
# build_union.py concatenates the assembler gtfs (3 tools x n groups) into one union with
# globally unique, origin-tagged ids. gffcompare classifies every transcript against the
# reference, and consensus_endaware.py collapses the novel candidates: two transcripts merge
# only when their intron chains agree within junction_wobble and their 3' ends within tes_tol,
# so alternative polyA isoforms stay distinct. With consensus_ignore_tss_multiexon (the
# default) the 5' criterion applies to mono-exons only, which have no chain to anchor them.
#

gffcompare = "/sc/arion/projects/ad-omics/data/software/gffcompare-0.12.10.Linux_x86_64/gffcompare"
python_bin = "/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/python"

min_consensus_tools = config.get("min_consensus_tools", 1)   # 1 keeps the union of all three tools
junction_wobble     = config.get("consensus_junction_wobble", 6)
tss_tol             = config.get("consensus_tss_tol", 100)
tes_tol             = config.get("consensus_tes_tol", 100)
ignore_tss_multi    = config.get("consensus_ignore_tss_multiexon", True)   # 5' ends are unreliable on ONT
n_tools = 3


def consensus_query_gtfs(wildcards):
    files = []
    for g in groups:
        files.append(out_folder + "bambu/" + g + "/" + data_code + "_" + g + "_extended_annotations.gtf")
        files.append(out_folder + "isoquant/" + g + "/" + data_code + "_" + g + "/" + data_code + "_" + g + ".extended_annotation.gtf")
        files.append(out_folder + "stringtie3/" + g + "/" + data_code + "_" + g + ".stringtie.gtf")
    return files


def union_gtf_args():
    """--gtf TOOL GROUP PATH triples for build_union.py, matching consensus_query_gtfs."""
    parts = []
    for g in groups:
        parts.append("--gtf bambu %s %sbambu/%s/%s_%s_extended_annotations.gtf" % (g, out_folder, g, data_code, g))
        parts.append("--gtf isoquant %s %sisoquant/%s/%s_%s/%s_%s.extended_annotation.gtf" % (g, out_folder, g, data_code, g, data_code, g))
        parts.append("--gtf stringtie %s %sstringtie3/%s/%s_%s.stringtie.gtf" % (g, out_folder, g, data_code, g))
    return " ".join(parts)


rule consensus_isoforms:
    input:
        gtfs = consensus_query_gtfs,
        ref  = ref_gtf
    output:
        gtf      = out_folder + "consensus/" + data_code + "_merged_consensus.gtf",
        venn_tsv = out_folder + "consensus/" + data_code + "_tool_overlap.tsv",
        stats    = out_folder + "consensus/" + data_code + "_consensus_stats.tsv",
        # kept rather than temp so per-isoform provenance can be regenerated later
        member   = out_folder + "consensus/" + data_code + "_tool_membership.tsv"
    params:
        prefix     = out_folder + "consensus/" + data_code + "_gffcmp",
        union      = out_folder + "consensus/" + data_code + "_union.gtf",
        union_args = union_gtf_args(),
        min_tools  = min_consensus_tools,
        wobble     = junction_wobble,
        tss        = tss_tol,
        tes        = tes_tol,
        ign_tss    = "--ignore-tss-multiexon" if ignore_tss_multi else "",
        exclude    = config.get("novel_exclude_codes", "=,c"),
        python     = python_bin
    shell:
        """
        conda activate isoseq-pipeline
        {params.python} scripts/build_union.py {params.union_args} --reference {input.ref} -o {params.union}
        {gffcompare} -r {input.ref} -o {params.prefix} {params.union}
        {params.python} scripts/consensus_endaware.py \
        --union {params.union} \
        --gffcmp-prefix {params.prefix} \
        --min-tools {params.min_tools} \
        --junction-wobble {params.wobble} \
        --tss-tol {params.tss} \
        --tes-tol {params.tes} {params.ign_tss} \
        --exclude-codes '{params.exclude}' \
        --overlap-out {output.venn_tsv} \
        --membership-out {output.member} \
        --stats-out {output.stats} \
        -o {output.gtf}
        rm -f {params.union} {params.prefix}.combined.gtf {params.prefix}.*.tmap {params.prefix}.*.refmap
        """


# the quantification reference: the full annotation verbatim plus the novel isoforms that
# passed SQANTI, renamed <prefix>_*. this rule is the only place renaming, gene inheritance
# and provenance happen
rule reference_with_novel:
    input:
        merged = out_folder + "consensus/" + data_code + "_merged_consensus.gtf",
        member = out_folder + "consensus/" + data_code + "_tool_membership.tsv",
        sqanti = out_folder + "sqanti/" + data_code + "_filter_sqanti_classification.tsv",
        ref    = ref_gtf
    output:
        gtf  = temp(out_folder + "consensus/" + data_code + ".isoforms.gtf"),
        prov = out_folder + "consensus/" + data_code + "_novel_provenance.tsv"
    params:
        keep    = out_folder + "consensus/" + data_code + "_sqanti_pass_ids.txt",
        prefix  = config.get("novel_prefix", "NOVEL"),
        source  = config.get("novel_source", "ONT"),
        exclude = config.get("novel_exclude_codes", "=,c"),
        script  = "scripts/build_reference_with_novel.py",
        python  = "/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/python"
    shell:
        # column 1 of the sqanti classification is the isoform id
        """
        tail -n +2 {input.sqanti} | cut -f1 > {params.keep}
        {params.python} {params.script} \
        --merged {input.merged} \
        --reference {input.ref} \
        --prefix {params.prefix} \
        --source {params.source} \
        --exclude-codes '{params.exclude}' \
        --keep-list {params.keep} \
        --membership {input.member} \
        --provenance-out {output.prov} \
        -o {output.gtf}
        """


# transcript fasta for downstream quantification. the genome .fai is already present so
# gffread does not need write access to the reference folder
rule reference_fasta:
    input:
        gtf    = out_folder + "consensus/" + data_code + ".isoforms.gtf",
        genome = referenceFa
    output:
        fa = out_folder + "consensus/" + data_code + ".isoforms.fa.gz"
    params:
        tmp = out_folder + "consensus/" + data_code + ".isoforms.fa"
    shell:
        """
        conda activate isoseq-pipeline
        gffread -w {params.tmp} -g {input.genome} {input.gtf}
        gzip -nc {params.tmp} > {output.fa}
        rm -f {params.tmp}
        """


rule compress_reference_gtf:
    input:
        gtf = out_folder + "consensus/" + data_code + ".isoforms.gtf"
    output:
        gz = out_folder + "consensus/" + data_code + ".isoforms.gtf.gz"
    shell:
        "gzip -nc {input.gtf} > {output.gz}"


#####
## Diagnostics
#####

sanity_genes = config.get("sanity_genes", ["STMN2", "UNC13A", "KCNQ2"])

# known cryptic events traced stage by stage. counting novel isoforms is not a valid
# recovery test on its own, since an event annotated in the reference is scored as 0 novel
sanity_events = config.get("sanity_events", [
    "STMN2:chr8:79616822-79617207",
    "UNC13A:chr19:17642414-17642541",
    "KCNQ2:chr20:63439709-63444658:skip",
])

# gffcmp.annotated.gtf and sqanti_pass_ids.txt are side products of earlier rules rather
# than declared outputs, so they are passed as params. ordering is still guaranteed by input.gtf
rule sanity_genes:
    input:
        gtf  = out_folder + "consensus/" + data_code + ".isoforms.gtf",
        ref  = ref_gtf,
        cons = out_folder + "consensus/" + data_code + "_merged_consensus.gtf"
    output:
        tsv    = out_folder + "consensus/" + data_code + "_sanity_genes.tsv",
        events = out_folder + "consensus/" + data_code + "_sanity_genes.tsv.events.tsv"
    params:
        genes  = ",".join(sanity_genes),
        events = ",".join(sanity_events),
        tol    = config.get("sanity_event_tol", 0),
        prefix = config.get("novel_prefix", "NOVEL"),
        source = config.get("novel_source", "ONT"),
        union  = out_folder + "consensus/" + data_code + "_gffcmp.annotated.gtf",
        sqc    = out_folder + "sqanti/" + data_code + "_classification.txt",
        sqp    = out_folder + "consensus/" + data_code + "_sqanti_pass_ids.txt",
        script = "scripts/check_genes.py",
        python = "/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/python"
    shell:
        """
        {params.python} {params.script} \
        --reference {input.ref} \
        --gtf {input.gtf} \
        --genes '{params.genes}' \
        --prefix {params.prefix} \
        --source {params.source} \
        --events '{params.events}' \
        --event-tol {params.tol} \
        --union-gtf {params.union} \
        --consensus-gtf {input.cons} \
        --sqanti-classification {params.sqc} \
        --sqanti-pass-ids {params.sqp} \
        -o {output.tsv}
        """


# novel candidates carrying at least one intron absent from the pooled short-read SJ.out.tab,
# i.e. isoforms that could not have been rebuilt from short reads alone
rule longread_specific:
    input:
        cons = out_folder + "consensus/" + data_code + "_merged_consensus.gtf"
    output:
        tsv = out_folder + "consensus/" + data_code + "_longread_specific.tsv"
    params:
        annotated = out_folder + "consensus/" + data_code + "_gffcmp.annotated.gtf",
        sj_folder = junction_folder,
        min_uniq  = config.get("longread_specific_sj_min_uniq", 1),
        exclude   = config.get("longread_specific_exclude_samples", ""),
        wobble    = config.get("longread_specific_wobble", 2),
        codes     = config.get("novel_exclude_codes", "=,c"),
        cohort    = data_code,
        python    = python_bin
    shell:
        """
        {params.python} scripts/longread_specific.py \
        --annotated {params.annotated} \
        --sj-folder {params.sj_folder} \
        --sj-min-uniq {params.min_uniq} \
        --exclude-samples '{params.exclude}' \
        --exclude-codes '{params.codes}' \
        --wobble {params.wobble} \
        --cohort {params.cohort} \
        -o {output.tsv}
        """


# the per-cohort funnel, from raw per-tool counts through to final novel. column
# definitions are in the script header
rule tracking_numbers:
    input:
        stats  = out_folder + "consensus/" + data_code + "_consensus_stats.tsv",
        member = out_folder + "consensus/" + data_code + "_tool_membership.tsv",
        sqanti = out_folder + "sqanti/" + data_code + "_filter_sqanti_classification.tsv",
        lrspec = out_folder + "consensus/" + data_code + "_longread_specific.tsv",
        ref    = ref_gtf
    output:
        tsv = out_folder + "consensus/" + data_code + "_tracking_numbers.tsv"
    params:
        exclude = config.get("novel_exclude_codes", "=,c"),
        cohort  = data_code,
        python  = python_bin
    shell:
        """
        {params.python} scripts/tracking_numbers.py \
        --stats {input.stats} \
        --membership {input.member} \
        --sqanti-classification {input.sqanti} \
        --longread-specific {input.lrspec} \
        --reference {input.ref} \
        --exclude-codes '{params.exclude}' \
        --cohort {params.cohort} \
        -o {output.tsv}
        """
