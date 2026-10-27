gffcompare = "/sc/arion/projects/ad-omics/data/software/gffcompare-0.12.10.Linux_x86_64/gffcompare"
tama_merge = "/sc/arion/projects/ad-omics/data/software/tama/tama_merge.py"
python2    = "/sc/arion/projects/als-omics/conda/envs/tama_env/bin/python"
cons       = out_folder + "consensus/" + data_code
# the two TDP-43 cryptic exons used as positive controls for recovery (hg38, 1-based exon coordinates)
sanity_gene_list  = "STMN2,UNC13A"
sanity_event_list = "STMN2:chr8:79616822-79617207,UNC13A:chr19:17642414-17642541"


def tool_gtfs(g):
    return {"bambu":     out_folder + "bambu/%s/%s_%s_extended_annotations.gtf" % (g, data_code, g),
            "isoquant":  out_folder + "isoquant/%s/%s_%s/%s_%s.extended_annotation.gtf" % (g, data_code, g, data_code, g),
            "stringtie": out_folder + "stringtie3/%s/%s_%s.stringtie.gtf" % (g, data_code, g)}


# every tool x group gtf, unioned; tama_merge keeps any model from any one tool (no vote)
rule consensus_isoforms:
    input:
        gtfs = [p for g in groups for p in tool_gtfs(g).values()],
        ref  = ref_gtf
    output:
        gtf      = cons + "_merged_consensus.gtf",
        venn_tsv = cons + "_tool_overlap.tsv",
        stats    = cons + "_consensus_stats.tsv",
        member   = cons + "_tool_membership.tsv"
    params:
        prefix     = cons + "_gffcmp",
        union      = cons + "_union.gtf",
        union_args = " ".join("--gtf %s %s %s" % (t, g, p) for g in groups for t, p in tool_gtfs(g).items()),
        workdir    = out_folder + "consensus/tama_tmp"
    shell:
        """
        python scripts/build_union.py {params.union_args} --reference {input.ref} -o {params.union}
        {gffcompare} -r {input.ref} -o {params.prefix} {params.union}
        python scripts/consensus_tama.py \
        --union {params.union} \
        --gffcmp-prefix {params.prefix} \
        --reference {input.ref} \
        --junction-wobble {config[consensus_junction_wobble]} \
        --tes-tol {config[consensus_tes_tol]} \
        --ends common_ends \
        --exclude-codes '=,c' \
        --tama-merge {tama_merge} \
        --python2 {python2} \
        --gffcompare {gffcompare} \
        --workdir {params.workdir} \
        --overlap-out {output.venn_tsv} \
        --membership-out {output.member} \
        --stats-out {output.stats} \
        -o {output.gtf}
        rm -rf {params.union} {params.prefix}.combined.gtf {params.prefix}.*.tmap {params.prefix}.*.refmap {params.workdir}
        """


# the reference verbatim plus the SQANTI-passed novels, renamed NOVEL_*, gene from SQANTI's associated_gene
rule reference_with_novel:
    input:
        merged = cons + "_merged_consensus.gtf",
        member = cons + "_tool_membership.tsv",
        keep   = out_folder + "sqanti/" + data_code + "_pass_isoforms.txt",
        sqanti = out_folder + "sqanti/" + data_code + "_classification.txt",
        ref    = ref_gtf
    output:
        gtf  = temp(cons + ".isoforms.gtf"),
        prov = cons + "_novel_provenance.tsv"
    shell:
        """
        python scripts/build_reference_with_novel.py \
        --merged {input.merged} \
        --reference {input.ref} \
        --sqanti-classification {input.sqanti} \
        --prefix NOVEL \
        --source ONT \
        --exclude-codes '=,c' \
        --keep-list {input.keep} \
        --membership {input.member} \
        --provenance-out {output.prov} \
        -o {output.gtf}
        """


rule reference_fasta:
    input:
        gtf    = cons + ".isoforms.gtf",
        genome = genome
    output:
        fa = cons + ".isoforms.fa.gz"
    params:
        tmp = cons + ".isoforms.fa"
    shell:
        """
        gffread -w {params.tmp} -g {input.genome} {input.gtf}
        gzip -nc {params.tmp} > {output.fa}
        rm -f {params.tmp}
        """


rule compress_reference_gtf:
    input:
        gtf = cons + ".isoforms.gtf"
    output:
        gz = cons + ".isoforms.gtf.gz"
    shell:
        "gzip -nc {input.gtf} > {output.gz}"


# traces the STMN2 and UNC13A cryptic exons through union -> consensus -> SQANTI -> final gtf, to show
# where a known event is kept or lost. An annotated event scores 0 novel, so novel counts alone are not a
# recovery test
rule sanity_genes:
    input:
        gtf  = cons + ".isoforms.gtf",
        ref  = ref_gtf,
        cons = cons + "_merged_consensus.gtf",
        sqp  = out_folder + "sqanti/" + data_code + "_pass_isoforms.txt"
    output:
        tsv    = cons + "_sanity_genes.tsv",
        events = cons + "_sanity_genes.tsv.events.tsv"
    params:
        union = cons + "_gffcmp.annotated.gtf",
        sqr   = out_folder + "sqanti/" + data_code + "_filtering_reasons.txt"
    shell:
        """
        python scripts/check_genes.py \
        --reference {input.ref} \
        --gtf {input.gtf} \
        --genes '{sanity_gene_list}' \
        --prefix NOVEL \
        --source ONT \
        --events '{sanity_event_list}' \
        --event-tol 0 \
        --union-gtf {params.union} \
        --consensus-gtf {input.cons} \
        --sqanti-pass-ids {input.sqp} \
        --sqanti-reasons {params.sqr} \
        -o {output.tsv}
        """


rule tracking_numbers:
    input:
        stats  = cons + "_consensus_stats.tsv",
        member = cons + "_tool_membership.tsv",
        sqanti = out_folder + "sqanti/" + data_code + "_RulesFilter_classification.txt",
        keep   = out_folder + "sqanti/" + data_code + "_pass_isoforms.txt",
        ref    = ref_gtf
    output:
        tsv = cons + "_tracking_numbers.tsv"
    shell:
        """
        python scripts/tracking_numbers.py \
        --stats {input.stats} \
        --membership {input.member} \
        --sqanti-classification {input.sqanti} \
        --pass-ids {input.keep} \
        --reference {input.ref} \
        --cohort {data_code} \
        -o {output.tsv}
        """
