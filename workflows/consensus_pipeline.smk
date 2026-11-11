gffcompare = "/sc/arion/projects/ad-omics/data/software/gffcompare-0.12.10.Linux_x86_64/gffcompare"     #Compare assembled transcripts with reference annotations
tama_merge = "/sc/arion/projects/ad-omics/data/software/tama/tama_merge.py"                             # Merge transcript models with similar structures
python2    = "/sc/arion/projects/als-omics/conda/envs/tama_env/bin/python"                              # Python interpreter used to run TAMA
cons       = out_folder + "consensus/" + data_code

sanity_gene_list  = "UNC13A,STMN2,KCNQ2,HYOU1,NUP93,XPNPEP1"
sanity_event_list = ("UNC13A:chr19:17641557-17642413+17642542-17642844," # because both junctions are on either end of the cryptic exon
                     "STMN2:chr8:79611215-79616821,"
                     "KCNQ2:chr20:63439709-63444658,"
                     "HYOU1:chr11:119052173-119052636,"
                     "NUP93:chr16:56839070-56841704,"
                     "XPNPEP1:chr10:109884149-109888048")


# the three assembler gtfs of one group, by tool (so tool_gtfs("control") would give the 2 control GTFs)
def tool_gtfs(g):
    return {"bambu":     out_folder + "bambu/%s/%s_%s_extended_annotations.gtf" % (g, data_code, g),
            "isoquant":  out_folder + "isoquant/%s/%s_%s/%s_%s.extended_annotation.gtf" % (g, data_code, g, data_code, g),
            "stringtie": out_folder + "stringtie3/%s/%s_%s.stringtie.gtf" % (g, data_code, g)}


# taking a union using tama_merge, which keeps any model from any one tool
rule consensus_isoforms:
    input:
        gtfs = [p for g in groups for p in tool_gtfs(g).values()],
        ref  = ref_gtf,
        # both scripts are inputs (consensus_tama.py runs tama_merge), so editing either rebuilds the consensus
        union_script = "scripts/build_union.py",
        tama_script  = "scripts/consensus_tama.py"
    output:
        gtf      = cons + "_merged_consensus.gtf",
        venn_tsv = cons + "_tool_overlap.tsv",
        stats    = cons + "_consensus_stats.tsv",
        member   = cons + "_tool_membership.tsv"
    params:
        prefix     = cons + "_gffcmp",
        union      = cons + "_union.gtf",
        # one --gtf <tool> <group> <path> argument per tool and group
        union_args = " ".join("--gtf %s %s %s" % (t, g, p) for g in groups for t, p in tool_gtfs(g).items()),
        workdir    = out_folder + "consensus/tama_tmp"
    shell:
        """
        python {input.union_script} {params.union_args} --reference {input.ref} -o {params.union}
        {gffcompare} -r {input.ref} -o {params.prefix} {params.union}
        python {input.tama_script} \
        --union {params.union} \
        --gffcmp-prefix {params.prefix} \
        --reference {input.ref} \
        --junction-wobble {config[consensus_junction_wobble]} \
        --tes-tol {config[consensus_tes_tol]} \
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
        ref    = ref_gtf,
        # the script is an input, so editing it re-runs this step
        script = "scripts/build_reference_with_novel.py"
    output:
        gtf  = temp(cons + ".isoforms.gtf"),
        prov = cons + "_novel_provenance.tsv"
    shell:
        """
        python {input.script} \
        --merged {input.merged} \
        --reference {input.ref} \
        --sqanti-classification {input.sqanti} \
        --exclude-codes '=,c' \
        --keep-list {input.keep} \
        --membership {input.member} \
        --provenance-out {output.prov} \
        -o {output.gtf}
        """


# transcript sequences of the cohort's final gtf
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


# finds the isoforms after SQANTI (GENCODE plus the novels that passed) that carry each control junction, and traces
# each event through union -> consensus -> SQANTI -> final gtf to show where it is kept or lost.
# STMN2, KCNQ2, NUP93 and XPNPEP1 are in GENCODE v50, so they always match in the final gtf: for those,
# assembled_by shows whether our data built the event. An annotated event scores 0 novel, so novel counts alone
# are not a recovery test
rule sanity_genes:
    input:
        gtf  = cons + ".isoforms.gtf",
        ref  = ref_gtf,
        cons = cons + "_merged_consensus.gtf",
        sqp  = out_folder + "sqanti/" + data_code + "_pass_isoforms.txt",
        sqr  = out_folder + "sqanti/" + data_code + "_filtering_reasons.txt",
        # the script is an input, so editing it re-runs this step
        script = "scripts/check_genes.py"
    output:
        tsv    = cons + "_sanity_genes.tsv",
        events = cons + "_sanity_genes.tsv.events.tsv"
    params:
        union = cons + "_gffcmp.annotated.gtf"
    shell:
        """
        python {input.script} \
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
        --sqanti-reasons {input.sqr} \
        -o {output.tsv}
        """


# the cohort's funnel numbers (raw -> consensus -> SQANTI -> final) as metric<TAB>value
rule tracking_numbers:
    input:
        stats  = cons + "_consensus_stats.tsv",
        member = cons + "_tool_membership.tsv",
        sqanti = out_folder + "sqanti/" + data_code + "_RulesFilter_classification.txt",
        keep   = out_folder + "sqanti/" + data_code + "_pass_isoforms.txt",
        ref    = ref_gtf,
        # the script is an input, so editing it re-runs this step
        script = "scripts/tracking_numbers.py"
    output:
        tsv = cons + "_tracking_numbers.tsv"
    shell:
        """
        python {input.script} \
        --stats {input.stats} \
        --membership {input.member} \
        --sqanti-classification {input.sqanti} \
        --pass-ids {input.keep} \
        --reference {input.ref} \
        --cohort {data_code} \
        -o {output.tsv}
        """
