# the pooled stage after discovery, included by cross_cohort/Snakefile, so the cohorts meet once on their raw assembler outputs
# every rule lists its script as an input, so editing a script re-runs that rule
cons = OUT + "consensus/" + NAME
sq   = OUT + "sqanti/" + NAME


# taking a union of all 12 assembler gtfs (3 tools x 2 groups x 2 cohorts) using isomatch merge, which keeps any model from any one tool
# gencode transcripts, = and c models, single-exon and unstranded models are removed before the merge
rule consensus_isoforms:
    input:
        gtfs         = [path for _tool, _cohort, _group, path in TOOL_GTFS],
        ref          = REF_GTF,
        union_script = SCRIPTS + "/build_union.py",
        merge_script = SCRIPTS + "/consensus_isomatch.py"
    output:
        merged     = cons + "_merged_consensus.gtf",
        overlap    = cons + "_tool_overlap.tsv",
        stats      = cons + "_consensus_stats.tsv",
        membership = cons + "_tool_membership.tsv",
        # declared because rule sanity_genes reads it
        annotated  = cons + "_gffcmp.annotated.gtf"
    params:
        gffcmp_prefix = cons + "_gffcmp",
        union         = cons + "_union.gtf",
        union_args    = " ".join("--gtf %s %s %s %s" % entry for entry in TOOL_GTFS),
        workdir       = OUT + "consensus/isomatch_tmp"
    shell:
        """
        python {input.union_script} {params.union_args} --reference {input.ref} -o {params.union}
        {GFFCOMPARE} -r {input.ref} -o {params.gffcmp_prefix} {params.union}
        python {input.merge_script} \
        --union {params.union} \
        --gffcmp-prefix {params.gffcmp_prefix} \
        --reference {input.ref} \
        --genome {GENOME} \
        --exclude-codes '=,c' \
        --isomatch {ISOMATCH} \
        --gffcompare {GFFCOMPARE} \
        --workdir {params.workdir} \
        --overlap-out {output.overlap} \
        --membership-out {output.membership} \
        --stats-out {output.stats} \
        -o {output.merged}
        rm -rf {params.union} {params.gffcmp_prefix}.combined.gtf {params.gffcmp_prefix}.*.tmap {params.gffcmp_prefix}.*.refmap {params.workdir}
        """


# the reference verbatim plus the SQANTI-passed novels, renamed XNOVEL_*, gene from SQANTI's associated_gene
# each novel keeps a cohorts attribute (sun, tdpkd or both) and one row in the provenance table
rule reference_with_novel:
    input:
        merged      = cons + "_merged_consensus.gtf",
        membership  = cons + "_tool_membership.tsv",
        sqanti_pass = sq + "_pass_isoforms.txt",
        sqanti      = sq + "_classification.txt",
        ref         = REF_GTF,
        script      = SCRIPTS + "/build_reference_with_novel.py"
    output:
        gtf        = temp(OUT + NAME + ".gtf"),
        provenance = OUT + "novel_isoform_discovery_by_tool_cohort_and_condition.tsv"
    shell:
        """
        python {input.script} \
        --merged {input.merged} \
        --reference {input.ref} \
        --sqanti-classification {input.sqanti} \
        --exclude-codes '=,c' \
        --keep-list {input.sqanti_pass} \
        --membership {input.membership} \
        --prefix XNOVEL \
        --provenance-out {output.provenance} \
        -o {output.gtf}
        """


rule compress_reference_gtf:
    input:
        OUT + NAME + ".gtf"
    output:
        OUT + NAME + ".gtf.gz"
    shell:
        "gzip -nc {input} > {output}"


# transcript sequences of the final gtf, made with the salmon-pipeline gffread wrapper
rule reference_fasta:
    input:
        gtf    = OUT + NAME + ".gtf.gz",
        script = FASTA_SH
    output:
        OUT + NAME + ".fa.gz"
    shell:
        "bash {input.script} {input.gtf} {GENOME} {output}"


# traces each positive-control junction through union -> consensus -> SQANTI -> final gtf
# STMN2, KCNQ2, NUP93 and XPNPEP1 are annotated in GENCODE v50 and score 0 novel, so assembled_by and assembled_in show whether our data built them
rule sanity_genes:
    input:
        gtf            = OUT + NAME + ".gtf",
        ref            = REF_GTF,
        merged         = cons + "_merged_consensus.gtf",
        union          = cons + "_gffcmp.annotated.gtf",
        sqanti_pass    = sq + "_pass_isoforms.txt",
        sqanti_reasons = sq + "_filtering_reasons.txt",
        script         = SCRIPTS + "/check_genes.py"
    output:
        tsv    = OUT + NAME + "_sanity_genes.tsv",
        events = OUT + NAME + "_sanity_genes.tsv.events.tsv"
    shell:
        """
        python {input.script} \
        --reference {input.ref} \
        --gtf {input.gtf} \
        --genes '{SANITY_GENES}' \
        --prefix XNOVEL \
        --source ONT \
        --events '{SANITY_EVENTS}' \
        --event-tol 0 \
        --union-gtf {input.union} \
        --consensus-gtf {input.merged} \
        --sqanti-pass-ids {input.sqanti_pass} \
        --sqanti-reasons {input.sqanti_reasons} \
        -o {output.tsv}
        """


# the build's funnel numbers (raw -> consensus -> SQANTI -> final), every metric prefixed pooled_
rule isoform_filtering_summary_pooled:
    input:
        stats       = cons + "_consensus_stats.tsv",
        membership  = cons + "_tool_membership.tsv",
        sqanti      = sq + "_RulesFilter_classification.txt",
        sqanti_pass = sq + "_pass_isoforms.txt",
        ref         = REF_GTF,
        script      = SCRIPTS + "/tracking_numbers.py"
    output:
        OUT + "isoform_filtering_summary.tsv"
    shell:
        """
        python {input.script} \
        --stats {input.stats} \
        --membership {input.membership} \
        --sqanti-classification {input.sqanti} \
        --pass-ids {input.sqanti_pass} \
        --reference {input.ref} \
        --cohort pooled \
        --metric-prefix pooled_ \
        -o {output}
        """


# one cohort's share of the funnel numbers
rule isoform_filtering_summary_cohort:
    input:
        stats       = cons + "_consensus_stats.tsv",
        membership  = cons + "_tool_membership.tsv",
        sqanti      = sq + "_RulesFilter_classification.txt",
        sqanti_pass = sq + "_pass_isoforms.txt",
        ref         = REF_GTF,
        script      = SCRIPTS + "/tracking_numbers.py"
    output:
        OUT + "isoform_filtering_summary_{cohort}.tsv"
    wildcard_constraints:
        cohort = "|".join(COHORT_NAMES)
    shell:
        """
        python {input.script} \
        --stats {input.stats} \
        --membership {input.membership} \
        --sqanti-classification {input.sqanti} \
        --pass-ids {input.sqanti_pass} \
        --reference {input.ref} \
        --cohort {wildcards.cohort} \
        --only {wildcards.cohort} \
        -o {output}
        """


# writes the funnel numbers into the project tracking table, skipped by builds without a tracking_numbers key (the T2T build)
rule upsert_tracking:
    input:
        pooled = OUT + "isoform_filtering_summary.tsv",
        cohort = expand(OUT + "isoform_filtering_summary_{cohort}.tsv", cohort=COHORT_NAMES),
        script = SCRIPTS + "/upsert_tracking_numbers.py"
    output:
        touch(OUT + "isoform_filtering_summary.done")
    params:
        cohort = " ".join("%s=%sisoform_filtering_summary_%s.tsv" % (c, OUT, c) for c in COHORT_NAMES)
    shell:
        """
        python {input.script} {config[tracking_numbers]} {config[run_code]} {input.pooled} {params.cohort}
        """
