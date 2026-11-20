# The pooled stage after discovery. It is included by cross_cohort/Snakefile, which defines OUT, NAME, REF_GTF,
# GENOME, TOOL_GTFS, COHORT_NAMES, SCRIPTS, the tool paths and the sanity events. The per-cohort Snakefile no longer
# includes it, so the cohorts meet once, on their raw assembler outputs.
cons = OUT + "consensus/" + NAME


# Combine the assembler outputs of every cohort, remove the transcripts that carry a GENCODE id, and merge the
# rest with tama_merge (Kuo et al. 2020). Each of the three tools ran on the pooled reads of each group of each
# cohort, so the inputs are 3 tools x 2 groups x 2 cohorts = 12 GTFs.
# Removal comes before the merge. Bambu and IsoQuant pass the whole of GENCODE through, and TAMA would otherwise
# fold a novel model that differs from a known transcript only at its 5' end or within 100 bp at its 3' end into
# that known transcript.
# gffcompare 0.12.10 (Pertea and Pertea 2020) then classifies the remaining union against GENCODE. Models of class
# = (same intron chain as a GENCODE transcript, which also catches a match that no tool labelled) and c (contained
# in one) are left out of the merge, as are single-exon and unstranded models.
# The merged models are classified against GENCODE again with gffcompare, because merging can change a model's
# chain or ends. TAMA keeps every model built by any one source, so this is a union with no tool vote.
rule consensus_isoforms:
    input:
        gtfs = [path for _tool, _cohort, _group, path in TOOL_GTFS],
        ref  = REF_GTF,
        # Both scripts are inputs, so editing either one rebuilds the consensus.
        union_script = SCRIPTS + "/build_union.py",
        tama_script  = SCRIPTS + "/consensus_tama.py"
    output:
        gtf      = cons + "_merged_consensus.gtf",
        venn_tsv = cons + "_tool_overlap.tsv",
        stats    = cons + "_consensus_stats.tsv",
        member   = cons + "_tool_membership.tsv",
        # gffcompare's classification of the union is kept and declared, because rule sanity_genes reads it.
        annotated = cons + "_gffcmp.annotated.gtf"
    params:
        prefix     = cons + "_gffcmp",
        union      = cons + "_union.gtf",
        # One --gtf <tool> <cohort> <group> <path> argument per assembler output.
        union_args = " ".join("--gtf %s %s %s %s" % entry for entry in TOOL_GTFS),
        workdir    = OUT + "consensus/tama_tmp"
    shell:
        """
        python {input.union_script} {params.union_args} --reference {input.ref} -o {params.union}
        {GFFCOMPARE} -r {input.ref} -o {params.prefix} {params.union}
        python {input.tama_script} \
        --union {params.union} \
        --gffcmp-prefix {params.prefix} \
        --reference {input.ref} \
        --junction-wobble {config[wobble]} \
        --tes-tol {config[tes_tol]} \
        --exclude-codes '=,c' \
        --tama-merge {TAMA} \
        --python2 {PYTHON2} \
        --gffcompare {GFFCOMPARE} \
        --workdir {params.workdir} \
        --overlap-out {output.venn_tsv} \
        --membership-out {output.member} \
        --stats-out {output.stats} \
        -o {output.gtf}
        rm -rf {params.union} {params.prefix}.combined.gtf {params.prefix}.*.tmap {params.prefix}.*.refmap {params.workdir}
        """


# Write GENCODE unchanged followed by the novel models that passed SQANTI, renamed XNOVEL_NNNNNNNN.
# Each novel model's gene is SQANTI's associated_gene, and a SQANTI novel gene becomes a new XNOVELG_ locus. One
# SQANTI run covers both cohorts, so every locus id is unique in the build.
# Each novel transcript keeps the cohorts attribute (sun, tdpkd or sun,tdpkd) of the cohorts whose assemblies built it.
# The provenance table holds one row per novel model with its gene, class code, tools, cohorts and groups.
rule reference_with_novel:
    input:
        merged = cons + "_merged_consensus.gtf",
        member = cons + "_tool_membership.tsv",
        keep   = OUT + "sqanti/" + NAME + "_pass_isoforms.txt",
        sqanti = OUT + "sqanti/" + NAME + "_classification.txt",
        ref    = REF_GTF,
        # The script is an input, so editing it re-runs this step.
        script = SCRIPTS + "/build_reference_with_novel.py"
    output:
        gtf  = temp(OUT + NAME + ".gtf"),
        prov = OUT + "cross_cohort_novel_provenance.tsv"
    shell:
        """
        python {input.script} \
        --merged {input.merged} \
        --reference {input.ref} \
        --sqanti-classification {input.sqanti} \
        --exclude-codes '=,c' \
        --keep-list {input.keep} \
        --membership {input.member} \
        --prefix XNOVEL \
        --provenance-out {output.prov} \
        -o {output.gtf}
        """


# Compress the final GTF.
rule compress_reference_gtf:
    input:
        gtf = OUT + NAME + ".gtf"
    output:
        gz = OUT + NAME + ".gtf.gz"
    shell:
        "gzip -nc {input.gtf} > {output.gz}"


# Write the transcript sequences of the final GTF with gffread. This FASTA is the input of the salmon index.
rule reference_fasta:
    input:
        gtf = OUT + NAME + ".gtf.gz",
        # The salmon-pipeline gffread wrapper is an input, so editing it rebuilds the FASTA.
        script = FASTA_SH
    output:
        fa = OUT + NAME + ".fa.gz"
    shell:
        "bash {input.script} {input.gtf} {GENOME} {output.fa}"


# Find the isoforms of the final GTF that carry each positive-control junction, and trace each event through
# union, consensus, SQANTI and final GTF to show where it is kept or lost, and which tools and cohorts built it.
# STMN2, KCNQ2, NUP93 and XPNPEP1 are in GENCODE v50, so they always match in the final GTF. For those,
# assembled_by and assembled_in show whether our data built the event. An annotated event scores 0 novel, so novel
# counts alone are not a recovery test.
rule sanity_genes:
    input:
        gtf  = OUT + NAME + ".gtf",
        ref  = REF_GTF,
        cons = cons + "_merged_consensus.gtf",
        sqp  = OUT + "sqanti/" + NAME + "_pass_isoforms.txt",
        sqr  = OUT + "sqanti/" + NAME + "_filtering_reasons.txt",
        union = cons + "_gffcmp.annotated.gtf",
        # The script is an input, so editing it re-runs this step.
        script = SCRIPTS + "/check_genes.py"
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
        --consensus-gtf {input.cons} \
        --sqanti-pass-ids {input.sqp} \
        --sqanti-reasons {input.sqr} \
        -o {output.tsv}
        """


# Write the funnel numbers (raw, consensus, SQANTI, final) of the whole build. Every metric name starts with
# pooled_, so these rows stay apart from each cohort's own rows in the project tracking table.
rule tracking_numbers_pooled:
    input:
        stats  = cons + "_consensus_stats.tsv",
        member = cons + "_tool_membership.tsv",
        sqanti = OUT + "sqanti/" + NAME + "_RulesFilter_classification.txt",
        keep   = OUT + "sqanti/" + NAME + "_pass_isoforms.txt",
        ref    = REF_GTF,
        # The script is an input, so editing it re-runs this step.
        script = SCRIPTS + "/tracking_numbers.py"
    output:
        tsv = OUT + "tracking_numbers_final_numbers.tsv"
    shell:
        """
        python {input.script} \
        --stats {input.stats} \
        --membership {input.member} \
        --sqanti-classification {input.sqanti} \
        --pass-ids {input.keep} \
        --reference {input.ref} \
        --cohort pooled \
        --metric-prefix pooled_ \
        -o {output.tsv}
        """


# Write one cohort's share of the funnel: its own assembler counts, and the merged models that its assemblies built.
rule tracking_numbers_cohort:
    input:
        stats  = cons + "_consensus_stats.tsv",
        member = cons + "_tool_membership.tsv",
        sqanti = OUT + "sqanti/" + NAME + "_RulesFilter_classification.txt",
        keep   = OUT + "sqanti/" + NAME + "_pass_isoforms.txt",
        ref    = REF_GTF,
        # The script is an input, so editing it re-runs this step.
        script = SCRIPTS + "/tracking_numbers.py"
    output:
        tsv = OUT + "tracking_numbers_{cohort}.tsv"
    wildcard_constraints:
        cohort = "|".join(COHORT_NAMES)
    shell:
        """
        python {input.script} \
        --stats {input.stats} \
        --membership {input.member} \
        --sqanti-classification {input.sqanti} \
        --pass-ids {input.keep} \
        --reference {input.ref} \
        --cohort {wildcards.cohort} \
        --only {wildcards.cohort} \
        -o {output.tsv}
        """


# Write this build's funnel numbers into the project tracking table. Each cohort gets its own share plus the pooled
# rows. A build without a tracking_numbers key in its config, such as the T2T comparison build, skips this step, so
# it cannot replace the hg38 rows.
rule upsert_tracking:
    input:
        pooled = OUT + "tracking_numbers_final_numbers.tsv",
        cohort = expand(OUT + "tracking_numbers_{cohort}.tsv", cohort=COHORT_NAMES),
        # The script is an input, so editing it re-runs this step.
        script = SCRIPTS + "/upsert_tracking_numbers.py"
    output:
        done = touch(OUT + "tracking_numbers.done")
    params:
        cohort = " ".join("%s=%stracking_numbers_%s.tsv" % (c, OUT, c) for c in COHORT_NAMES)
    shell:
        """
        python {input.script} {config[tracking_numbers]} {config[run_code]} {input.pooled} {params.cohort}
        """
