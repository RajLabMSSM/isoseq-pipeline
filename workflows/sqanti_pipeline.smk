# SQANTI3 on the pooled novel models (Pardo-Palacios et al. 2024). This file is included by cross_cohort/Snakefile,
# which defines OUT, NAME, REF_GTF, GENOME, SQANTI_SJ and the shared genome settings.
sq = OUT + "sqanti/" + NAME


# Write only the novel consensus models as SQANTI's input.
# Models with gffcompare class code = (same intron chain as a reference transcript) or c (contained in one) after
# the merge are left out.
# The models keep their merged ids, so SQANTI's pass list maps straight back onto the merged consensus.
rule novel_consensus_gtf:
    input:
        merged = cons + "_merged_consensus.gtf",
        ref    = REF_GTF,
        # The script is listed as an input, so editing it re-runs this step.
        script = SCRIPTS + "/build_reference_with_novel.py"
    output:
        gtf = sq + "_novel_consensus.gtf"
    shell:
        """
        python {input.script} \
        --merged {input.merged} \
        --reference {input.ref} \
        --exclude-codes '=,c' \
        --novel-only --no-rename \
        -o {output.gtf}
        """


# Run SQANTI3 quality control on the pooled novel models.
# The -c option gives every junction its short-read support. SQANTI counts, per junction, the STAR SJ.out.tab files
# with at least one uniquely mapped read across it (sample_with_cov). A model's min_sample_cov is the lowest of these
# counts over its junctions. The files are those of both cohorts (10 tdpkd and 8 sun libraries), so a model built in
# one cohort can reach min_sample_cov >= 2 on short reads of the other cohort.
# The CAGE peaks and the polyA motif list only add information to the QC report. No filter rule reads them.
# Old genePred files are removed first, because SQANTI otherwise reuses a stale one.
# One thread is enough. The input is a GTF, so SQANTI aligns nothing and has no work to split.
# SQANTI also writes _corrected.gtf and _corrected.fasta files. Nothing reads them, so they are not declared as outputs.
rule sqanti_qc:
    input:
        gtf = sq + "_novel_consensus.gtf",
        ref = REF_GTF,
        fa  = GENOME,
        sj  = SQANTI_SJ
    output:
        classification = sq + "_classification.txt"
    params:
        outdir = OUT + "sqanti/",
        sj     = ",".join(SQANTI_SJ)
    threads: 1
    shell:
        """
        conda activate sqanti3_v6
        rm -f {params.outdir}{NAME}_corrected.genePred {params.outdir}refAnnotation_{NAME}.genePred
        sqanti3_qc.py \\
        --isoforms {input.gtf} \\
        --refGTF {input.ref} \\
        --refFasta {input.fa} \\
        -c '{params.sj}' \\
        --CAGE_peak {SHARED[sqanti_cage_peak]} --polyA_motif_list {SHARED[sqanti_polya_motifs]} \\
        -d {params.outdir} \\
        -o {NAME} \\
        -t {threads}
        """


# filter with SQANTI3's default rules (Pardo-Palacios et al. 2024) plus a third rest rule in misc/sqanti_filter_rules.json
# the third rule keeps a model flagged for RT template switching only when every junction has >= 3 short reads (min_cov)
# -e is only a safeguard, because single-exon models never reach SQANTI
# the filtering reasons file is read by rule sanity_genes
rule sqanti_filter:
    input:
        classification = sq + "_classification.txt",
        rules          = ISOSEQ_DIR + "/misc/sqanti_filter_rules.json"
    output:
        classification = sq + "_RulesFilter_classification.txt",
        pass_ids       = sq + "_pass_isoforms.txt",
        reasons        = sq + "_filtering_reasons.txt"
    params:
        outdir = OUT + "sqanti/"
    shell:
        """
        conda activate sqanti3_v6
        sqanti3_filter.py rules \\
        --sqanti_class {input.classification} \\
        -j {input.rules} \\
        -e \\
        -d {params.outdir} \\
        -o {NAME}
        """
