import glob
import os

# Path to the SQANTI rules filter definition, which is kept in the pipeline repo.
sqanti_rules = os.path.join(workflow.basedir, "misc", "sqanti_rules.json")
# Collect the STAR junction files that SQANTI uses to measure short-read support for each junction.
# Samples listed in exclude_samples are left out. For tdpkd these are TDP-5 (aberrant differentiation)
# and CTRL-4 (a low-depth outlier).
# The sample id is the part of the file name before the first dot, as in <sample>.SJ.out.tab.
sqanti_sj =[f for f in sorted(glob.glob(os.path.join(junction_folder, "*.SJ.out.tab")))
             if os.path.basename(f).split(".")[0] not in excluded_samples]


# Write only the novel consensus models as SQANTI's input.
# Models with gffcompare class code = (same intron chain as a reference transcript) or c (contained in one)
# are left out.
# The models keep their merged ids, so SQANTI's pass list maps straight back onto the merged consensus.
rule novel_consensus_gtf:
    input:
        merged = cons + "_merged_consensus.gtf",
        ref    = ref_gtf,
        # The script is listed as an input, so editing it re-runs this step.
        script = "scripts/build_reference_with_novel.py"
    output:
        gtf = out_folder + "sqanti/" + data_code + "_novel_consensus.gtf"
    shell:
        """
        python {input.script} \
        --merged {input.merged} \
        --reference {input.ref} \
        --exclude-codes '=,c' \
        --novel-only --no-rename \
        -o {output.gtf}
        """


# Run SQANTI3 quality control on the novel models.
# The -c option gives every junction its short-read support (min_sample_cov), which the rules filter uses.
# The CAGE peaks and the polyA motif list only add information to the QC report. No filter rule reads them.
# Old genePred files are removed first, because SQANTI otherwise reuses a stale one.
# One thread is enough. The input is a GTF, so SQANTI aligns nothing and has no work to split.
# SQANTI also writes _corrected.gtf and .fasta files. Nothing reads them, so they are not declared as outputs.
rule sqanti_qc:
    input:
        gtf = out_folder + "sqanti/" + data_code + "_novel_consensus.gtf",
        ref = ref_gtf,
        fa  = genome,
        sj  = sqanti_sj
    output:
        classification = out_folder + "sqanti/" + data_code + "_classification.txt"
    params:
        outdir = out_folder + "sqanti/",
        sj     = ",".join(sqanti_sj)
    threads: 1
    shell:
        """
        conda activate sqanti3_v6
        rm -f {params.outdir}{data_code}_corrected.genePred {params.outdir}refAnnotation_{data_code}.genePred
        sqanti3_qc.py \\
        --isoforms {input.gtf} \\
        --refGTF {input.ref} \\
        --refFasta {input.fa} \\
        -c '{params.sj}' \\
        --CAGE_peak {config[sqanti_cage_peak]} --polyA_motif_list {config[sqanti_polya_motifs]} \\
        -d {params.outdir} \\
        -o {data_code} \\
        -t {threads}
        """


# Filter the SQANTI calls with the rules in sqanti_rules.json.
# Full splice matches (FSM) are checked for intra-priming only.
# Every other model must also have only canonical junctions, each supported by short reads in at least
# 2 samples (min_sample_cov >= 2).
# The RTS_stage rule is left out, because it removes the UNC13A cryptic exon.
# The -e option (--filter_mono_exonic) is a safeguard. The consensus step already drops single-exon models,
# so it removes nothing today, but it would catch them if that step changed.
# The FSM rule is the same kind of safeguard. No FSM model reaches this step, because models with class
# code = are left out above.
# The next steps only need the pass list, so no filtered GTF or FASTA is requested.
# The filtering reasons file is read by rule sanity_genes.
rule sqanti_filter:
    input:
        classification = out_folder + "sqanti/" + data_code + "_classification.txt"
    output:
        classification = out_folder + "sqanti/" + data_code + "_RulesFilter_classification.txt",
        pass_ids       = out_folder + "sqanti/" + data_code + "_pass_isoforms.txt",
        reasons        = out_folder + "sqanti/" + data_code + "_filtering_reasons.txt"
    params:
        outdir = out_folder + "sqanti/"
    shell:
        """
        conda activate sqanti3_v6
        sqanti3_filter.py rules \\
        --sqanti_class {input.classification} \\
        -j {sqanti_rules} \\
        -e \\
        -d {params.outdir} \\
        -o {data_code}
        """
