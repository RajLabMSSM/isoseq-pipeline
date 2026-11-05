import glob
import os

# the rules filter definition kept in the pipeline repo
sqanti_rules = os.path.join(workflow.basedir, "misc", "sqanti_rules.json")
# STAR junctions for SQANTI's short-read junction coverage. exclude_samples are dropped (tdpkd: TDP-5
# aberrant differentiation, CTRL-4 low-depth outlier); the sample id is the file name before the
# first dot (<sample>.SJ.out.tab)
sqanti_sj =[f for f in sorted(glob.glob(os.path.join(junction_folder, "*.SJ.out.tab")))
             if os.path.basename(f).split(".")[0] not in excluded_samples]


# the novel consensus models (class codes = and c left out) under their merged ids, so SQANTI's pass
# list maps straight back onto merged_consensus
rule novel_consensus_gtf:
    input:
        merged = cons + "_merged_consensus.gtf",
        ref    = ref_gtf,
        # the script is an input, so editing it re-runs this step
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


# -c gives every junction its short-read support (min_sample_cov) for the rules filter.
# the genePreds are removed first, sqanti reuses a stale one otherwise
rule sqanti_qc:
    input:
        gtf = out_folder + "sqanti/" + data_code + "_novel_consensus.gtf",
        ref = ref_gtf,
        fa  = genome,
        sj  = sqanti_sj
    output:
        classification = out_folder + "sqanti/" + data_code + "_classification.txt",
        corrected_gtf  = out_folder + "sqanti/" + data_code + "_corrected.gtf",
        corrected_fa   = out_folder + "sqanti/" + data_code + "_corrected.fasta"
    params:
        outdir = out_folder + "sqanti/",
        sj     = ",".join(sqanti_sj)
    threads: 12
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


# rules: FSM intra-priming only; the rest also all-canonical + min_sample_cov >= 2. RTS_stage is left out, it drops the UNC13A CE.
# -e (--filter_mono_exonic) also drops every single-exon model
rule sqanti_filter:
    input:
        classification = out_folder + "sqanti/" + data_code + "_classification.txt",
        fasta          = out_folder + "sqanti/" + data_code + "_corrected.fasta",
        gtf            = out_folder + "sqanti/" + data_code + "_corrected.gtf"
    output:
        classification = out_folder + "sqanti/" + data_code + "_RulesFilter_classification.txt",
        pass_ids       = out_folder + "sqanti/" + data_code + "_pass_isoforms.txt",
        gtf            = out_folder + "sqanti/" + data_code + ".filtered.gtf",
        fasta          = out_folder + "sqanti/" + data_code + ".filtered.fasta"
    params:
        outdir = out_folder + "sqanti/"
    threads: 4
    shell:
        """
        conda activate sqanti3_v6
        sqanti3_filter.py rules \\
        --sqanti_class {input.classification} \\
        --filter_gtf {input.gtf} \\
        --filter_isoforms {input.fasta} \\
        -j {sqanti_rules} \\
        -e \\
        -d {params.outdir} \\
        -o {data_code} \\
        -c {threads}
        """
