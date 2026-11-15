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
# The _corrected.gtf file is declared because rule sqanti_browser draws the models from it.
# SQANTI also writes a _corrected.fasta file. Nothing reads it, so it is not declared as an output.
rule sqanti_qc:
    input:
        gtf = out_folder + "sqanti/" + data_code + "_novel_consensus.gtf",
        ref = ref_gtf,
        fa  = genome,
        sj  = sqanti_sj
    output:
        classification = out_folder + "sqanti/" + data_code + "_classification.txt",
        corrected_gtf  = out_folder + "sqanti/" + data_code + "_corrected.gtf"
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


# Build a UCSC Genome Browser track hub from the SQANTI calls with SQANTI-browser (Paniagua et al. 2026, bioRxiv,
# doi:10.64898/2026.05.25.727625).
# The main track colours each model by its SQANTI structural category, and one hidden track per category follows it.
# The input is the rules-filter classification, which is the QC classification plus the filter_result column.
# SQANTI-browser keeps every classification column as a field of the model's details page. A model kept by the
# filter therefore reads Isoform there and a removed one reads Artifact, so both can be inspected side by side.
# Only novel models reach SQANTI, so the hub holds no GENCODE transcript. The browser's own GENCODE track serves
# for comparison.
# The CAGE track is the same peak file that sqanti_qc used for dist_to_CAGE_peak, so it shows that evidence.
# Chromosome sizes are taken from the genome index. Without them the tool guesses each length from the furthest feature.
# GENCODE names unplaced scaffolds like GL000195.1, while UCSC hg38 names them like chrUn_GL000195v1. Models on these
# scaffolds stay in the bigBed files, but UCSC hg38 cannot display them (83 of 47,410 models in tdpkd v9).
# The hub name is the data code, so hubs of different cohorts can be loaded together and told apart.
# hubCheck then checks the hub and opens every bigBed file, and the job fails if it reports a problem.
# hubCheck gets hub.txt by its absolute path, because it looks for the track files in the wrong folder from a relative one.
# UCSC reads a hub only over public HTTP, so the hub folder must be copied to a public web server by hand (see README).
rule sqanti_browser:
    input:
        gtf            = out_folder + "sqanti/" + data_code + "_corrected.gtf",
        classification = out_folder + "sqanti/" + data_code + "_RulesFilter_classification.txt",
        fai            = genome + ".fai",
        cage           = config['sqanti_cage_peak']
    output:
        hub         = directory(out_folder + "sqanti/" + data_code + "_ucsc_hub"),
        chrom_sizes = temp(out_folder + "sqanti/" + data_code + "_chrom.sizes")
    params:
        hub_txt = os.path.abspath(out_folder + "sqanti/" + data_code + "_ucsc_hub/hub.txt")
    threads: 1
    shell:
        """
        conda activate sqanti-browser
        cut -f1,2 {input.fai} > {output.chrom_sizes}
        python -m sqanti_browser \\
        --gtf {input.gtf} \\
        --classification {input.classification} \\
        --CAGE-peak {input.cage} \\
        --chrom-sizes {output.chrom_sizes} \\
        --genome hg38 \\
        --hub-name {data_code} \\
        --output {output.hub}
        hubCheck {params.hub_txt}
        """
