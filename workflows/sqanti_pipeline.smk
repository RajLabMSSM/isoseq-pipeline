# SQANTI3 6.0.1 QC and rules filter, applied to the novel consensus isoforms only
# this is the precision filter for the catalogue. surviving ids feed reference_with_novel
# QC runs in the sqanti3_v6 env, the filter script needs R from bambu_env

sqanti_threads = 12


# drop transfrags that already match the reference, keeping the original TCONS ids so
# the surviving ids map straight back onto merged_consensus
rule novel_consensus_gtf:
    input:
        merged = out_folder + "consensus/" + data_code + "_merged_consensus.gtf",
        ref    = ref_gtf
    output:
        gtf = out_folder + "sqanti/" + data_code + "_novel_consensus.gtf"
    params:
        exclude = config.get("novel_exclude_codes", "=,c"),
        script  = "scripts/build_reference_with_novel.py",
        python  = "/sc/arion/projects/als-omics/conda/envs/isoseq-pipeline/bin/python"
    shell:
        """
        {params.python} {params.script} \
        --merged {input.merged} \
        --reference {input.ref} \
        --exclude-codes '{params.exclude}' \
        --novel-only --no-rename \
        -o {output.gtf}
        """


# no short reads are supplied, so this is a sequence-only run
# both genePreds are removed first, sqanti reuses a stale one otherwise
rule sqanti_qc:
    input:
        gtf = out_folder + "sqanti/" + data_code + "_novel_consensus.gtf",
        ref = ref_gtf,
        fa  = referenceFa
    output:
        classification = out_folder + "sqanti/" + data_code + "_classification.txt",
        corrected_gtf  = out_folder + "sqanti/" + data_code + "_corrected.gtf",
        corrected_fa   = out_folder + "sqanti/" + data_code + "_corrected.fasta"
    params:
        outdir = out_folder + "sqanti/",
        prefix = data_code
    threads: sqanti_threads
    shell:
        """
        conda activate sqanti3_v6
        rm -f {params.outdir}{params.prefix}_corrected.genePred {params.outdir}refAnnotation_{params.prefix}.genePred
        sqanti3_qc.py \
        --isoforms {input.gtf} \
        --refGTF {input.ref} \
        --refFasta {input.fa} \
        -d {params.outdir} \
        -o {params.prefix} \
        -t {threads}
        """


rule sqanti_filter:
    input:
        classification = out_folder + "sqanti/" + data_code + "_classification.txt",
        fasta          = out_folder + "sqanti/" + data_code + "_corrected.fasta",
        gff            = out_folder + "sqanti/" + data_code + "_corrected.gtf"
    output:
        gtf            = out_folder + "sqanti/" + data_code + "_filter_sqanti.cds.gtf",
        classification = out_folder + "sqanti/" + data_code + "_filter_sqanti_classification.tsv",
        fasta          = out_folder + "sqanti/" + data_code + "_filter_sqanti.fasta"
    params:
        script     = "scripts/filter_sqanti.R",
        out_prefix = out_folder + "sqanti/" + data_code
    shell:
        """
        conda activate bambu_env
        Rscript {params.script} \
        --sqanti {input.classification} \
        --fasta {input.fasta} \
        --gff {input.gff} \
        --output {params.out_prefix}
        """
