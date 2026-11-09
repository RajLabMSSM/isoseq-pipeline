# isoquant discovery on the pooled long-read bam of one group (Prjibelski et al. 2023)
rule run_isoquant:
    input:
        bams = lambda wc: pooled_bam(wc.group),
        bai  = lambda wc: pooled_bam(wc.group) + ".bai",
        ref  = genome,
        gtf  = ref_gtf
    output:
        gtf = temp(out_folder + "isoquant/{group}/" + data_code + "_{group}/" + data_code + "_{group}.extended_annotation.gtf")
    params:
        outdir = out_folder + "isoquant/{group}/",
        prefix = data_code + "_{group}"
    threads: 8
# no short reads: --illumina_bam only reaches reads outside annotated genes (alignment_processor.py process_intergenic). 
# --clean_start = starts without reusing previous IsoQuant run state
# --large_output = disables large optional output files
# --model_construction_strategy	= controls how transcript models are constructed
# after the run the converted gene database (*.db) is deleted to save space
    shell:
        """
        conda activate isoquant313
        isoquant --data_type nanopore \
        --bam {input.bams} \
        --reference {input.ref} \
        --complete_genedb --genedb {input.gtf} --clean_start --large_output none \
        --model_construction_strategy {config[isoquant_model_strategy]} \
        --output {params.outdir} --prefix {params.prefix} \
        --threads {threads}
        rm -f {params.outdir}*.db
        """
