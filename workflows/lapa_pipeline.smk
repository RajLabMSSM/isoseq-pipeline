# lapa_pipeline.smk
# APA (alternative polyadenylation) branch — LAPA 0.0.5 (Mortazavi lab, ENCODE4).
# DOWNSTREAM of the isoform reference: LAPA calls poly(A)-site clusters de novo from the
# LONG-read 3' ends (pooled across all samples in ONE run = joint discovery, the recommended
# mode), annotates them against OUR novel isoform reference (GENCODE + novels), and keeps
# per-sample counts for a differential-APA test between the two groups. Long reads give the
# true 3' cleavage site, so this is the APA layer; short-read Salmon handles expression/DTU.
#
# Pooled discovery is achieved by handing LAPA the PER-SAMPLE bams in one invocation (a CSV
# sample,dataset,path) -- NOT a pre-merged bam, which LAPA neither accepts nor needs: it pools
# internally for cluster calling and keeps per-sample identity for the cross-sample replication
# filter (min_replication_rate 0.95) + the differential. polyA tails were retained upstream
# (pychopper removes primers only, no polyA trim), so LAPA's internal-priming model has real
# tail signal. All intermediates are temp(); the persisted deliverables sit next to the final
# reference GTF in consensus/: <dc>_polyA_clusters.bed and <dc>_apa_differential.tsv.

run_lapa      = config.get("run_lapa", False)
lapa_dir      = out_folder + "lapa/"
lapa_env      = config.get("lapa_env", "lapa_env")
lapa_counting = config.get("lapa_counting_method", "end")   # 'end' (3' end) or 'tail' (needs polyA tail)


if prep == "nanopore_cdna" and run_lapa and len(groups) >= 2:

    rule lapa_bam:
        input:
            cram = out_folder + "{sample}/alignment/{sample}.aligned.cram"
        output:
            bam = temp(lapa_dir + "bams/{sample}.bam"),
            bai = temp(lapa_dir + "bams/{sample}.bam.bai")
        threads: 4
        shell:
            "conda activate {lapa_env}; "
            "samtools view -@ {threads} -b -T {referenceFa} {input.cram} > {output.bam}; "
            "samtools index {output.bam}"

    rule lapa_chrom_sizes:
        input:
            fai = referenceFa + ".fai"
        output:
            temp(lapa_dir + "chrom_sizes.txt")
        shell:
            "cut -f1,2 {input.fai} > {output}"

    rule lapa_annotation:
        input:
            gz = out_folder + "consensus/" + data_code + ".isoforms.gtf.gz"
        output:
            temp(lapa_dir + data_code + "_annotation.gtf")
        shell:
            "gzip -dc {input.gz} > {output}"

    rule lapa_sample_config:
        input:
            bams = expand(lapa_dir + "bams/{sample}.bam", sample=samples),
            bais = expand(lapa_dir + "bams/{sample}.bam.bai", sample=samples)
        output:
            csv = temp(lapa_dir + "sample_config.csv")
        run:
            with open(output.csv, "w") as fh:
                fh.write("sample,dataset,path\n")
                for g in groups:
                    for s in sample_groups[g]:
                        fh.write("%s,%s,%sbams/%s.bam\n" % (s, g, lapa_dir, s))

    # joint (pooled) poly(A)-site discovery + per-sample quantification, annotated on our reference
    rule lapa_run:
        input:
            csv   = lapa_dir + "sample_config.csv",
            bams  = expand(lapa_dir + "bams/{sample}.bam", sample=samples),
            fasta = referenceFa,
            annot = lapa_dir + data_code + "_annotation.gtf",
            chrom = lapa_dir + "chrom_sizes.txt"
        output:
            clusters = out_folder + "consensus/" + data_code + "_polyA_clusters.bed",
            resdir   = temp(directory(lapa_dir + "results"))
        params:
            counting = lapa_counting,
            env      = lapa_env
        threads: 8
        shell:
            "conda activate {params.env}; "
            "lapa --alignment {input.csv} --fasta {input.fasta} "
            "--annotation {input.annot} --chrom_sizes {input.chrom} "
            "--counting_method {params.counting} --output_dir {output.resdir}; "
            "cp {output.resdir}/polyA_clusters.bed {output.clusters}"

    # differential APA between the two groups from LAPA's per-sample cluster counts
    # (per-gene poly(A)-site usage shift; Fisher + delta-PAU + BH). Runs in lapa_env
    # (lapa_diff.py imports lapa). config CSV carries the sample->group mapping.
    rule lapa_differential:
        input:
            resdir   = lapa_dir + "results",
            csv      = lapa_dir + "sample_config.csv",
            clusters = out_folder + "consensus/" + data_code + "_polyA_clusters.bed"
        output:
            diff = out_folder + "consensus/" + data_code + "_apa_differential.tsv"
        params:
            group_a = groups[0],
            group_b = groups[1],
            script  = "scripts/lapa_diff.py",
            env     = lapa_env
        shell:
            "conda activate {params.env}; "
            "python {params.script} --lapa-dir {input.resdir} --config {input.csv} "
            "--group-a {params.group_a} --group-b {params.group_b} -o {output.diff}"
