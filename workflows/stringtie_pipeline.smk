stringtie = "/sc/arion/projects/ad-omics/data/software/stringtie-3.0.3.Linux_x86_64/stringtie"


rule stringtie_pool_short_bams:
    input:
        bams = lambda wc: [metadata_dict[s]["short_read_bam_path"] for s in sample_groups[wc.group]]
    output:
        bam = temp(out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.short.bam"),
        bai = temp(out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.short.bam.bai")
    threads: 4
    shell:
        """
        samtools merge -f -@ {threads} --reference {genome} {output.bam} {input.bams}
        samtools index {output.bam}
        """


# hybrid mode: --mix takes the short-read bam first. -t is off, with it the HYOU1 exon-10 skip is lost
rule run_stringtie:
    input:
        gtf       = ref_gtf,
        long_bam  = lambda wc: pooled_bam(wc.group),
        long_bai  = lambda wc: pooled_bam(wc.group) + ".bai",
        short_bam = lambda wc: short_bam(wc.group)
    output:
        gtf = out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf"
    threads: 4
    shell:
        "{stringtie} -p {threads} --mix -G {input.gtf} -o {output.gtf} {input.short_bam} {input.long_bam}"
