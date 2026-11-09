stringtie = "/sc/arion/projects/ad-omics/data/software/stringtie-3.0.3.Linux_x86_64/stringtie"

# pool the short-read bams of one group for stringtie --mix. a library in
# exclude_samples stays out of the pool, even when its long reads are used (tdpkd CTRL4: long reads in,
# short reads out), so the same short-read libraries feed the hints, this pool and the SQANTI coverage
rule stringtie_pool_short_bams:
    input:
        bams = lambda wc: [metadata_dict[s]["short_read_bam_path"] for s in sample_groups[wc.group]
                           if os.path.basename(metadata_dict[s]["short_read_bam_path"]).split(".")[0] not in excluded_samples]
    output:
        bam = temp(out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.short.bam"),
        bai = temp(out_folder + "stringtie3/{group}/" + data_code + "_{group}.pooled.short.bam.bai")
    threads: 4
    shell:
        """
        samtools merge -f -@ {threads} --reference {genome} {output.bam} {input.bams}
        samtools index {output.bam}
        """


# hybrid assembly (Shumate et al. 2022): --mix reads the short-read bam first, the long-read bam second; -G guides it
# with GENCODE. in this mode a novel model needs coverage 4.75 aka the average number of reads per base of a trasncript (the -s default), 
# and it can rest on short reads alone.
# i'm not using -t so that i can see the HYOU1 exon-10 skip; using -t turns trimming OFF for the ends of a transcript where coverage falls,
# so i'm not sure why it messes up HYOU1 despite great coverage in 1340 short reads and 10 long reads. oh well
# --rf is off even though the short reads are stranded: it also sets the strand of unspliced long reads, wrongly
rule run_stringtie:
    input:
        gtf       = ref_gtf,
        long_bam  = lambda wc: pooled_bam(wc.group),
        short_bam = lambda wc: short_bam(wc.group)
    output:
        gtf = out_folder + "stringtie3/{group}/" + data_code + "_{group}.stringtie.gtf"
    threads: 4
    shell:
        "{stringtie} -p {threads} --mix -G {input.gtf} -o {output.gtf} {input.short_bam} {input.long_bam}"
