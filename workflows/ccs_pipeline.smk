# For taking raw PacBio data and making FLNC BAM files

# CCS
rule CCS:
    input:
        raw_folder + "{sample}/{sample}.subreadset.xml"
    output:
        bam =  "{sample}/isoseq3-ccs/{sample}.ccs.bam",
        report =  "{sample}/isoseq3-ccs/{sample}.ccs.report.txt"
    shell:
        "ccs -j 0 {input} {output.bam} --report-file {output.report}"

# primer removal and demultiplexing
rule isoseq_lima:
    input:
         "{sample}/isoseq3-ccs/{sample}.ccs.bam",
    output:
         "{sample}/isoseq3-lima/{sample}.fl.bam"
    shell:
        "lima --isoseq --different --min-passes 1 --split-bam-named --dump-clips --dump-removed -j 0 {input} {primers} {output}"

# trimming of polya tails and removal of concatemers
rule isoseq_refine:
    input:
         "{sample}/isoseq3-lima/{sample}.fl.bam"
    output:
         "{sample}/isoseq3-refine/{sample}.flnc.bam"
    shell:
        "isoseq3 refine --require-polya {input} {primers} {output}"

