# PacBio subreads to FLNC: ccs, then lima to remove primers, then isoseq3 refine
# only active when prep is pacbio and pacbio_input is subreads, inert otherwise.
# the default pacbio_input of flnc reads flnc_bam_path straight from the metadata
# tool versions: ccs 6.0.0, lima 2.2.0, isoseq3 3.4.0

if prep == "pacbio" and pacbio_input == "subreads":

    def subreads_xml(wildcards):
        return metadata_dict[wildcards.sample]["subreads_xml"]

    rule ccs:
        output:
            bam    = out_folder + "{sample}/ccs/{sample}.ccs.bam",
            report = out_folder + "{sample}/ccs/{sample}.ccs.report.txt"
        threads: 16
        run:
            xml = subreads_xml(wildcards)
            shell("ccs -j {threads} %s {output.bam} --report-file {output.report}" % xml)

    # assumes one sample per subreadset. multiplexed cells would need --split-bam-named
    # and a per-barcode demultiplex step
    rule isoseq_lima:
        input:
            bam     = out_folder + "{sample}/ccs/{sample}.ccs.bam",
            primers = pacbio_primers
        output:
            bam = out_folder + "{sample}/lima/{sample}.fl.bam"
        threads: 16
        shell:
            "lima --isoseq --peek-guess -j {threads} {input.bam} {input.primers} {output.bam}"

    rule isoseq_refine:
        input:
            bam     = out_folder + "{sample}/lima/{sample}.fl.bam",
            primers = pacbio_primers
        output:
            bam = out_folder + "{sample}/flnc/{sample}.flnc.bam"
        threads: 16
        shell:
            "isoseq3 refine --require-polya -j {threads} {input.bam} {input.primers} {output.bam}"
