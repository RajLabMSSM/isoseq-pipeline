# run bambu transcript discovery and quantification on one group of long-read bams (Chen et al. 2023), optparse means you don't have to hardcode anything
library(optparse)
# point Bioconductor at the shared bambu config so it does not look for a per-user one
options(BIOCONDUCTOR_CONFIG_FILE = "/sc/arion/projects/ad-omics/data/software/bambu/config.yaml")

# read the options; the bams come in as positional arguments
arguments <- parse_args(OptionParser(option_list = list(
    make_option(c("--anno", "-a"), default = ""),
    make_option(c("--prefix", "-p"), default = ""),
    make_option(c("--fasta", "-f"), default = ""),
    make_option(c("--cores", "-c"), default = 1),
    make_option(c("--NDR", "-n"), default = NULL)
)), positional_arguments = TRUE)
opt <- arguments$opt

library(bambu)
# the prepared annotation .RData holds one object called anno (from prepareAnnotations)
load(opt$anno)
# a FaFile, not a path: this bambu install stops with an error on a path string
genome <- Rsamtools::FaFile(opt$fasta)
# no --NDR means NULL, which makes bambu use its own recommended NDR
ndr <- if (is.null(opt$NDR)) NULL else as.numeric(opt$NDR)
# one pooled bam per group, so pooling has already happened. NDR = novel discovery rate: the share of unannotated
# models among all candidates ranked at or above a model's score (Chen et al. 2023). with one bam, ncore has no effect
res <- bambu(reads = arguments$args, annotations = anno, genome = genome,
             NDR = ndr, ncore = as.numeric(opt$cores))
# write the gtf and count tables, each file name starting with "<prefix>_"
writeBambuOutput(se = res, path = dirname(opt$prefix), prefix = paste0(basename(opt$prefix), "_"))
# keep the full result; the object inside the .RData is named res, so the name stays
save(res, file = paste0(opt$prefix, "_bambu.RData"))
