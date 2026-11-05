# run bambu transcript discovery and quantification on one group of long-read bams (Chen et al. 2023)
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
# a FaFile, not a path: a path string makes bambu query BiocManager over the network
genome <- Rsamtools::FaFile(opt$fasta)
# all of the group's bams together: joint discovery; NDR = novel discovery rate, the expected
# share of false models among the novel ones bambu keeps (Chen et al. 2023)
res <- bambu(reads = arguments$args, annotations = anno, genome = genome,
             NDR = as.numeric(opt$NDR), ncore = as.numeric(opt$cores))
# write the gtf and count tables, each file name starting with "<prefix>_"
writeBambuOutput(se = res, path = dirname(opt$prefix), prefix = paste0(basename(opt$prefix), "_"))
# keep the full result; the object inside the .RData is named res, so the name stays
save(res, file = paste0(opt$prefix, "_bambu.RData"))
