library(optparse)
options(BIOCONDUCTOR_CONFIG_FILE = "/sc/arion/projects/ad-omics/data/software/bambu/config.yaml")

arguments <- parse_args(OptionParser(option_list = list(
    make_option(c("--anno", "-a"), default = ""),
    make_option(c("--prefix", "-p"), default = ""),
    make_option(c("--fasta", "-f"), default = ""),
    make_option(c("--cores", "-c"), default = 1),
    make_option(c("--NDR", "-n"), default = NULL)
)), positional_arguments = TRUE)
opt <- arguments$opt

library(bambu)
load(opt$anno)
# a FaFile, not a path: a path string makes bambu query BiocManager over the network
genome <- Rsamtools::FaFile(opt$fasta)
# all of the group's bams together: joint discovery
res <- bambu(reads = arguments$args, annotations = anno, genome = genome,
             NDR = as.numeric(opt$NDR), ncore = as.numeric(opt$cores))
writeBambuOutput(se = res, path = dirname(opt$prefix), prefix = paste0(basename(opt$prefix), "_"))
save(res, file = paste0(opt$prefix, "_bambu.RData"))
