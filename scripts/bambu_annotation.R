library(optparse)
opt <- parse_args(OptionParser(option_list = list(
    make_option(c("--input", "-i"), default = ""),
    make_option(c("--output", "-o"), default = "")
)))
suppressPackageStartupMessages(library(bambu))
anno <- prepareAnnotations(opt$input)
save(anno, file = opt$output)
