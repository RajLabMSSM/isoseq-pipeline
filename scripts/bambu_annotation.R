# turn the reference gtf into bambu's annotation object once, so each group's run can load it (Chen et al. 2023), optparse means you don't have to hardcode anything
library(optparse)
# read the input gtf path and the output .RData path
opt <- parse_args(OptionParser(option_list = list(
    make_option(c("--input", "-i"), default = ""),
    make_option(c("--output", "-o"), default = "")
)))
suppressPackageStartupMessages(library(bambu))
# the object must be called anno, because bambu_run.R uses that name after load()
anno <- prepareAnnotations(opt$input)
save(anno, file = opt$output)
