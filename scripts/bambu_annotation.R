# turn GENCODE gtf into bambu's annotation object once, so each group's run can load it (Chen et al. 2023) - as opposed to making it every time
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
