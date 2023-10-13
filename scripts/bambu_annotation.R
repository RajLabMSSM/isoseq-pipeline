message(" * BAMBU - annotation")
library(optparse)

option_list <- list(
    make_option(c('--output', '-o' ), help='the annotation RData file produced by bambu', default = ""),
    make_option(c('--input', '-i'), help = "the GTF annotation file", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


input <- opt$input
output <- opt$output

input_file <- opt$input
output_file <- opt$output

message(" * preparing annotations for bambu")

message(paste0(" * using ", input_file) )

suppressPackageStartupMessages(library(bambu))
anno <- prepareAnnotations(input_file)

message(paste0(" * writing to ", output_file))

save(anno, file = output_file)

