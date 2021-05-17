# convert CSV TPM matrix to TPM 
# set transcripts to rownames and write out without a column name for transcripts, just for samples
# also gsub the ends of the sample IDs
library(optparse)

option_list <- list(
    make_option(c('--input', '-i'), help = "the name of the input file", default = ""),
    make_option(c('--output', '-o'), help = "the name of the output file", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

output <- opt$output
input <- opt$input

stopifnot(file.exists(input) )

suppressPackageStartupMessages(library(tidyverse))

df <- read_csv(input) %>% as.data.frame()

rownames(df) <- df$id

df$id <- NULL

names(df) <- gsub(".aligned.md" ,"", names(df) )

write.table(df, file = output,  sep = "\t", quote = FALSE, row.names = TRUE) 
