library(tidyverse)

# collates lists of read lengths into a single table for plotting
# give option of downsampling

library(optparse)

set.seed(1234)
option_list <- list(
        make_option(c('--output', "-o"), help='', default = "example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

output <- opt$output

mapped_files <- list.files(pattern = "\\.mapped.readlengths.txt", recursive = TRUE, full.names = TRUE)

unmapped_files <- list.files(pattern = "\\.unmapped.readlengths.txt", recursive = TRUE, full.names = TRUE)

message( "* ", length(mapped_files), " samples found with mapped read lengths")

# read in lists of read lengths
# subsample to 100,000
# make nice dataframe for plotting
read_lengths <- function(file, n = 1e5){
    sample_id <- gsub("\\.(unmapped|mapped).readlengths.txt" , "", basename(file) )
    type <- ifelse(grepl("unmapped", file), "unmapped", "mapped")

    x <- as.numeric(readLines(file))
    if(n < Inf){
        if( n < length(x) ){ 
            x <- x[sample(n, replace=FALSE)]
        }
    }
    df <- data.frame(read_length = x, sample = sample_id, type = type)
    return(df)
}

mapped_df <- purrr::map_df( mapped_files, read_lengths, n = 1e5)
unmapped_df <- purrr::map_df( unmapped_files, read_lengths, n = 1e5)

all_df <- dplyr::bind_rows(mapped_df, unmapped_df)

readr::write_tsv(all_df, output)


    

