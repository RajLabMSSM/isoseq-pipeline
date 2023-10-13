library(readr)
library(purrr)

options <- commandArgs(trailingOnly = TRUE)

outFile <- options[1]

files <- list.files(pattern = "prepend.cluster_report.csv", recursive = TRUE, full.names = TRUE)

all <- map_df(files, read_tsv, col_types = "ccc")

write_tsv(all, path = outFile)
