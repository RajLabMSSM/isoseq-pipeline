library(tidyverse)
library(rtracklayer)

# parse arguments
# create outFile names from input names


df <- vroom::vroom("all_samples/cupcake/all_samples.demux_fl_count.csv")

df <- column_to_rownames(df, var = "id")

# remove any transcript that is not present in at least 80% of samples
clean <- df[ rowSums(df != 0) >= floor(0.8 * ncol(df) ), ]

clean <- rownames_to_column(clean, var = "id")

readr::write_csv(clean, path = "all_samples/cupcake_filtered/all_samples.demux_fl_count_filtered.csv")

# filter GFF
gff_file <- "all_samples/cupcake/all_samples.cupcake.collapsed.gff"

gff <- import.gff(gff_file)

clean_gff <- gff[gff$transcript_id %in% clean$id]

export(clean_gff, con = "all_samples/cupcake_filtered/all_samples.cupcake.collapsed.filtered.gtf", format = "GTF")
