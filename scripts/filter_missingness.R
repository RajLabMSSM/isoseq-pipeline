suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))
library(optparse)

option_list <- list(
    make_option(c('--counts', '-c' ), help='the count matrix', default = ""),
    make_option(c('--gff', '-g' ), help = "the GTF or GFF", default = ""),
    make_option(c('--prefix', '-p'), help = "the prefix to the output files", default = ""),
    make_option(c('--min_samples', '-s'), help = "the minimum proportion or minimum number of samples", default = 2),
    make_option(c('--min_reads', '-r'), help = "the minimum read count per sample", default = 2),
    make_option(c('--remove_monoexons', '-m'), help = "whether to remove mono-exonic transcripts", action="store_true", default = FALSE)
)

#df <- vroom::vroom("all_samples/cupcake/all_samples.demux_fl_count.csv")
#gff_file <- "all_samples/cupcake/all_samples.cupcake.collapsed.gff"

#input <- "all_samples/cupcake/all_samples"
#output <- "all_samples/cupcake_filtered/all_samples"

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

counts_file <- opt$counts
gff_file <- opt$gff
output <- opt$prefix
min_samples <- as.numeric(opt$min_samples)
min_reads <- as.numeric(opt$min_reads)
remove_monoexons <- opt$remove_monoexons

stopifnot(min_samples > 0 & min_samples < ncol(df) )
stopifnot(min_reads > 0)

# make sure input and output are valid 
#stopifnot(dir.exists(basename(counts_file) ) )
#stopifnot(dir.exists(basename(output) ) )

# parse arguments
#counts_file <- paste0(input, ".demux_fl_count.csv" )
#gff_file <- paste0(input, ".cupcake.collapsed.gff")

stopifnot(file.exists(counts_file) )
stopifnot(file.exists(gff_file) )

counts_out <- paste0(output, "_filtered_strong.csv")
gff_out <- paste0(output, "_filtered_strong.gtf")

message(" * reading in demultiplexed transcript counts")

df <- vroom::vroom(counts_file)

# detect type of abundance file from file name
    if ( grepl("talon", counts_file ) ){
    df <- column_to_rownames(df, var = "annot_transcript_id")

    # remove all other columns
    df <- df[,11:ncol(df) ]

}

if( grepl("collapsed", counts_file) ){

df <- column_to_rownames(df, var = "id")
}
# add parameter min_samples for min_samples filtering
# if min_samples between 0 and <1 - keep transcripts present in min_samples * total sample size
# if min_samples = 1 or greater - keep transcripts present in at least min_samples samples
if( min_samples > 0 & min_samples < 1){
    message( paste0(" * Keeping only transcripts with support in at least ", min_reads, " reads in at least ", min_samples * 100, "% of samples"))
    clean <- df[ rowSums(df >= min_reads) >= floor(min_reads * ncol(df) ), ]
}

if( min_samples >= 1 ){
    message( paste0(" * Keeping transcripts with support in at least ",  min_reads, " reads in at least ", min_samples, " samples"))
    clean <- df[ rowSums(df >= min_reads) >= min_samples , ]
}


message( paste0( " * ", nrow(clean) , " out of ", nrow(df), " transcripts kept!" ) )


# filter GFF - using default gffread cleaned GFF3 files now
message(" * importing GFF")
gff_file_ext <- tools::file_ext(gff_file)
gff <- import(gff_file, format = gff_file_ext)

save.image("debug.RData")

clean_gff <- gff[gff$transcript_id %in% row.names(clean)]
stopifnot(length(clean_gff) > 0 )

if( remove_monoexons == TRUE){
    message(" * removing mono-exonic transcripts")
    tx_tally <- clean_gff %>% as.data.frame()  %>% group_by(transcript_id) %>% summarise(n_tx = n() )
    
    multiexons <- filter(tx_tally, n_tx > 2)

    clean_gff <- clean_gff[ clean_gff$transcript_id %in% multiexons$transcript_id ]
    message(" * kept ", nrow(multiexons), " multi-exonic transcripts!" )

    clean <- clean[ multiexons$transcript_id,]

    stopifnot(nrow(clean) == length(unique(clean_gff$transcript_id) ) ) 
    
}    


message(" * writing filtered counts")

clean <- rownames_to_column(clean, var = "id")
readr::write_csv(clean, path = counts_out)


message(" * writing filtered GFF")
export(clean_gff, con = gff_out, format = "GTF")

