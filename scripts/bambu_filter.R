library(optparse)

option_list <- list(
    make_option(c('--matrix', '-m' ), help='the count matrix', default = ""),
    make_option(c('--gff', '-g' ), help = "the GTF or GFF", default = ""),
    make_option(c('--prefix', '-p'), help = "the prefix to the output files", default = ""),
    make_option(c('--min_samples', '-s'), help = "the minimum proportion or minimum number of samples", default = 2),
    make_option(c('--min_reads', '-r'), help = "the minimum read count per sample", default = 2),
    make_option(c('--counts', '-c'), help = "filter on raw counts", action="store_true", default = FALSE),
    make_option(c('--tpm', '-t'), help = "filter on transcripts per million (TPM)", action="store_true", default = FALSE),
    make_option(c('--remove_monoexons', '-e'), help = "whether to remove mono-exonic transcripts", action="store_true", default = FALSE)
)

#df <- vroom::vroom("all_samples/cupcake/all_samples.demux_fl_count.csv")
#gff_file <- "all_samples/cupcake/all_samples.cupcake.collapsed.gff"

#input <- "all_samples/cupcake/all_samples"
#output <- "all_samples/cupcake_filtered/all_samples"

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

counts_file <- opt$matrix
gff_file <- opt$gff
output <- opt$prefix
min_samples <- as.numeric(opt$min_samples)
min_reads <- as.numeric(opt$min_reads)

# mode - filter either by TPM or raw counts
stopifnot( opt$tpm == TRUE | opt$counts == TRUE )

if(opt$tpm == TRUE){
    mode <- "tpm"
}
if(opt$counts == TRUE){
    mode <- "counts"
}

stopifnot(mode %in% c("counts", "tpm" ) )

remove_monoexons <- opt$remove_monoexons

stopifnot(min_samples > 0 & min_samples < ncol(df) )
#stopifnot(min_reads > 0)

# make sure input and output are valid 
#stopifnot(dir.exists(basename(counts_file) ) )
#stopifnot(dir.exists(basename(output) ) )

# parse arguments
#counts_file <- paste0(input, ".demux_fl_count.csv" )
#gff_file <- paste0(input, ".cupcake.collapsed.gff")

stopifnot(file.exists(counts_file) )
stopifnot(file.exists(gff_file) )

counts_out <- paste0(output, "_miss_counts.csv")
tpm_out <- paste0(output, "_miss_tpm.csv")
gff_out <- paste0(output, "_miss.gtf")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))

message(" * reading in transcript counts")

df <- vroom::vroom(counts_file)

# detect type of abundance file from file name
if ( grepl("talon", counts_file ) ){
    message(" * TALON input ")
    df <- column_to_rownames(df, var = "annot_transcript_id")

    # remove all other columns
    df <- df[,11:ncol(df) ]

}
# CUPCAKE
if( grepl("collapsed", counts_file) ){
    message(" * Cupcake input")
    df <- column_to_rownames(df, var = "id")
}
# BAMBU
if( grepl("counts_transcript.txt", counts_file) ){
    message(" * Bambu input")
    df$GENEID <- NULL
    df <- column_to_rownames(df, var = "TXNAME")
    # round to integers
    df <- round(df)
}
#ISOQUANT
if( grepl("transcript_model_grouped_counts", counts_file) ){
    message(" * isoquant input")
    df <- column_to_rownames(df, var = "#feature_id")
}


# add parameter min_samples for min_samples filtering
# if min_samples between 0 and <1 - keep transcripts present in min_samples * total sample size
# if min_samples = 1 or greater - keep transcripts present in at least min_samples samples

# filter on TPM rather than counts
tpm_df <- edgeR::cpm(df)


# deprecated for now
if(FALSE){
if( min_samples > 0 & min_samples < 1){
    
    message( paste0(" * Keeping only transcripts with support in at least ", min_reads, " ", mode,  " in at least ", min_samples * 100, "% of samples"))
    
    if( mode == "counts"){
        clean <- df[ rowSums(df >= min_reads) >= floor(min_reads * ncol(df) ), ]
    }
    if( mode == "tpm"){
        clean <- df[ rowSums(tpm_df >= min_reads) >= floor(min_reads * ncol(tpm_df) ), ]
    }
}

if( min_samples >= 1 ){
    message( paste0(" * Keeping transcripts with support in at least ",  min_reads, " ", mode, " in at least ", min_samples, " samples"))
    if( mode == "counts" ){
        clean <- df[ rowSums(df >= min_reads) >= min_samples , ]
    }
    if( mode == "tpm" ){
        clean <- df[ rowSums(tpm_df >= min_reads) >= min_samples , ]
    }
}
}


# new for 2023: keep all transcripts found at least once
# keep all annotated transcripts
# only keep novel if seen at least twice
clean <- df[
    ( grepl("ENST", row.names(df) ) & rowSums(df > 0) > 0 ) |
    ( !grepl("ENST",row.names(df) ) & rowSums(df > 0) > 1 )
,]

message( paste0( " * ", nrow(clean) , " out of ", nrow(df), " transcripts kept!" ) )



message( paste0( " * ", nrow(clean) , " out of ", nrow(df), " transcripts kept!" ) )


# filter GFF - using default gffread cleaned GFF3 files now
message(" * importing GFF")
gff_file_ext <- tools::file_ext(gff_file)
gff <- import(gff_file, format = gff_file_ext)

# remove duplicate entries - bug with isoquant?
gff <- unique(gff)

#save.image("debug.RData")

clean_gff <- gff[gff$transcript_id %in% row.names(clean)]
stopifnot(length(clean_gff) > 0 )

# get rid of mono-exons if user requests
if( remove_monoexons == TRUE){
    message(" * removing mono-exonic transcripts")
    tx_tally <- clean_gff %>% as.data.frame()  %>% group_by(transcript_id) %>% summarise(n_tx = n() )
    
    multiexons <- filter(tx_tally, n_tx > 2)

    clean_gff <- clean_gff[ clean_gff$transcript_id %in% multiexons$transcript_id ]
    message(" * kept ", nrow(multiexons), " multi-exonic transcripts!" )

    clean <- clean[ multiexons$transcript_id,]

    stopifnot(nrow(clean) == length(unique(clean_gff$transcript_id) ) ) 
    
}    


# remove transcripts without assigned strand - found in small number in Bambu GTF
message( paste0(" * removing ", length( clean_gff[strand(clean_gff) == "*" ] ), " transcripts with ambiguous stranding" ) )
clean_gff <- clean_gff[ strand(clean_gff) != "*" ]

clean <- clean[ unique(clean_gff$transcript_id), ]

message(" * writing filtered counts and TPM")

clean <- rownames_to_column(clean, var = "id")
readr::write_csv(clean, path = counts_out)


tpm_clean <- as.data.frame(tpm_df[clean$id,])
tpm_clean <- rownames_to_column(tpm_clean, var = "id")

readr::write_csv(tpm_clean, path = tpm_out)

message(" * writing filtered GFF")
export(clean_gff, con = gff_out, format = "GTF")

