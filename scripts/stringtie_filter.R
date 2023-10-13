# assemble stringtie matrix and filter
library(optparse)

option_list <- list(
    make_option(c('--inFolder', '-i' ), help='the input directory', default = ""),
    make_option(c('--runCode', '-r' ), help="the run code", default = "stringtie"),
    make_option(c('--gff', '-g' ), help = "the GTF or GFF", default = ""),
    make_option(c('--prefix', '-p'), help = "the prefix to the output files", default = ""),
    make_option(c('--min_samples', '-s'), help = "the minimum proportion or minimum number of samples", default = 2),
    make_option(c('--min_reads', '-m'), help = "the minimum read count per sample", default = 2),
    make_option(c('--counts', '-c'), help = "filter on raw counts", action="store_true", default = FALSE),
    make_option(c('--fpkm', '-t'), help = "filter on FPKM", action="store_true", default = FALSE),
    make_option(c('--remove_monoexons', '-e'), help = "whether to remove mono-exonic transcripts", action="store_true", default = FALSE)
)

#df <- vroom::vroom("all_samples/cupcake/all_samples.demux_fl_count.csv")
#gff_file <- "all_samples/cupcake/all_samples.cupcake.collapsed.gff"

#input <- "all_samples/cupcake/all_samples"
#output <- "all_samples/cupcake_filtered/all_samples"

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

in_folder <- opt$inFolder
gff_file <- opt$gff
output <- opt$prefix
min_samples <- as.numeric(opt$min_samples)
min_reads <- as.numeric(opt$min_reads)
run_code <- opt$runCode
# mode - filter either by FPKM or coverage
mode <- "fpkm"

remove_monoexons <- opt$remove_monoexons

stopifnot( min_samples < ncol(df) )
stopifnot(min_reads > 0)

# make sure input and output are valid 
#stopifnot(dir.exists(basename(counts_file) ) )
#stopifnot(dir.exists(basename(output) ) )

# parse arguments
#counts_file <- paste0(input, ".demux_fl_count.csv" )
#gff_file <- paste0(input, ".cupcake.collapsed.gff")

#stopifnot(file.exists(counts_file) )
stopifnot(file.exists(gff_file) )

counts_out <- paste0(output, "_filter_fpkm.csv")
gff_out <- paste0(output, "_filter.gtf")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))

count_files <- list.files(in_folder, pattern = "t_data.ctab", full.names = TRUE, recursive = TRUE)
count_files <- count_files[file.size(count_files) > 0]
count_files <- count_files[ grepl(run_code, count_files) ]
message(" * reading in transcript counts from ", length(count_files), " files")
count_res <- map(count_files, read_tsv)
names(count_res) <- gsub("sample_", "", basename(dirname(count_files)))
save.image("filter_stringtie_debug.RData")

df <- map2(count_res, names(count_res), ~{
        x <- data.frame(.x$FPKM, row.names = .x$t_name)
        names(x) <- .y
        return(x) }) %>% purrr::reduce(cbind)

# remove monoexons
# get rid of mono-exons if user requests
if( remove_monoexons == TRUE){
    anno <- count_res[[1]]
    multi <- filter(anno, num_exons > 1)
    message( " * removing ", nrow(anno) - nrow(multi), " monoexons")
    df <- df[multi$t_name,]
}

# keep all transcripts found at least once
# keep all annotated transcripts
# only keep novel if seen at least twice
clean <- df[ 
    ( grepl("ENST", row.names(df) ) & rowSums(df > 0) >= min_samples ) |
    ( grepl("MSTRG",row.names(df) ) & rowSums(df > 0) >= min_samples ) 
,] 

message( paste0( " * ", nrow(clean) , " out of ", nrow(df), " transcripts kept!" ) )


# filter GFF - using default gffread cleaned GFF3 files now
message(" * importing GFF")
gff_file_ext <- tools::file_ext(gff_file)
gff <- import(gff_file, format = gff_file_ext)


#save.image("debug.RData")

clean_gff <- gff[gff$transcript_id %in% row.names(clean)]
stopifnot(length(clean_gff) > 0 )

  

# remove transcripts without assigned strand - found in small number in Bambu GTF
message( paste0(" * removing ", length( clean_gff[strand(clean_gff) == "*" ] ), " transcripts with ambiguous stranding" ) )
clean_gff <- clean_gff[ strand(clean_gff) != "*" ]

clean <- clean[ unique(clean_gff$transcript_id), ]

message(" * writing filtered counts and TPM")

clean <- rownames_to_column(clean, var = "id")
readr::write_csv(clean, path = counts_out)

message(" * writing filtered GFF")
export(clean_gff, con = gff_out, format = "GTF")

