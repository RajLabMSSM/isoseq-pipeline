library(optparse)

option_list <- list(
    make_option(c('--sqanti', '-s' ), help = "the SQANTI classifications", default = ""),
    make_option(c('--fasta', '-f' ), help = "the FASTA sequences", default = ""),
    make_option(c('--gff', '-g'), help = "the CDS GFF file from SQANTI", default = ""),
    make_option(c('--input', '-i'), help = "the prefix to the input files", default = ""),
    make_option(c('--output', '-o'), help = "the prefix to the output files", default = ""),
    make_option(c('--counts', '-c'), help = "the counts file", default = "")

)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

output <- opt$output
input <- opt$input
sqanti_file <- opt$sqanti
fasta_file <- opt$fasta
counts_file <- opt$counts
#counts_file <- paste0(input, "_miss_counts.csv")
#tpm_file <- paste0(input, "_miss_tpm.csv")
gff_file <- opt$gff

stopifnot(file.exists(counts_file) )
stopifnot(file.exists(gff_file) )
stopifnot(file.exists(sqanti_file))


sqanti_out <- paste0(output, "_filter_sqanti_classification.tsv")
counts_out <- paste0(output, "_filter_sqanti_counts.csv")
#tpm_out <- paste0(output, "_filter_sqanti_tpm.csv")
gff_out <- paste0(output, "_filter_sqanti.cds.gtf")
fasta_out <- paste0(output, "_filter_sqanti.fasta")

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))


# read in SQANTI classifications
pre <-  read_tsv(sqanti_file)

pre$structural_category <- 
  gsub("-|_", " ", pre$structural_category)

# run of A
calculate_polyA_run <- function(seqs){
  # for vector of sequences
  a_runs <- strsplit(seqs, "") %>% 
    purrr::map_dbl(~rle(.x)$lengths[1])
  a_runs[ !grepl("^A", seqs)] <- 0
  return(a_runs)
}

pre$N_A_downstream_TTS <- calculate_polyA_run(pre$seq_A_downstream_TTS)

# Roy-Chanfreau score
calculate_rc_score <- function(seqs){
  # for vector of sequences
  a_start <- stringr::str_count(stringr::str_sub(seqs,start = 1, end = 6), "A")
  a_total <- stringr::str_count(stringr::str_sub(seqs,start = 1, end = 19), "A")

  a_score <- a_start + a_total
  return(a_score)
}

pre$rc_score <- calculate_rc_score(pre$seq_A_downstream_TTS)


# read in counts, TPM and GTF
counts <- read_csv(counts_file)
#tpm <- read_csv(tpm_file)
gff <- import(gff_file, format = "GFF")

# filter transcripts using SQANTI settings
# now done using diff_to_gene_TTS
#fsm_polya <- 
#  filter(pre, structural_category == "full splice match" ) %>% dplyr::select(associated_gene, seq_A_downstream_TTS) %>%
#  inner_join(pre)

#pre$fsm_TTS_match <- pre$isoform %in% fsm_polya$isoform

post <- pre %>% 
  mutate( filter_pass = case_when(
    structural_category == "full splice match" ~ TRUE,
    structural_category != "full splice match" & 
      RTS_stage == FALSE & # no RT-switching junction
      min_sample_cov >= 5 & min_cov >= 50 & # each junction supported by at least 5 reads in at least 10 short read samples.
      ( 
        diff_to_gene_TTS == 0 | # has an annotated TTS - in GENCODE v30
        #fsm_TTS_match == TRUE | # has an annotated TTS
        (N_A_downstream_TTS < 6 & perc_A_downstream_TTS < 60 & rc_score <= 15) # or if not, TTS is low A content
      ) ~ TRUE,
    TRUE ~ FALSE
  ))


post <- filter(post, filter_pass == TRUE)
message( " * filtering transcripts based on SQANTI annotation" )
message( " * kept ", nrow(post), " transcripts!" )

reasons <- post %>%
  filter(filter_pass == FALSE) %>%
  transmute( isoform, type = structural_category, RTS = RTS_stage == TRUE,
          coverage = min_sample_cov < 5 | min_cov < 50 | is.na(min_sample_cov) | is.na(min_cov),
          N_A = N_A_downstream_TTS >= 6,
          perc_A = perc_A_downstream_TTS >= 60,
          rc = rc_score > 15)


# filter output files
counts <- counts[ counts$id %in% post$isoform,]
#tpm <- tpm[ tpm$id %in% post$isoform,] 
gff <- gff[ gff$transcript_id %in% post$isoform ]

# read in FASTA and filter
message(" * reading FASTA")
fasta <- readDNAStringSet(fasta_file)

fasta <- fasta[ post$isoform ]


message(" * writing to", sqanti_out)
write_tsv(post, file = sqanti_out)
message( " * writing to ", counts_out )
write_csv(counts, file = counts_out)
#message(" * writing to ", tpm_out)
#write_csv(tpm, file = tpm_out)
message(" * writing to ", gff_out)
export(gff, con = gff_out, format = "GTF")
message(" * writing to ", fasta_out )
writeXStringSet(fasta, fasta_out)
