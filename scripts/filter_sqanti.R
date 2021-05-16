library(optparse)

option_list <- list(
    make_option(c('--sqanti', '-s' ), help = "the SQANTI classifications", default = ""),
    make_option(c('--fasta', '-f' ), help = "the FASTA sequences", default = ""),
    make_option(c('--input', '-i'), help = "the prefix to the input files", default = ""),
    make_option(c('--output', '-o'), help = "the prefix to the output files", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

output <- opt$output
input <- opt$input
sqanti_file <- opt$sqanti
fasta_file <- opt$fasta

counts_file <- paste0(input, "_miss_counts.csv")
tpm_file <- paste0(input, "_miss_tpm.csv")
gtf_file <- paste0(input, "_miss.gtf")

stopifnot(file.exists(counts_file) )
stopifnot(file.exists(gtf_file) )
stopifnot(file.exists(sqanti_file))


sqanti_out <- paste0(output, "_miss_sqanti_classification.tsv")
counts_out <- paste0(output, "_miss_sqanti_counts.csv")
tpm_out <- paste0(output, "_miss_sqanti_tpm.csv")
gtf_out <- paste0(output, "_miss_sqanti.gtf")
fasta_out <- paste0(output, "_miss_sqanti.fasta")

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
tpm <- read_csv(tpm_file)
gtf <- import(gtf_file, format = "GTF")

# filter transcripts using SQANTI settings
fsm_polya <- 
  filter(pre, structural_category == "full splice match" ) %>% dplyr::select(associated_gene, seq_A_downstream_TTS) %>%
  inner_join(pre)

pre$fsm_TTS_match <- pre$isoform %in% fsm_polya$isoform

post <- pre %>% 
  mutate( filter_pass = case_when(
    structural_category == "full splice match" ~ TRUE,
    structural_category != "full splice match" & 
      RTS_stage == FALSE &
      ( 
        (N_A_downstream_TTS < 6 & perc_A_downstream_TTS < 60 & rc_score <= 15) |
        fsm_TTS_match == TRUE
      ) ~ TRUE,
    TRUE ~ FALSE
  ))


post <- filter(post, filter_pass == TRUE)

# filter output files
counts <- counts[ counts$id %in% post$isoform,]
tpm <- tpm[ tpm$id %in% post$isoform,] 
gtf <- gtf[ gtf$transcript_id %in% post$isoform ]

# read in FASTA and filter
message(" * reading FASTA")
fasta <- readDNAStringSet(fasta_file)

fasta <- fasta[ post$isoform ]




message(" * writing to", sqanti_out)
write_tsv(post, file = sqanti_out)
message( " * writing to ", counts_out )
write_csv(counts, file = counts_out)
message(" * writing to ", tpm_out)
write_csv(tpm, file = tpm_out)
message(" * writing to ", gtf_out)
export(gtf, con = gtf_out, format = "GTF")
message(" * writing to ", fasta_out )
writeXStringSet(fasta, fasta_out)
