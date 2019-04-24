# extract length of polyA tails

# when reads are mapped by minimap2, the polyA tails are softclipped.
# this information is put into the CIGAR string of the BAM file
# for a read mapping sense to the genome, the polyA comes at the end of the read
# for a read mapping antisense, the polyA comes at the start
# Softclipping is represented by "nS" in the CIGAR - n is the number of bases that were softclipped - the polyA length
# FLAGs -F 272 = 

library(dplyr)
library(ggplot2)
library(ggthemes) # for tufte box plots

## FUNCTIONS
parseBam <- function(bamFile, strand = "sense"){
  if(strand == "sense"){
    SAMflag <- "-F 16 "
  }else{
    SAMflag <- "-f 16 "
  }
    
  cmd <- paste0("samtools view ", SAMflag, bamFile, " |  bioawk -t -c sam \'{print $qname \":\" $cigar  \":\" $seq }\' ")
  bam <- system(cmd, intern = TRUE )
  split <- stringr::str_split_fixed(bam, "\\:", 3)
  read_id <- split[,1]
  CIGAR <- split[,2]
  read <- split[,3]
  
  return( data.frame(read_id, CIGAR, read))
}

# duplicate reads should have same softclipping, right?
#table(duplicated(read_ids))

sensePolyA <- function(bamFile){
  read_table <- parseBam(bamFile, strand = "sense")
  read_table$strand <- "+"
  # split off reads with no softclipping
  no_softclip <- read_table[ !grepl("S$|^[0-9]+S", read_table$CIGAR), ]
  
  mapped <- dplyr::filter(read_table, CIGAR != "*", read != "*")

  cigars_split <- stringr::str_split(mapped$CIGAR, pattern = "[MIDNSHPX=]")
  
  mapped$n_softclip <- 
    purrr::map_dbl(cigars_split, ~{
      len <- length(.x)
      as.numeric(.x[len-1])
    })

  mapped$softclip_seq <- 
    purrr::map2_chr(.x = mapped$read, .y = mapped$n_softclip, ~{
      len <- stringr::str_length(.x)
      softclip <- .y
      seq <- stringr::str_sub(string = .x, start = len - (softclip + 1), end = len)
      return(seq)
    })
  
  mapped$percent <- 
    purrr::map_dbl(.x = mapped$softclip_seq, ~{
      len <- stringr::str_length(.x)
      count <- stringr::str_count(string = .x,  pattern = "A")
      prop <- count / len
      return(prop)
    })
  
   mapped$polyA_pass <-
     ifelse(mapped$n_softclip >=10 & mapped$percent > 0.8, TRUE, FALSE)
  
   return(mapped)
}


antisensePolyA <- function(bamFile){
  read_table <- parseBam(bamFile, strand = "antisense")
  read_table$strand <- "-"
  
  no_softclip <- read_table[ !grepl("S$|^[0-9]+S", read_table$CIGAR), ]
  
  mapped <- dplyr::filter(read_table, CIGAR != "*", read != "*")
  
  cigars_split <- stringr::str_split(mapped$CIGAR, pattern = "[MIDNSHPX=]")
  
  mapped$n_softclip <- 
    purrr::map_dbl(cigars_split, ~{
      as.numeric(.x[1])
    })
  
  mapped$softclip_seq <- 
    purrr::map2_chr(.x = mapped$read, .y = mapped$n_softclip, ~{
      softclip <- .y
      seq <- stringr::str_sub(string = .x, start = 1, end = .y)
      return(seq)
    })
  
  mapped$percent <- 
    purrr::map_dbl(.x = mapped$softclip_seq, ~{
      len <- stringr::str_length(.x)
      count <- stringr::str_count(string = .x,  pattern = "T")
      prop <- count / len
      return(prop)
    })
  
  mapped$polyA_pass <-
    ifelse(mapped$n_softclip >=10 & mapped$percent > 0.8, TRUE, FALSE)
  
  
  return(mapped)
}

# bamFile <- "bams/sorted/MNP_Iso_BS6_Atailprep_hq_transcripts_sorted.bam"
# sense <- sensePolyA(bamFile)
# antisense <- antisensePolyA(bamFile)
# all <- rbind(sense,antisense)
# 

bamFiles <- list.files(path = "bams/sorted", full.names = TRUE, pattern = "*sorted.bam$")

all_polya <- 
  purrr::map( bamFiles, ~{
  bamFile <- .x
  sense <- sensePolyA(bamFile)
  antisense <- antisensePolyA(bamFile)
  all <- rbind(sense,antisense)
  all$bam <- bamFile
  return(all)
} ) %>% 
  purrr::reduce(rbind) %>%
  filter(polyA_pass == TRUE)

if( ! dir.exists("polya/")) {
  dir.create("polya/")
}
readr::write_tsv(all_polya, path = "polya/all_polyas.tsv" )

