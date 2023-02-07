# combine novel transcripts with a set of annotated transcripts
# Jack Humphrey

# Inputs
# - annotated transcripts in GTF and FASTA format
# - set of long-read derived transcripts in GTF and FASTA format

# Options
# - whether to include novel genes or not

# Outputs
# - combined annotated + novel transcripts in GTF and FASTA format

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c('--annoGTF', '-a' ), help='annotated transcripts in GTF format', default = ""),
    make_option(c('--annoFASTA', '-b' ), help="annotated transcript in FASTA format", default = ""),
    make_option(c('--novelGTF', '-c' ), help = "long read transcripts in GTF format", default = ""),
    make_option(c('--novelFASTA', '-d'), help = "the longread transcripts in FASTA format", default = ""),
    make_option(c('--out', '-o'), help = "the stem of the output", default = "test"),
    make_option(c('--novelGenes', '-n'), help = "if used, remove novel genes from combined reference", action="store_true", default = FALSE)
)

read_gtf <- function(gff_file){
    gff_file_ext <- tools::file_ext(gff_file)
    gff <- rtracklayer::import(gff_file, format = gff_file_ext)
    return(gff)
}

read_fasta <- function(file){
    fasta <- Biostrings::readDNAStringSet(filepath = file)
    return(fasta)
}

write_gtf <- function(gff, gff_file){
   rtracklayer::export(gff, con = gff_file, format = "GTF") 
}

write_fasta <- function(fasta, fasta_file){
    Biostrings::writeXStringSet(fasta, format = "FASTA", filepath = fasta_file)
}

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

anno_gtf <- opt$annoGTF
anno_fasta <- opt$annoFASTA
novel_gtf <- opt$novelGTF
novel_fasta <- opt$novelFASTA
out <- opt$out
novel_genes <- opt$novelGenes

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Biostrings))

# read in GTF and FASTA files
message(" * reading in annotated GTF from ", anno_gtf)
anno_g <- read_gtf(anno_gtf)
message(" * reading in annotated FASTA from ", anno_fasta)
anno_f <- read_fasta(anno_fasta)
message(" * reading in novel GTF from ", novel_gtf)
novel_g <- read_gtf(novel_gtf)
message(" * reading in novel FASTA from ", novel_fasta)
novel_f <- read_fasta(novel_fasta)

message(" * Annotated reference contains ", length(anno_f), " transcripts" )
message(" * Novel reference contains ", length(novel_f), " transcripts")

# extract novel transcripts from novel GTF and add to annotated GTF
# must be agnostic to how novel transcripts are named
novel_tx_id <- setdiff(novel_g$transcript_id, anno_g$transcript_id)
novel_gene_id <- setdiff(novel_g$gene_id, anno_g$gene_id)

tx2gene <- data.frame(gene_id = novel_g$gene_id, transcript_id = novel_g$transcript_id) %>% distinct()

if(!novel_genes){
    message(" * removing novel genes")
    tx_id_to_add <- filter(tx2gene, !gene_id %in% novel_gene_id & transcript_id %in% novel_tx_id) %>% pull(transcript_id)
}else{
    message(" * keeping novel genes")
    tx_id_to_add <- filter(tx2gene, transcript_id %in% novel_tx_id) %>% pull(transcript_id)
}
message(" * adding ", length(tx_id_to_add), " novel transcripts to the annotation")

to_add_g <- novel_g[ novel_g$transcript_id %in% tx_id_to_add ]

to_add_f <- novel_f[ tx_id_to_add ]

# combine
combined_g <- c(anno_g, to_add_g)
combined_f <- c(anno_f, to_add_f)

# write out
out_g <- paste0(out, ".gtf")
out_f <- paste0(out, ".fa")

message(" * writing GTF to ", out_g)
write_gtf(combined_g, out_g)

message(" * writing FASTA to ", out_f)
write_fasta(combined_f, out_f)








