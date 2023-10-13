## Add genes to GTF
## Jack Humphrey
## 2021

# gene entries needed for collapsing GTFs for TensorQTL, probably other stuff too
library(rtracklayer)


in_gtf <- "/sc/arion/projects/als-omics/microglia_isoseq/isoseq-pipeline/results/bambu/current_gtf/raj_roussos_miss_sqanti.cds.gtf.sorted.gtf"
out_gtf <- "/sc/arion/projects/als-omics/microglia_isoseq/isoseq-pipeline/results/bambu/current_gtf/raj_roussos_miss_sqanti.cds.gtf.sorted.gene.gtf"

gtf <- import(in_gtf)

# testing
#gtf <- gtf[1:1000]

gtf$transcript_name <- gtf$transcript_id
gtf$transcript_type <- "isoseq"

# split into genes
gene_split <- split(gtf, gtf$gene_id)

names(gene_split) <- NULL

# for a GRanges of transcripts for a gene, construct gene entry
write_gene_record <- function(x){
    start <- min( start(x) )
    end <- max(end(x) )
    gene_id <- unique(x$gene_id)
    chr <- unique(seqnames(x) )
    strand <- unique(strand(x) )

    gene <- GenomicRanges::GRanges(
        seqnames = chr, 
        ranges = IRanges::IRanges(start, end),
        strand = strand,
        source = unique(x$source),
        type = "gene",
        score = unique(x$score),
        phase = unique(x$phase), 
        gene_id = gene_id,
        gene_name = gene_id,
        gene_type = "isoseq"
    )
    record <- c(gene, x)
    return(record)
}

# endoapply - applies the function to each GRList element 
# and keeps the output in the same class
# lapply returns a list, not a GRList
genes <- unlist(endoapply(gene_split, write_gene_record) )  

genes$gene_name <- genes$gene_id
genes$gene_type <- "isoseq"
genes$level <- 0

export(genes, out_gtf, format = "GTF")

