library(optparse)

option_list <- list(
    make_option(c('--anno', '-a' ), help='the annotation RData file produced by bambu', default = ""),
    make_option(c('--prefix', '-p'), help = "the prefix to the output files", default = ""),
    make_option(c('--fasta', '-f'), help = "the genome build", default = ""),
    make_option(c('--cores', '-c'), help = "the number of cores", default = 1)
)

option.parser <- OptionParser(option_list=option_list)
# allow for unlimited positional arguments
arguments <- parse_args(option.parser,positional_arguments = TRUE)

bam_files <- arguments$args

opt <- arguments$opt

anno <- opt$anno
prefix <- opt$prefix
fasta <- opt$fasta
cores <- as.numeric(opt$cores)

message(" * genome: ", fasta )

message(" * processing the following BAM files:" )

print(bam_files)

message(" * running bambu with", cores, "cores" )

load(anno)

suppressPackageStartupMessages(library(bambu))

# run bambu
res <- bambu(
    reads = bam_files,
    annotations = anno,
    genome = fasta,
    stranded = FALSE,
    ncore = cores,
    rcOutDir = dirname(prefix)
    )

# save GTF and counts
writeBambuOutput(se = res, path = dirname(prefix), prefix = paste0(basename(prefix), "_" ) )

# save bambu object

save( res, file = paste0(prefix, "_bambu.RData") )

