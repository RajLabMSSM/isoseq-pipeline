library(rtracklayer)
library(Biostrings)
library(optparse)

#options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE)
options(BIOCONDUCTOR_CONFIG_FILE = "/sc/arion/projects/ad-omics/data/software/bambu/config.yaml")

option_list <- list(
    make_option(c('--anno', '-a' ), help='the annotation RData file produced by bambu', default = ""),
    make_option(c('--prefix', '-p'), help = "the prefix to the output files", default = ""),
    make_option(c('--fasta', '-f'), help = "the genome build", default = ""),
    make_option(c('--cores', '-c'), help = "the number of cores", default = 1),
    # NDR = Novel Discovery Rate threshold (0-1). Default NULL lets bambu auto-recommend
    # based on annotation completeness (current behaviour). HIGHER NDR (towards 1.0) =
    # MORE sensitive = more novel transcripts kept -> better recall of rare isoforms for
    # Join & Call discovery (Chen et al. 2023). Set e.g. --NDR 0.5 when you want to tune.
    make_option(c('--NDR', '-n'), help = "bambu NDR threshold 0-1 (default NULL = auto)", default = NULL)
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
# NDR: NULL -> bambu auto-recommends; otherwise use the supplied numeric threshold
ndr <- if (is.null(opt$NDR)) NULL else as.numeric(opt$NDR)

message(" * genome: ", fasta )

message(" * processing the following BAM files:" )

message(bam_files)

message(" * running bambu with ", cores, " cores" )

load(anno)

#suppressPackageStartupMessages(
library(bambu)#)

# run bambu
# Join & Call: all bam_files are passed together so discovery is joint across the
# group's replicates. NDR defaults to NULL (bambu auto-recommends); pass --NDR to tune.
res <- bambu(
    reads = bam_files,
    annotations = anno,
    genome = fasta,
    #stranded = FALSE,
    NDR = ndr,
    ncore = cores
    #rcOutDir = dirname(prefix)
    #lowMemory = TRUE
    )

# save GTF and counts
writeBambuOutput(se = res, path = dirname(prefix), prefix = paste0(basename(prefix), "_" ) )

# save bambu object

save( res, file = paste0(prefix, "_bambu.RData") )

