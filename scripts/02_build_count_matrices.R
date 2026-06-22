#!/usr/bin/env Rscript
# 02_build_count_matrices.R
#
# Builds haplotype-split count matrices from aggregated ASE and ASJ tables:
#   - ASE: H1/H2 count matrix per gene for DEXSeq
#   - ASJ: present/absent count matrices per junction per haplotype for edgeR
#
# Usage:
#   Rscript 02_build_count_matrices.R \
#     --aggregated-dir aggregated \
#     --metadata metadata.tsv \
#     --out-dir matrices
#
# Expects outputs from 01_aggregate_results.R in --aggregated-dir

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--aggregated-dir", type = "character",
              help = "Directory containing ase_aggregated.tsv and asj_aggregated.tsv"),
  make_option("--metadata",       type = "character",
              help = "TSV with sample and condition columns"),
  make_option("--out-dir",        type = "character", default = "matrices",
              help = "Output directory [default: matrices]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`aggregated-dir`) || is.null(opt$metadata)) {
  stop("--aggregated-dir and --metadata are required")
}

dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

# ── Load metadata ─────────────────────────────────────────────────────────────

metadata <- read_tsv(opt$metadata, show_col_types = FALSE)
stopifnot(all(c("sample", "condition") %in% colnames(metadata)))

# Sample order: sort by condition then sample name for consistency
sample_order <- metadata %>% arrange(condition, sample) %>% pull(sample)

# ── Load aggregated tables ────────────────────────────────────────────────────

ase_long <- read_tsv(file.path(opt$`aggregated-dir`, "ase_aggregated.tsv"),
                     show_col_types = FALSE)
asj_long <- read_tsv(file.path(opt$`aggregated-dir`, "asj_aggregated.tsv"),
                     show_col_types = FALSE)

# ── ASE count matrix for DEXSeq ───────────────────────────────────────────────
#
# DEXSeq models haplotype imbalance by treating H1 and H2 as two "sub-samples"
# per donor. The matrix has one column per haplotype per sample:
#   <sample>_H1, <sample>_H2
# Rows are genes. The interaction term haplotype:condition tests for
# differential allelic imbalance between conditions.
#
# Accompanying coldata has columns:
#   sample_id  : e.g. CTRL1_H1
#   sample     : e.g. CTRL1
#   haplotype  : H1 or H2
#   condition  : scramble or TDP43_KD

message("Building ASE count matrix for DEXSeq...")

ase_h1 <- ase_long %>%
  select(Gene_name, sample, H1) %>%
  mutate(col = paste0(sample, "_H1")) %>%
  select(Gene_name, col, count = H1)

ase_h2 <- ase_long %>%
  select(Gene_name, sample, H2) %>%
  mutate(col = paste0(sample, "_H2")) %>%
  select(Gene_name, col, count = H2)

ase_matrix <- bind_rows(ase_h1, ase_h2) %>%
  pivot_wider(names_from = col, values_from = count) %>%
  column_to_rownames("Gene_name")

# Reorder columns: all H1 then H2, within each sorted by sample order
col_order <- c(paste0(sample_order, "_H1"), paste0(sample_order, "_H2"))
col_order <- col_order[col_order %in% colnames(ase_matrix)]
ase_matrix <- ase_matrix[, col_order]

# DEXSeq coldata
ase_coldata <- tibble(
  sample_id = col_order,
  sample    = sub("_H[12]$", "", col_order),
  haplotype = sub(".*_(H[12])$", "\\1", col_order)
) %>%
  left_join(metadata, by = "sample")

message("  ASE matrix: ", nrow(ase_matrix), " genes x ", ncol(ase_matrix), " columns")

# ── ASJ count matrices for edgeR ──────────────────────────────────────────────
#
# For junction-level allelic splicing, edgeR models the present:absent ratio
# per haplotype. We build four matrices (present/absent x H1/H2), all with
# the same rows (feature_id = Gene_name:Junction) and columns (samples).
#
# The edgeR model tests:
#   present_counts ~ offset(log(total)) + haplotype + condition + haplotype:condition
# where the interaction captures differential allelic splicing.

message("Building ASJ count matrices for edgeR...")

make_asj_matrix <- function(data, value_col, sample_order) {
  data %>%
    select(feature_id, sample, value = all_of(value_col)) %>%
    pivot_wider(names_from = sample, values_from = value) %>%
    column_to_rownames("feature_id") %>%
    .[, sample_order[sample_order %in% colnames(.)]]
}

asj_H1_present <- make_asj_matrix(asj_long, "Hap1_present", sample_order)
asj_H1_absent  <- make_asj_matrix(asj_long, "Hap1_absent",  sample_order)
asj_H2_present <- make_asj_matrix(asj_long, "Hap2_present", sample_order)
asj_H2_absent  <- make_asj_matrix(asj_long, "Hap2_absent",  sample_order)

# edgeR coldata — one row per sample
asj_coldata <- metadata %>%
  filter(sample %in% colnames(asj_H1_present)) %>%
  arrange(match(sample, sample_order))

message("  ASJ matrices: ", nrow(asj_H1_present), " junctions x ",
        ncol(asj_H1_present), " samples (4 matrices)")

# ── Junction annotation table ─────────────────────────────────────────────────
# Useful for downstream annotation of significant hits

junction_annot <- asj_long %>%
  distinct(feature_id, Gene_name, Junction, Strand, Junction_set, Novel, GT_AG)

# ── Write outputs ─────────────────────────────────────────────────────────────

write_tsv(as.data.frame(ase_matrix) %>% tibble::rownames_to_column("Gene_name"),
          file.path(opt$`out-dir`, "ase_counts_dexseq.tsv"))
write_tsv(ase_coldata,
          file.path(opt$`out-dir`, "ase_coldata_dexseq.tsv"))

write_tsv(as.data.frame(asj_H1_present) %>% tibble::rownames_to_column("feature_id"),
          file.path(opt$`out-dir`, "asj_H1_present.tsv"))
write_tsv(as.data.frame(asj_H1_absent) %>% tibble::rownames_to_column("feature_id"),
          file.path(opt$`out-dir`, "asj_H1_absent.tsv"))
write_tsv(as.data.frame(asj_H2_present) %>% tibble::rownames_to_column("feature_id"),
          file.path(opt$`out-dir`, "asj_H2_present.tsv"))
write_tsv(as.data.frame(asj_H2_absent) %>% tibble::rownames_to_column("feature_id"),
          file.path(opt$`out-dir`, "asj_H2_absent.tsv"))
write_tsv(asj_coldata,
          file.path(opt$`out-dir`, "asj_coldata_edger.tsv"))
write_tsv(junction_annot,
          file.path(opt$`out-dir`, "junction_annotation.tsv"))

message("\nWritten to ", opt$`out-dir`, ":")
message("  ase_counts_dexseq.tsv     — gene x (sample_H1/H2) count matrix")
message("  ase_coldata_dexseq.tsv    — DEXSeq sample metadata")
message("  asj_H1_present.tsv        — junction x sample H1 present counts")
message("  asj_H1_absent.tsv         — junction x sample H1 absent counts")
message("  asj_H2_present.tsv        — junction x sample H2 present counts")
message("  asj_H2_absent.tsv         — junction x sample H2 absent counts")
message("  asj_coldata_edger.tsv     — edgeR sample metadata")
message("  junction_annotation.tsv   — junction feature annotations")
