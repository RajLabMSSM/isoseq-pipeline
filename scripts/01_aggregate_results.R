#!/usr/bin/env Rscript
# 01_aggregate_results.R
#
# Aggregates longcallR ASE and ASJ output files across all samples into
# long-format tables for QC and consistency checks.
#
# Usage:
#   Rscript 01_aggregate_results.R \
#     --results-dir results \
#     --metadata metadata.tsv \
#     --out-dir aggregated
#
# Metadata file should have columns: sample, condition

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

# ── Arguments ────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--results-dir", type = "character", help = "Root results directory"),
  make_option("--metadata",    type = "character", help = "TSV with sample and condition columns"),
  make_option("--out-dir",     type = "character", default = "aggregated",
              help = "Output directory [default: aggregated]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`results-dir`) || is.null(opt$metadata)) {
  stop("--results-dir and --metadata are required")
}

dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

# ── Load metadata ─────────────────────────────────────────────────────────────

metadata <- read_tsv(opt$metadata, show_col_types = FALSE)
stopifnot(all(c("sample", "condition") %in% colnames(metadata)))
message("Loaded metadata for ", nrow(metadata), " samples")

# ── Helper: discover files ────────────────────────────────────────────────────

discover_files <- function(results_dir, samples, suffix) {
  map_dfr(samples, function(s) {
    path <- file.path(results_dir, s, "phased", paste0(s, suffix))
    if (!file.exists(path)) {
      warning("File not found, skipping: ", path)
      return(NULL)
    }
    tibble(sample = s, path = path)
  })
}

# ── Load ASE files ────────────────────────────────────────────────────────────

message("\nLoading ASE files...")

ase_files <- discover_files(opt$`results-dir`, metadata$sample, ".ase.tsv")
message("  Found ", nrow(ase_files), "/", nrow(metadata), " ASE files")

read_ase <- function(sample, path) {
  read_tsv(path, show_col_types = FALSE,
           col_names = c("Gene_name", "Chr", "PS", "H1", "H2", "P_value", "logFC"),
           comment = "#") %>%
    mutate(sample = sample,
           total  = H1 + H2,
           ratio  = H1 / total)
}

ase_long <- pmap_dfr(ase_files, read_ase) %>%
  left_join(metadata, by = "sample") %>%
  select(Gene_name, Chr, PS, sample, condition, H1, H2, total, ratio, logFC, P_value)

message("  Loaded ", nrow(ase_long), " ASE records across ",
        n_distinct(ase_long$sample), " samples and ",
        n_distinct(ase_long$Gene_name), " genes")

# ── Load ASJ files ────────────────────────────────────────────────────────────

message("\nLoading ASJ files...")

asj_files <- discover_files(opt$`results-dir`, metadata$sample, ".asj.tsv")
message("  Found ", nrow(asj_files), "/", nrow(metadata), " ASJ files")

read_asj <- function(sample, path) {
  read_tsv(path, show_col_types = FALSE,
           col_names = c("Junction", "Strand", "Junction_set", "Phase_set",
                         "Hap1_absent", "Hap1_present", "Hap2_absent", "Hap2_present",
                         "P_value", "SOR", "Novel", "GT_AG", "Gene_name"),
           comment = "#") %>%
    mutate(
      sample     = sample,
      feature_id = paste0(Gene_name, ":", Junction),
      H1_total   = Hap1_absent + Hap1_present,
      H2_total   = Hap2_absent + Hap2_present,
      H1_ratio   = Hap1_present / H1_total,
      H2_ratio   = Hap2_present / H2_total,
      logFC      = log2((Hap1_present + 0.5) / (H1_total + 0.5)) -
                   log2((Hap2_present + 0.5) / (H2_total + 0.5))
    )
}

asj_long <- pmap_dfr(asj_files, read_asj) %>%
  left_join(metadata, by = "sample") %>%
  select(feature_id, Gene_name, Junction, Strand, Junction_set, Phase_set,
         sample, condition, Novel, GT_AG,
         Hap1_present, Hap1_absent, H1_total, H1_ratio,
         Hap2_present, Hap2_absent, H2_total, H2_ratio,
         logFC, P_value, SOR)

message("  Loaded ", nrow(asj_long), " ASJ records across ",
        n_distinct(asj_long$sample), " samples and ",
        n_distinct(asj_long$feature_id), " unique junctions")

# ── Write outputs ─────────────────────────────────────────────────────────────

ase_out <- file.path(opt$`out-dir`, "ase_aggregated.tsv")
asj_out <- file.path(opt$`out-dir`, "asj_aggregated.tsv")

write_tsv(ase_long, ase_out)
write_tsv(asj_long, asj_out)

message("\nWritten:")
message("  ", ase_out)
message("  ", asj_out)

# ── Basic consistency summary ─────────────────────────────────────────────────

message("\n── ASE consistency summary ──────────────────────────────────────────")
ase_long %>%
  group_by(condition) %>%
  summarise(
    n_samples       = n_distinct(sample),
    n_genes         = n_distinct(Gene_name),
    median_total    = median(total),
    median_H1_ratio = median(ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

message("\n── ASJ consistency summary ──────────────────────────────────────────")
asj_long %>%
  group_by(condition) %>%
  summarise(
    n_samples        = n_distinct(sample),
    n_junctions      = n_distinct(feature_id),
    pct_novel        = mean(Novel) * 100,
    median_H1_ratio  = median(H1_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()
