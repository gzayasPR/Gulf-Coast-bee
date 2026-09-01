#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

rm(list = ls())

# ----------------- Args -----------------
# 1) out.dir    : output directory (also working dir)
# 2) het.path   : combined .het file => ${out_dir}/${name}.het
# 3) meta.path  : metadata CSV (first col is BioSample)
args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript Hetero.V3.R <out.dir> <het.path> <meta.path>\n")
}
out.dir   <- args[1]
het.path  <- args[2]
meta.path <- args[3]

setwd(out.dir)

# ----------------- Load data -----------------
# vcftools --het header is: INDV  O(HOM)  E(HOM)  N_SITES  F
het_data <- read.table(het.path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
# Standardize sample ID column name to 'BioSample'
# (vcftools column is 'INDV'; if you've already changed it, this keeps it consistent)
if ("INDV" %in% names(het_data)) names(het_data)[names(het_data) == "INDV"] <- "BioSample"
if ("BioSample" %in% names(het_data) == FALSE) stop("Could not find 'INDV' or 'BioSample' in het file.")

# Compute observed/expected heterozygosity (%)
# Columns from vcftools are typically named with punctuation: 'O(HOM)' -> 'O.HOM.' when read with check.names=TRUE
# We used check.names=FALSE, so access with backticks to be safe.
required_cols <- c("N_SITES", "O(HOM)", "E(HOM)")
if (!all(required_cols %in% names(het_data))) {
  stop("Expected columns not present in het file: ", paste(setdiff(required_cols, names(het_data)), collapse = ", "))
}
het_data <- het_data %>%
  mutate(
    Observed_Heterozygosity = ((N_SITES - `O(HOM)`) / N_SITES) * 100,
    Expected_Heterozygosity = ((N_SITES - `E(HOM)`) / N_SITES) * 100
  )

# Inbreeding coefficient column 'F'
if (!"F" %in% names(het_data)) stop("Column 'F' not found in het file.")

# ----------------- Metadata -----------------
metadata <- read.csv(meta.path, stringsAsFactors = FALSE, check.names = FALSE)
# Standardize metadata sample ID column to 'BioSample'
if (!"BioSample" %in% names(metadata)) names(metadata)[1] <- "BioSample"

# Optional: recode location labels
if ("Location" %in% names(metadata)) {
  metadata <- metadata %>%
    mutate(Location = dplyr::recode(Location,
                                    "USA-AL_Baldwin_co" = "AL",
                                    "USA-FL_Escambia_co" = "FL"))
}

# ----------------- Merge & filter -----------------
merged_data <- het_data %>%
  inner_join(metadata, by = "BioSample")

# Keep only females if 'Sex' is present
if ("Sex" %in% names(merged_data)) {
  merged_data <- merged_data %>% filter(Sex == "female")
}

# ----------------- Plots -----------------
# Observed heterozygosity by location
if ("Location" %in% names(merged_data)) {
  p_obs <- ggplot(merged_data, aes(x = Location, y = Observed_Heterozygosity, fill = Location)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.3, alpha = 0.7) +
    labs(title = "Observed Heterozygosity by Location",
         x = "Location", y = "Observed Heterozygosity (%)") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          legend.position = "none") +
    scale_fill_manual(values = c("AL" = "#4CAF50", "FL" = "#2196F3")) +
    theme(panel.background = element_rect(fill = "white", color = "white"))
  ggsave(file.path(out.dir, "Observed_Heterozygosity_byLocation.png"), p_obs, width = 11, height = 6, dpi = 300)

  # Inbreeding coefficient (F) by location
  p_F <- ggplot(merged_data, aes(x = Location, y = F, fill = Location)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.3, alpha = 0.7) +
    labs(title = "Inbreeding Coefficient (F) by Location",
         x = "Location", y = "F (Wright's inbreeding coefficient)") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          legend.position = "none") +
    scale_fill_manual(values = c("AL" = "#4CAF50", "FL" = "#2196F3")) +
    theme(panel.background = element_rect(fill = "white", color = "white"))
  ggsave(file.path(out.dir, "InbreedingCoefficient_byLocation.png"), p_F, width = 11, height = 6, dpi = 300)
} else {
  message("No 'Location' column found in metadata; saving overall distributions instead.")
  p_obs <- ggplot(merged_data, aes(x = Observed_Heterozygosity)) +
    geom_histogram(bins = 30) +
    theme_bw(base_size = 15) +
    labs(title = "Observed Heterozygosity", x = "Observed Heterozygosity (%)", y = "Count")
  ggsave(file.path(out.dir, "Observed_Heterozygosity_hist.png"), p_obs, width = 8, height = 5, dpi = 300)

  p_F <- ggplot(merged_data, aes(x = F)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
    theme_bw(base_size = 15) +
    labs(title = "Inbreeding Coefficient (F)", x = "F", y = "Count")
  ggsave(file.path(out.dir, "InbreedingCoefficient_hist.png"), p_F, width = 8, height = 5, dpi = 300)
}

# ----------------- Save table -----------------
write.csv(merged_data, file.path(out.dir, "Heterozygosity_merged.csv"), row.names = FALSE)
cat("Saved:\n - Observed_Heterozygosity_byLocation.png\n - InbreedingCoefficient_byLocation.png (if Location present)\n - Heterozygosity_merged.csv\n")
