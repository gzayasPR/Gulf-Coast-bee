rm(list = ls())
# Set the working directory
proj.dir <- "C:/Users/gzayas97/Documents/USDA Internship/"
setwd(proj.dir)
data.dir <- paste(proj.dir, "/data/", sep="")
results.dir <- paste(proj.dir, "/results/", sep="")

# Install necessary packages if you haven't already
if (!requireNamespace("vcfR", quietly = TRUE)) {
  install.packages("vcfR")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# Load the packages
library(vcfR)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(openxlsx)


# Read the VCF file
vcf.path <- paste(data.dir, "/Hesperapis_oraria.vcf", sep="")
vcf <- read.vcfR(vcf.path)

# Extract the genotype (GT) and depth (DP) information
meta.path <- paste(data.dir, "/Samples.Metadata.csv", sep="")
meta.data <- read.csv(meta.path)
females.id <- meta.data$BioSample[meta.data$Sex == "female"]

# Get the column indices of the individuals
individual_indices <- which(colnames(vcf@gt) %in% females.id)

# Subset the VCF object
females.vcf <- vcf[, c(1, individual_indices)]

# Define the depth thresholds
depth_thresholds <- c(3,5,7,9, 11, 5,9,13,15)
# Function to apply filters and perform analysis
analyze_vcf <- function(vcf, depth_threshold, min_called_fraction = 0.25) {
  # Filter by depth
  dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  gt_matrix <- extract.gt(vcf, element = "GT")
  
  # Mask genotypes with depth less than the threshold
  gt_matrix[dp_matrix < depth_threshold] <- NA
 # gt_matrix[dp_matrix > 35] <- NA
  
  # Remove variants with less than the minimum called fraction of individuals
  called_fraction <- rowMeans(is.na(gt_matrix))
  valid_variants <- called_fraction <= min_called_fraction
  gt_matrix <- gt_matrix[valid_variants, ]
  dp_matrix <- dp_matrix[valid_variants, ]
  
  # Calculate depth statistics
  avg_depth_per_individual <- colMeans(dp_matrix, na.rm = TRUE)
  median_depth_per_individual <- apply(dp_matrix, 2, median, na.rm = TRUE)
  na_counts <- colSums(is.na(dp_matrix))
  na_percentage <- na_counts / nrow(dp_matrix)
  
  # Calculate observed heterozygosity (Ho) per individual
  het_counts_per_individual <- apply(gt_matrix, 2, function(geno) {
    sum(geno %in% c("0/1", "1/0", "0|1", "1|0"), na.rm = TRUE)
  })
  hom_ref_counts_per_individual <- apply(gt_matrix, 2, function(geno) {
    sum(geno %in% c("0|0", "0/0"), na.rm = TRUE)
  })
  hom_alt_counts_per_individual <- apply(gt_matrix, 2, function(geno) {
    sum(geno %in% c("1/1", "1|1"), na.rm = TRUE)
  })
  n_alleles_per_individual <- colSums(!is.na(gt_matrix), na.rm = TRUE)
  ho_per_individual <- het_counts_per_individual / n_alleles_per_individual
  
  # Calculate expected heterozygosity (He) per individual
  allele_freqs_per_site <- apply(gt_matrix, 1, function(geno) {
    af <- table(factor(unlist(strsplit(geno, "[/|]")), levels = c("0", "1")))
    af / sum(af)
  })
  
  he_per_individual <- apply(gt_matrix, 2, function(ind_geno) {
    af_per_site <- allele_freqs_per_site[, !is.na(ind_geno)]
    he_per_site <- 2 * af_per_site["0", ] * af_per_site["1", ]
    mean(he_per_site, na.rm = TRUE)
  })
  
  # Number of markers retained after filtering
  num_markers_retained <- nrow(gt_matrix)
  
  # Summary of results per individual
  summary_results_per_individual <- data.frame(
    individual = colnames(dp_matrix),
    avg_depth = avg_depth_per_individual,
    median_depth = median_depth_per_individual,
    observed_heterozygosity = ho_per_individual,
    expected_heterozygosity = he_per_individual,
    hom_ref_counts = hom_ref_counts_per_individual,
    het_counts = het_counts_per_individual,
    hom_alt_counts = hom_alt_counts_per_individual,
    na_counts = na_counts,
    na_percentage = na_percentage,
    num_markers = num_markers_retained,
    min_depth = depth_threshold
  )
  
  return(summary_results_per_individual)
}

# List to store results
results_list <- list()

# Iterate over each depth threshold and store results
for (depth_threshold in depth_thresholds) {
  cat("Analyzing with minimum depth:", depth_threshold, "\n")
  results <- analyze_vcf(females.vcf, depth_threshold)
  results_list[[paste0("depth_", depth_threshold)]] <- results
}

# Combine results into a single data frame
combined_results <- bind_rows(results_list)

# Plot all results together
plot1 <- ggplot(combined_results, aes(x = avg_depth, y = observed_heterozygosity, 
                                      color = factor(min_depth))) +
  geom_point() +
  geom_smooth(aes(x = avg_depth, y = observed_heterozygosity),
              method = "gam",inherit.aes = F) +
  scale_x_continuous(
    breaks = seq(0, 30, by = 5) )+
  labs(title = "Observed Heterozygosity vs Average Read Depth",
       x = "Average Read Depth",
       y = "Observed Heterozygosity",
       color = "Min Depth")

print(plot1)

plot2 <- ggplot(combined_results, aes(x = median_depth, y = observed_heterozygosity, 
                        color = factor(min_depth))) +
  geom_point()+
  geom_smooth(aes(x = median_depth, y = observed_heterozygosity) ,
              method = "gam",inherit.aes = F) + 
  scale_x_continuous(
  breaks = seq(0, 30, by = 5) )+
  labs(title = "Observed Heterozygosity vs Median Read Depth",
       x = "Median Read Depth",
       y = "Observed Heterozygosity",
       color = "Min Depth")
print(plot2)
plot3 <- ggplot(combined_results, aes(x = avg_depth, y = expected_heterozygosity,
                                      color = factor(min_depth))) +
  geom_point()+
  geom_smooth(aes(x = median_depth, y = expected_heterozygosity) ,
              method = "gam",inherit.aes = F) + 
  labs(title = "Expected Heterozygosity vs Average Depth",
       x = "Average Read Depth",
       y = "Expected Heterozygosity",
       color = "Min Read Depth")
print(plot3)
plot4 <- ggplot(combined_results, aes(x = median_depth, y = expected_heterozygosity, color = factor(min_depth))) +
  geom_point() +
  geom_smooth(aes(x = median_depth, y = expected_heterozygosity) ,
              method = "gam",inherit.aes = F) + 
  labs(title = "Expected Heterozygosity vs Median Depth",
       x = "Median Read Depth",
       y = "Expected Heterozygosity",
       color = "Min Read Depth")
print(plot4)
# Arrange the plots in a 2x2 grid
combined_plot <- grid.arrange(plot1, plot2, ncol = 1, nrow = 2)


setwd(results.dir)
dir.create("DepthvsHet",showWarnings = F)
setwd("DepthvsHet")
ggsave("Ohet vs depth .png",combined_plot)

write.csv(combined_results,"test.csv",row.names = F,quote = F)
