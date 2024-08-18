# Set the working directory
setwd("D:/beenome100/variants")

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

# Load the packages
library(vcfR)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Read the VCF file
vcf <- read.vcfR("fixed_missing.25.recode.vcf")

# Extract the genotype (GT) and depth (DP) information
meta.data <- read.csv("../Samples.Metadata.csv")
females.id <- meta.data$BioSample[meta.data$Sex == "female"]

# Get the column indices of the individuals
individual_indices <- which(colnames(vcf@gt) %in% females.id)

# Subset the VCF object
females.vcf <- vcf[, c(1, individual_indices)]

# Calculate depth (DP) per individual
dp_matrix <- extract.gt(females.vcf, element = "DP", as.numeric = TRUE)
avg_depth_per_individual <- colMeans(dp_matrix, na.rm = TRUE)
median_depth_per_individual <- apply(dp_matrix, 2, median, na.rm = TRUE)

# Calculate observed heterozygosity (Ho) per individual
gt_matrix <- extract.gt(females.vcf, element = "GT")
het_counts_per_individual <- apply(gt_matrix, 2, function(geno) {
  sum(geno %in% c("0/1", "1/0", "0|1", "1|0"), na.rm = TRUE)
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

# Summary of results per individual
summary_results_per_individual <- data.frame(
  individual = colnames(dp_matrix),
  avg_depth = avg_depth_per_individual,
  median_depth = median_depth_per_individual,
  observed_heterozygosity = ho_per_individual,
  expected_heterozygosity = he_per_individual
)

# Plot observed heterozygosity vs average depth
plot1 <- ggplot(summary_results_per_individual, aes(x = avg_depth, y = observed_heterozygosity)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Observed Heterozygosity vs Average Depth",
       x = "Average Depth",
       y = "Observed Heterozygosity")

# Plot observed heterozygosity vs median depth
plot2 <- ggplot(summary_results_per_individual, aes(x = median_depth, y = observed_heterozygosity)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Observed Heterozygosity vs Median Depth",
       x = "Median Depth",
       y = "Observed Heterozygosity")

# Plot expected heterozygosity vs average depth
plot3 <- ggplot(summary_results_per_individual, aes(x = avg_depth, y = expected_heterozygosity)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Expected Heterozygosity vs Average Depth",
       x = "Average Depth",
       y = "Expected Heterozygosity")

# Plot expected heterozygosity vs median depth
plot4 <- ggplot(summary_results_per_individual, aes(x = median_depth, y = expected_heterozygosity)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Expected Heterozygosity vs Median Depth",
       x = "Median Depth",
       y = "Expected Heterozygosity")

# Arrange the plots in a 2x2 grid
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)

# Display summary results per individual
print(summary_results_per_individual)




