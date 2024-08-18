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
meta.path <- paste(data.dir,"/Samples.Metadata.csv",sep="")
meta.data <- read.csv(meta.path)
females.id <- meta.data$BioSample[meta.data$Sex == "female"]

# Get the column indices of the individuals
individual_indices <- which(colnames(vcf@gt) %in% females.id)

# Subset the VCF object
females.vcf <- vcf[, c(1, individual_indices)]

# Define the depth thresholds
depth_thresholds <- c(3,7,11,15,20,25,30)
#depth_thresholds <- c(7)
# Function to apply filters and perform analysis
analyze_vcf <- function(vcf, depth_threshold, min_called_fraction = 0.25) {
  # Filter by depth
  dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  gt_matrix <- extract.gt(vcf, element = "GT")
  
  # Mask genotypes with depth less than the threshold
  gt_matrix[dp_matrix < depth_threshold] <- NA
  gt_matrix[dp_matrix > 1000] <- NA  
  # Remove variants with less than the minimum called fraction of individuals
  called_fraction <- rowMeans(is.na(gt_matrix))
  valid_variants <- called_fraction >= min_called_fraction
  gt_matrix <- gt_matrix[-valid_variants, ]
  dp_matrix <- dp_matrix[-valid_variants, ]
  
  # Calculate depth statistics
  avg_depth_per_individual <- colMeans(dp_matrix, na.rm = TRUE)
  median_depth_per_individual <- apply(dp_matrix, 2, median, na.rm = TRUE)
  na_counts <- colSums(is.na(dp_matrix))
  na_percentage <- na_counts/nrow(dp_matrix)
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
    hom_ref_counts = hom_ref_counts_per_individual ,
    het_counts = het_counts_per_individual ,
    hom_alt_counts = hom_alt_counts_per_individual,
    na_counts =  na_counts,
    na_per = na_percentage,
    num_markers = num_markers_retained
  )
  
  return(summary_results_per_individual)
}

# List to store results
results_list <- list()

# Iterate over each depth threshold
for (depth_threshold in depth_thresholds) {
  cat("Analyzing with minimum depth:", depth_threshold, "\n")
  
 results <- analyze_vcf(females.vcf, depth_threshold)
  results_list[[paste0("depth_", depth_threshold)]] <- results
  
  # Plot observed heterozygosity vs average depth
  plot1 <- ggplot(results, aes(x = avg_depth, y = observed_heterozygosity)) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    labs(title = paste("Observed Heterozygosity vs Average Depth (Min Depth:", depth_threshold, ")"),
         x = "Average Depth",
         y = "Observed Heterozygosity")
  
  # Plot observed heterozygosity vs median depth
  plot2 <- ggplot(results, aes(x = median_depth, y = observed_heterozygosity)) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    labs(title = paste("Observed Heterozygosity vs Median Depth (Min Depth:", depth_threshold, ")"),
         x = "Median Depth",
         y = "Observed Heterozygosity")
  
  # Plot expected heterozygosity vs average depth
  plot3 <- ggplot(results, aes(x = avg_depth, y = expected_heterozygosity)) +
    geom_point() +
    geom_smooth(method = "lm", col = "red") +
    labs(title = paste("Expected Heterozygosity vs Average Depth (Min Depth:", depth_threshold, ")"),
         x = "Average Depth",
         y = "Expected Heterozygosity")
  
  # Plot expected heterozygosity vs median depth
  plot4 <- ggplot(results, aes(x = median_depth, y = expected_heterozygosity)) +
    geom_point() +
    geom_smooth(method = "lm", col = "red") +
    labs(title = paste("Expected Heterozygosity vs Median Depth (Min Depth:", depth_threshold, ")"),
         x = "Median Depth",
         y = "Expected Heterozygosity")
  
  # Arrange the plots in a 2x2 grid
 grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)
}

# Function to create plots for a given set of results
create_plots <- function(results, depth_threshold) {
  # Plot observed heterozygosity vs average depth
  plot1 <- ggplot(results, aes(x = avg_depth, y = observed_heterozygosity)) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    labs(title = paste("Observed Heterozygosity vs Average Depth (Min Depth:", depth_threshold, ")"),
         x = "Average Depth",
         y = "Observed Heterozygosity")
  
  # Plot observed heterozygosity vs median depth
  plot2 <- ggplot(results, aes(x = median_depth, y = observed_heterozygosity)) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    labs(title = paste("Observed Heterozygosity vs Median Depth (Min Depth:", depth_threshold, ")"),
         x = "Median Depth",
         y = "Observed Heterozygosity")
  
  # Plot expected heterozygosity vs average depth
  plot3 <- ggplot(results, aes(x = avg_depth, y = expected_heterozygosity)) +
    geom_point() +
    geom_smooth(method = "lm", col = "red") +
    labs(title = paste("Expected Heterozygosity vs Average Depth (Min Depth:", depth_threshold, ")"),
         x = "Average Depth",
         y = "Expected Heterozygosity")
  
  # Plot expected heterozygosity vs median depth
  plot4 <- ggplot(results, aes(x = median_depth, y = expected_heterozygosity)) +
    geom_point() +
    geom_smooth(method = "lm", col = "red") +
    labs(title = paste("Expected Heterozygosity vs Median Depth (Min Depth:", depth_threshold, ")"),
         x = "Median Depth",
         y = "Expected Heterozygosity")
  
  return(list(plot1, plot2, plot3, plot4))
}

# Create a list to store all the plots
all_plots <- list()
median_Ohet_plots <- list()
average_Ohet_plots <- list()
median_Ehet_plots <- list()
average_Ehet_plots <- list()

# Iterate over each depth threshold in the results list and create plots
for (depth_threshold in names(results_list)) {
  results <- results_list[[depth_threshold]]
  plots <- create_plots(results, gsub("depth_", "", depth_threshold))
  median_Ohet_plots <- c(median_Ohet_plots, list(plots[[2]]))
  average_Ohet_plots <- c(average_Ohet_plots, list(plots[[1]]))
  median_Ehet_plots <- c(median_Ehet_plots, list(plots[[4]]))
  average_Ehet_plots <- c(average_Ehet_plots, list(plots[[3]]))
  all_plots <- c(all_plots, plots)
}



# Arrange all the plots into a single figure
plot_all <- do.call(grid.arrange, c(all_plots, ncol = 4, nrow = length(depth_thresholds)))

# Arrange each set of plots into separate figures
plot_median_Ohet<- do.call(grid.arrange, c(median_Ohet_plots, ncol = 1, nrow = length(depth_thresholds)))
plot_mean_Ohet <- do.call(grid.arrange, c(average_Ohet_plots, ncol =1, nrow = length(depth_thresholds)))
plot_median_Ehet  <- do.call(grid.arrange, c(median_Ehet_plots, ncol = 1, nrow = length(depth_thresholds)))
plot_mean_Ehet  <-do.call(grid.arrange, c(average_Ehet_plots, ncol = 1, nrow = length(depth_thresholds)))

setwd(results.dir)
dir.create("DepthvsHet",showWarnings = F)
setwd("DepthvsHet")

# Assuming results_list is your list of data frames with names
names(results_list) <- paste("Minimum Depth of",depth_thresholds,sep=" ")

# Create a new workbook
wb <- createWorkbook()

# Loop through the list and add each data frame to a new sheet with names from the list
for (name in names(results_list)) {
  addWorksheet(wb, name)
  writeData(wb, sheet = name, results_list[[name]])
}

# Save the workbook to a file
saveWorkbook(wb, "results_list.xlsx", overwrite = TRUE)


ggsave("median_Ohet.png",plot_median_Ohet)
ggsave("mean_Ohet.png",plot_mean_Ohet)

ggsave("median_Ehet.png",plot_median_Ehet)
ggsave("mean_Ehet.png",plot_mean_Ehet)