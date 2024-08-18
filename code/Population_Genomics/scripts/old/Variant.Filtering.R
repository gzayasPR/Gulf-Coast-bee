# Install and load necessary packages if not already installed
if (!requireNamespace("vcfR", quietly = TRUE)) install.packages("vcfR")
if (!requireNamespace("genetics", quietly = TRUE)) install.packages("genetics")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("bedr", quietly = TRUE)) install.packages("bedr")
if (!requireNamespace("gaston", quietly = TRUE)) install.packages("gaston")
library(vcfR)
library(genetics)
library(parallel)
library(bedr)
library(gaston)

rm(list=ls())

# Set project directory and paths
proj.dir <- "C:/Users/gzayas97/Documents/USDA Internship/"
setwd(proj.dir)
data.dir <- paste(proj.dir, "/data/", sep="")
vcf.path <- paste(data.dir, "/Hesperapis_oraria.vcf", sep="")

# Function to calculate MAF for a single locus
calculate_maf <- function(genotypes) {
  ref_count <- 0
  alt_count <- 0
  total_alleles <- 0
  
  for (genotype in genotypes) {
    if (!is.na(genotype)) {
      if (genotype == "0/0" || genotype == "0|0") {
        ref_count <- ref_count + 2
        total_alleles <- total_alleles + 2
      } else if (genotype == "0/1" || genotype == "1/0" || genotype == "0|1" || genotype == "1|0") {
        ref_count <- ref_count + 1
        alt_count <- alt_count + 1
        total_alleles <- total_alleles + 2
      } else if (genotype == "1/1" || genotype == "1|1") {
        alt_count <- alt_count + 2
        total_alleles <- total_alleles + 2
      } else if (genotype == "0") {
        ref_count <- ref_count + 1
        total_alleles <- total_alleles + 1
      } else if (genotype == "1") {
        alt_count <- alt_count + 1
        total_alleles <- total_alleles + 1
      }
    }
  }
  
  ref_freq <- ref_count / total_alleles
  alt_freq <- alt_count / total_alleles
  maf <- min(ref_freq, alt_freq)
  
  return(maf)
}
# Load necessary library
library(vcfR)




# Read the VCF file
vcf <- read.vcfR(vcf.path)
# Initialize a dataframe to save results
results <- data.frame(
  min.dp = numeric(), max.dp = numeric(),
  min.missing.indiv = numeric(), min.missing.locus = numeric(),
  maf_threshold = numeric(), order_of_filtering = character(),
  retained_markers = integer(), retained_individuals = integer(),
  stringsAsFactors = FALSE
)

# Extract genotype and depth matrices from VCF
genotype_matrix <- extract.gt(vcf)
genotype_matrix[genotype_matrix == '.'] <- NA
depth_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
 
parameters <- read.csv("./data/parameters.csv") 

# Function to calculate MAF for a single locus
calculate_maf <- function(genotypes) {
  ref_count <- 0
  alt_count <- 0
  total_alleles <- 0
  
  for (genotype in genotypes) {
    if (!is.na(genotype)) {
      if (genotype == "0/0" || genotype == "0|0") {
        ref_count <- ref_count + 2
        total_alleles <- total_alleles + 2
      } else if (genotype == "0/1" || genotype == "1/0" || genotype == "0|1" || genotype == "1|0") {
        ref_count <- ref_count + 1
        alt_count <- alt_count + 1
        total_alleles <- total_alleles + 2
      } else if (genotype == "1/1" || genotype == "1|1") {
        alt_count <- alt_count + 2
        total_alleles <- total_alleles + 2
      } else if (genotype == "0") {
        ref_count <- ref_count + 1
        total_alleles <- total_alleles + 1
      } else if (genotype == "1") {
        alt_count <- alt_count + 1
        total_alleles <- total_alleles + 1
      }
    }
  }
  
  ref_freq <- ref_count / total_alleles
  alt_freq <- alt_count / total_alleles
  maf <- min(ref_freq, alt_freq)
  
  return(maf)
}
# Loop through each combination of parameters in the dataframe
for (i in 1:nrow(parameters)) {
  param_row <- parameters[i,]
  param_set <- list(min.dp = param_row$min.dp, max.dp = param_row$max.dp,
                    min.missing.indiv = param_row$min.missing.indiv, min.missing.locus = param_row$min.missing.locus,
                    maf_threshold = param_row$maf_threshold, order_of_filtering = param_row$order_of_filtering)
  
  # Check if this parameter combination is already in the results dataframe
  if (nrow(subset(results, min.dp == param_set$min.dp & max.dp == param_set$max.dp &
                  min.missing.indiv == param_set$min.missing.indiv & min.missing.locus == param_set$min.missing.locus &
                  maf_threshold == param_set$maf_threshold & order_of_filtering == param_set$order_of_filtering)) == 0) {
    
    # Extract genotype and depth matrices from VCF
    gt_matrix  <- genotype_matrix 
    dp_matrix <- depth_matrix 
    gt_matrix[dp_matrix < (param_set$min.dp - 1)] <- NA
    gt_matrix[dp_matrix > param_set$max.dp] <- NA

    # Apply order of filtering
    if (param_set$order_of_filtering == "locus_first") {
      # Filter loci with more than min.missing.locus fraction of samples missing
      missing_fraction_locus <- rowMeans(is.na(gt_matrix))
      loci_filter <- which(missing_fraction_locus <= param_set$min.missing.locus)
      # Apply individual and locus filters
      gt_matrix <- gt_matrix[loci_filter,]
      # Filter individuals with more than min.missing.indiv fraction of loci missing
      missing_fraction_indiv <- colMeans(is.na(gt_matrix))
      indiv_filter <- which(missing_fraction_indiv <= param_set$min.missing.indiv)
      gt_matrix <- gt_matrix[,indiv_filter]
    } else {
      # Filter individuals with more than min.missing.indiv fraction of loci missing
      missing_fraction_indiv <- colMeans(is.na(gt_matrix))
      indiv_filter <- which(missing_fraction_indiv <= param_set$min.missing.indiv)
      gt_matrix <- gt_matrix[,indiv_filter]
      # Filter loci with more than min.missing.locus fraction of samples missing
      missing_fraction_locus <- rowMeans(is.na(gt_matrix))
      loci_filter <- which(missing_fraction_locus <= param_set$min.missing.locus)
      gt_matrix <- gt_matrix[loci_filter,]
    }
    
    # Calculate MAF for each locus
    maf_values <- apply(gt_matrix, 1, calculate_maf)
    
    # Filter loci based on MAF
    maf_filter <- which(maf_values >= param_set$maf_threshold)
    gt_matrix <- gt_matrix[ maf_filter ,]
    
    # Save the results
    retained_markers <- nrow(gt_matrix)
    retained_individuals <- ncol(gt_matrix)
    results <- rbind(results, data.frame(min.dp = param_set$min.dp, max.dp = param_set$max.dp, 
                                         min.missing.indiv = param_set$min.missing.indiv, min.missing.locus = param_set$min.missing.locus,
                                         maf_threshold = param_set$maf_threshold, order_of_filtering = param_set$order_of_filtering,
                                         retained_markers = retained_markers, retained_individuals = retained_individuals,
                                         stringsAsFactors = FALSE))
  }
}

# Print the results
print(results)

# Optionally, save the results to a CSV file
results_path <- paste(proj.dir, "/results/filtering_results.csv", sep="")
write.csv(results, file = results_path, row.names = FALSE)
