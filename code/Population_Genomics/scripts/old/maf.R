# Set the working directory

proj.dir <- "C:/Users/gzayas97/Documents/USDA Internship/"
setwd(proj.dir)
data.dir <- paste(proj.dir, "/data/", sep="")
data.dir <- paste(proj.dir, "/results/", sep="")
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
min.dp <- 4
# Read the VCF file
vcf.path <- paste(data.dir, "/Hesperapis_oraria.vcf", sep="")
vcf <- read.vcfR(vcf.path)

# Extract genotype and depth matrices from VCF
genotype_matrix <- extract.gt(vcf, as.numeric = FALSE)
depth_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

# Filter genotypes based on minimum depth
genotype_matrix[depth_matrix < min.dp ] <- NA

# Preserve FORMAT column and update the VCF object with NA genotypes
format_column <- vcf@gt[, 1, drop = FALSE]
vcf@gt <- cbind(format_column, genotype_matrix)

# Write the VCF with NA genotypes to a file
vcf_with_na_path <- paste(data.dir ,"/Hesperapis_oraria_with_na.vcf",sep="")
write.vcf(vcf, file = vcf_with_na_path)

# Reload the VCF file with NA genotypes (optional step for clarity)
vcf_with_na <- read.vcfR(vcf_with_na_path)

# Extract the updated genotype matrix
genotype_matrix <- extract.gt(vcf_with_na,  as.numeric = FALSE)

# Filter markers with more than 25% of samples missing that marker
missing_fraction <- rowMeans(is.na(genotype_matrix))
min.missing <- 0.25
head(genotype_matrix)
filtered_genotype_matrix <- genotype_matrix[missing_fraction <= min.missing, ]

nrow(filtered_genotype_matrix)
# Calculate MAF for each locus
maf_values <- apply(filtered_genotype_matrix, 1, calculate_maf)

# Define MAF threshold for filtering
maf_threshold <- 0.0001

# Filter loci based on MAF
filtered_loci <- which(maf_values >= maf_threshold)
final_filtered_genotype_matrix <- filtered_genotype_matrix[filtered_loci, ]
final_filtered_vcf <- vcf_with_na[filtered_loci,]

nrow(filtered_genotype_matrix)
final_filtered_vcf
head(extract.gt(final_filtered_vcf,  as.numeric = FALSE))
# Write the final filtered VCF to a file
filtered_vcf_path <-  paste(data.dir ,"/Hesperapis_oraria_filtered.vcf",sep="")
write.vcf(final_filtered_vcf, file = filtered_vcf_path)

# Save the R objects
rdata_geno_file <- paste(data.dir ,"D:/beenome100/variants/geno.Rdata")
save(final_filtered_genotype_matrix, filtered_depth_matrix, final_filtered_vcf, vcf, file = rdata_geno_file)
