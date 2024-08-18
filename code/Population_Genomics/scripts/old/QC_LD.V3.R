# Install and load necessary packages if not already installed

if (!requireNamespace("vcfR", quietly = TRUE)) install.packages("vcfR",repos="http://cran.r-project.org")
if (!requireNamespace("genetics", quietly = TRUE)) install.packages("genetics",repos="http://cran.r-project.org")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel",repos="http://cran.r-project.org")
if (!requireNamespace("bedr", quietly = TRUE)) install.packages("bedr",repos="http://cran.r-project.org")
if (!requireNamespace("gaston", quietly = TRUE)) install.packages("gaston",repos="http://cran.r-project.org")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table",repos="http://cran.r-project.org")
if (!requireNamespace("adegenet", quietly = TRUE)) install.packages("adegenet")
library(vcfR)
library(genetics)
library(parallel)
library(bedr)
library(gaston)
library(data.table)
library(adegenet)
rm(list=ls())

# Set project directory and paths
proj.dir <- "/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/"
setwd(proj.dir)
data.dir <- paste(proj.dir, "/data/variants/", sep="")

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

# Function to calculate MAC for a single locus
calculate_mac <- function(genotypes) { 
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
  mac <- min(ref_count, alt_count )
  return(mac)
}

min.dp <- 11
min.missing.locus <- 0
min.missing.indiv <- 0
min.MAC <- 1

# Read the VCF file
vcf <- read.vcfR(vcf.path)

# Extract genotype and depth matrices from VCF
genotype_matrix <- extract.gt(vcf)
genotype_matrix[genotype_matrix == '.'] <- NA
genotype_matrix[genotype_matrix == './.'] <- NA
depth_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
info.matrix <- vcf@gt[,-1]

# Convert info.matrix to data.table for faster processing
info_dt <- as.data.table(info.matrix)
remove_gt <- function(x) {
  sub("^[^:]*", "", x)
}

# Apply the function to all columns in data.table
info_dt <- info_dt[, lapply(.SD, remove_gt)]

# Convert back to matrix if needed
info.matrix <- as.matrix(info_dt)

# Filter genotypes based on minimum and maximum depth
gt_matrix <- genotype_matrix
dp_matrix <- depth_matrix
min.dp <- min.dp -1

gt_matrix[dp_matrix < min.dp] <- NA
dp2_matrix <- dp_matrix

dp_rows <- rowMeans(dp2_matrix, na.rm = TRUE)
mean_dp <- mean(dp_rows, na.rm = TRUE)
print(mean_dp)
dp.max.filter <- which(dp_rows <= mean_dp*10)

gt_matrix[-dp.max.filter,] <- NA
dp_matrix[-dp.max.filter,] <- NA

# Combine gt_matrix and info.matrix
gt_matrix2 <- gt_matrix
gt_matrix2[is.na(gt_matrix2)] <- "."

combined_matrix <- paste0(gt_matrix2, "", info.matrix)

combined_matrix <- matrix(combined_matrix, nrow = nrow(gt_matrix2), ncol = ncol(gt_matrix2))

vcf@gt[,-1] <- combined_matrix

missing_fraction_locus <- rowMeans(is.na(gt_matrix))
loci_filter <- which(missing_fraction_locus <= min.missing.locus)
gt_matrix <- gt_matrix[loci_filter,]
# Filter individuals with more than min.missing.indiv fraction of loci missing
missing_fraction_indiv <- colMeans(is.na(gt_matrix))
indiv_filter <- which(missing_fraction_indiv <= min.missing.indiv) 
gt_matrix <- gt_matrix[,indiv_filter]
indiv_filter <- indiv_filter +1
na_vcf <- vcf[loci_filter,c(1,indiv_filter)]

mac_values <- apply(gt_matrix, 1, calculate_mac)
mac_filter <- which(as.numeric(mac_values) >= as.numeric(min.MAC ))

gt_matrix <- gt_matrix[ mac_filter ,]

# Save the results
retained_markers <- nrow( gt_matrix)
retained_individuals <- ncol(gt_matrix)
filtered_vcf <- na_vcf [ mac_filter ,]

na.only_filtered_vcf_path <-  paste(data.dir ,"/Hesperapis_oraria_na_only.filtered.vcf",sep="")
filtered_vcf_path <-  paste(data.dir ,"/Hesperapis_oraria_filtered.vcf",sep="")

vcfR::write.vcf(na_vcf , file = na.only_filtered_vcf_path)
vcfR::write.vcf(filtered_vcf, file = filtered_vcf_path)
filtered_vcf <- read.vcfR(filtered_vcf_path)
filtered_bed <- gaston::read.vcf(filtered_vcf_path ,convert.chr = FALSE)
filtered_bed@snps$dist <- filtered_bed@snps$pos
ld.y <- LD.thin(filtered_bed, 0.4, max.dist =  500e3, 
                dist.unit = c("bases"))
Snps.keep <- ld.y@snps$id

y <- LD( ld.y, lim = c(1, ncol(ld.y)) )

non_pruned_loci <- which(Snps.keep %in%filtered_vcf@fix[,3])
pruned_vcf <- filtered_vcf[non_pruned_loci,]
# Load metadata
meta.path <- paste(data.dir, "/Samples.Metadata.csv", sep="")
meta_data <- read.csv(meta.path)

# Convert filtered VCF to genlight object
genlight_obj <- vcfR2genlight(filtered_vcf)
genlight_obj_pruned <- vcfR2genlight(pruned_vcf)
# Ensure the order of samples matches between genlight object and metadata
id.order <- genlight_obj@ind.names
meta_data <- meta_data[match(id.order, meta_data$BioSample),]

# Assign ploidy based on metadata
meta_data$ploidy <- ifelse(meta_data$Sex == "female", yes = as.integer(2), no = as.integer(1))
genlight_obj@ploidy <- meta_data$ploidy
genlight_obj_pruned@ploidy <- meta_data$ploidy
# Assign populations based on metadata
genlight_obj@pop <- as.factor(meta_data$Location)
print(genlight_obj)
print(genlight_obj_pruned)
pruned.vcf.path <- paste(data.dir, "/Hesperapis_oraria.pruned.vcf.gz", sep="")
vcfR::write.vcf(pruned_vcf,pruned.vcf.path )
# Save the genlight object and filtered VCF to an R data file
rdata.geno.file <- paste(data.dir, "/geno.Rdata", sep="")
save(genlight_obj,genlight_obj_pruned,  filtered_vcf,pruned_vcf, vcf, file = rdata.geno.file)

