library(vcfR)
library(genetics)
library(parallel)
library(bedr)
library(gaston)
library(data.table)
library(adegenet)
rm(list=ls())

args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
vcf.path <- args[2]
meta.path <- args[3]


# Read the VCF file
vcf <- read.vcfR(vcf.path)
meta_data <- read.csv(meta.path)
meta_data$ploidy <- ifelse(meta_data$Sex == "female", yes = as.integer(2), no = as.integer(1))

females.id <- meta_data$BioSample[meta_data$ploidy  == 2]

females.pos <- c(1,which(colnames(vcf@gt) %in% females.id))

vcf_females <-  vcf[,females.pos]
rm(vcf)
vcf <- vcf_females

gt_matrix <- extract.gt(vcf, element = "GT")
gt_matrix[gt_matrix == '.'] <- NA
gt_matrix[gt_matrix == './.'] <- NA
na_counts <- colSums(is.na(gt_matrix))
na_counts <- colSums(is.na(gt_matrix))
na_percentage <- na_counts / nrow(gt_matrix)
  
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

 summary_results_per_individual <- data.frame(
    individual = colnames(gt_matrix),
    observed_heterozygosity = ho_per_individual,
    expected_heterozygosity = he_per_individual,
    hom_ref_counts = hom_ref_counts_per_individual,
    het_counts = het_counts_per_individual,
    hom_alt_counts = hom_alt_counts_per_individual,
    na_counts = na_counts,
    na_percentage = na_percentage,
    num_markers = num_markers_retained
  )
setwd(out.dir)
write.csv(summary_results_per_individual,"het.results.csv")