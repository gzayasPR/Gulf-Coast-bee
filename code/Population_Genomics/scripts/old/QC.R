proj.dir <- "C:/Users/gzayas97/Documents/USDA Internship/"
setwd(proj.dir)
data.dir <- paste(proj.dir,"/data/",sep="")

#install.packages("adegenet")
#install.packages("vcfR")
library(vcfR)
library(adegenet)
meta.path <- paste(data.dir,"/Samples.Metadata.csv",sep="")
meta_data <- read.csv(meta.path)
vcf.path <- paste(data.dir,"/Hesperapis_oraria.vcf",sep="")
vcf <- read.vcfR(vcf.path )

# Function to calculate MAF for a single locus
calculate_maf <- function(genotypes) {
  # Initialize counts
  ref_count <- 0
  alt_count <- 0
  total_alleles <- 0
  
  # Iterate over each genotype
  for (genotype in genotypes) {
    if (!is.na(genotype)) {
      if (genotype == "0/0") {
        ref_count <- ref_count + 2
        total_alleles <- total_alleles + 2
      } else if (genotype == "0/1") {
        ref_count <- ref_count + 1
        alt_count <- alt_count + 1
        total_alleles <- total_alleles + 2
      } else if (genotype == "1/1") {
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
  
  # Calculate allele frequencies
  ref_freq <- ref_count / total_alleles
  alt_freq <- alt_count / total_alleles
  
  # Minor allele frequency is the smaller of the two frequencies
  maf <- min(ref_freq, alt_freq)
  
  return(maf)
}
# Extract genotype matrix from VCF
genotype_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# Calculate MAF for each locus
maf_values <- apply(genotype_matrix, 1, calculate_maf)

# Define MAF threshold for filtering
maf_threshold <- 0.05

# Filter loci based on MAF
filtered_loci <- which(maf_values >= maf_threshold)
filtered_vcf <- vcf[filtered_loci,]

genlight_obj <- vcfR2genlight(filtered_vcf )
id.order <- genlight_obj@ind.names
meta_data <- meta_data[order(meta_data$BioSample,id.order),]
meta_data$ploidy <- ifelse(meta_data$Sex == "female",yes = as.integer(2),as.integer(1))
genlight_obj@ploidy <- meta_data$ploidy
genlight_obj@pop <- as.factor(meta_data$Location)
rdata.geno.file <- paste(data.dir,"/geno.Rdata",sep="")
save(genlight_obj,filtered_vcf,vcf,file=rdata.geno.file)
