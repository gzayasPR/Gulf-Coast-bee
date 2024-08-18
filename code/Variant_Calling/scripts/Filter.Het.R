library(vcfR)
rm(list=ls())
args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
out.dir <- "/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/data/variants/DeepVariant"
vcf.path <- args[2]
vcf.path <- "/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/data/variants/DeepVariant/Hesperapis_oraria.vcf" 
meta.path <- args[3]
meta.path <-  "/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/data/variants/DeepVariant/Samples.Metadata.csv" 
name <- args[4]
name<-"Hesperapis_oraria"

# Read the VCF file
vcf <- read.vcfR(vcf.path)
# Load metadata
meta_data <- read.csv(meta.path)
# Extract genotype and depth matrices from VCF

meta_data$ploidy <- ifelse(meta_data$Sex == "female", yes = as.integer(2), no = as.integer(1))
males.id <- meta_data$BioSample[meta_data$ploidy  == 1]
males.pos <- c(1,which(colnames(vcf@gt) %in% males.id))
vcf_males <-  vcf[,males.pos]
genotype_matrix <- extract.gt(vcf_males)
genotype_matrix[genotype_matrix == '.'] <- NA
genotype_matrix[genotype_matrix == './.'] <- NA
HET.calls <- apply(genotype_matrix, 1, function(r) any(r %in% c("0/1","1/0","0|1","1|0") )) 
no_het.calls <- which(HET.calls == FALSE)
filtered_vcf <- vcf [ no_het.calls ,]
filtered_vcf_path <-  paste(out.dir,paste("/",name,"_no.het.vcf.gz",sep=""),sep="")
vcfR::write.vcf(filtered_vcf, file = filtered_vcf_path)
filtered_vcf <- read.vcfR(filtered_vcf_path)
