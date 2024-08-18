library(SNPRelate)

rm(list=ls())

args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
vcf.file <- args[2]
meta.path <- args[3]

vcf.fn <- vcf.file 
gds.file <- paste(out.dir,"/filtered_data.gds",sep="")
gds.fn <- gds.file 

# Convert VCF to GDS
snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")

# Open GDS file
genofile <- snpgdsOpen(gds.fn)

# Calculate inbreeding coefficient
ibd <- snpgdsIBDMoM(genofile)
inbreeding <- snpgdsIndInb(genofile)

# Calculate kinship
kinship <- snpgdsIBDSelection(ibd, kinship = TRUE)

# Save results
inbreeding.file <- paste(out.dir,"/inbreeding_coeff.csv",sep="")
kinship.file <- paste(out.dir,"/kinship_coeff.csv",sep="")
write.csv(inbreeding, file=inbreeding.file )
write.csv(kinship, file=kinship.file)

snpgdsClose(genofile)
