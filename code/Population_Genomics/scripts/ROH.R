
rm(list=ls())

# Set project directory and paths
proj.dir <- "/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/"
setwd(proj.dir)
results.dir <- paste(proj.dir, "/results/Population_Genomics/ROH/DeepVariant/", sep="")
setwd(results.dir)


snpsinruns <- read.table("plink.hom.summary",header=T)
snpsinruns$index <- nrow(snpsinruns)
snpsinruns$PERCENTAGE <- snpsinruns$UNAFF/7

png("manhattan.png")
plot(snpsinruns$index,snpsinruns$PERCENTAGE)
dev.off()
