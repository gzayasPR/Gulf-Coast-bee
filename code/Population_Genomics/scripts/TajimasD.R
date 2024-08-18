rm(list=ls())
args<-commandArgs(TRUE)
# Set project directory and paths
out.dir <- args[1]
meta.path <- args[2]
D.file <- args[3]

setwd(out.dir)
taj.all <- read.table(D.file ,header=T)
png("TD.Hist.png")
hist(taj.all$TajimaD,br=20)
dev.off()
